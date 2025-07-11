#!/usr/bin/env python3
"""
Unit tests for mob_merger_v3.py
Tests the mobile genetic element merger functionality
"""

import pytest
import tempfile
import os
from unittest.mock import patch, mock_open
import sys

# Import the functions from the script
# You'll need to save the script as mob_merger_v3.py in the same directory
from mob_merger_v3 import (
    gff_parser,
    momo_parser,
    promge_parser,
    mapper,
    to_print,
    merger
)


class TestGffParser:
    """Test the GFF line parser function"""
    
    def test_gff_parser_basic(self):
        """Test basic GFF line parsing"""
        gff_line = "contig_1\tMAP\tinsertion_sequence\t1000\t2000\t.\t+\t.\tID=IS_1;mobile_element_type=IS5"
        result = gff_parser(gff_line)
        expected = ["contig_1", "insertion_sequence", "1000", "2000", "IS_1"]
        assert result == expected
    
    def test_gff_parser_with_complex_attributes(self):
        """Test GFF parsing with complex attributes"""
        gff_line = "contig_2\tproMGE\tphage\t5000\t15000\t.\t-\t.\tID=phage_1;mge=phage:0.8;mgeR=tyrosine:0.9"
        result = gff_parser(gff_line)
        expected = ["contig_2", "phage", "5000", "15000", "phage_1"]
        assert result == expected


class TestMomoParser:
    """Test the MAP/momo GFF parser"""
    
    def test_momo_parser_with_sample_data(self):
        """Test momo parser with sample GFF data"""
        sample_gff = """##gff-version 3
contig_1	MAP	insertion_sequence	1000	2000	.	+	.	ID=IS_contig_1_1;mobile_element_type=IS5
contig_1	MAP	terminal_inverted_repeat_element	950	970	.	+	.	ID=TIR_1;flanking_site=IS_1
contig_1	MAP	terminal_inverted_repeat_element	2020	2040	.	+	.	ID=TIR_2;flanking_site=IS_2
contig_2	MAP	integron	3000	4000	.	+	.	ID=INT_contig_2_1;mobile_element_type=class1
"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gff') as tmp_file:
            tmp_file.write(sample_gff)
            tmp_file.flush()
            
            try:
                momofy_dict, irdr_dict, momo_subtypes = momo_parser(tmp_file.name)
                
                # Test momofy_dict structure
                assert 'contig_1' in momofy_dict
                assert 'contig_2' in momofy_dict
                assert len(momofy_dict['contig_1']) == 1
                assert len(momofy_dict['contig_2']) == 1
                
                # Test element data
                contig1_element = momofy_dict['contig_1'][0]
                assert contig1_element[0] == 'contig_1'  # contig
                assert contig1_element[1] == 'insertion_sequence'  # type
                assert contig1_element[2] == '1000'  # start
                assert contig1_element[3] == '2000'  # end
                
                # Test subtypes
                assert len(momo_subtypes) > 0
                
            finally:
                os.unlink(tmp_file.name)
    
    def test_momo_parser_type_conversion(self):
        """Test that certain types are converted correctly"""
        sample_gff = """##gff-version 3
contig_1	MAP	viral_sequence	1000	2000	.	+	.	ID=VS_1;mobile_element_type=lambda
contig_1	MAP	prophage	3000	4000	.	+	.	ID=PP_1;mobile_element_type=mu
"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gff') as tmp_file:
            tmp_file.write(sample_gff)
            tmp_file.flush()
            
            try:
                momofy_dict, _, _ = momo_parser(tmp_file.name)
                
                # Both should be converted to 'phage'
                assert momofy_dict['contig_1'][0][1] == 'phage'
                assert momofy_dict['contig_1'][1][1] == 'phage'
                
            finally:
                os.unlink(tmp_file.name)


class TestPromgeParser:
    """Test the proMGE GFF parser"""
    
    def test_promge_parser_basic(self):
        """Test basic proMGE parsing"""
        sample_gff = """##gff-version 3
contig_1	proMGE	mobile_genetic_element	1000	3000	.	+	.	ID=MGE_contig_1_1;mge=is_tn:0.8,phage:0.6;mgeR=tyrosine:0.9
contig_2	proMGE	mobile_genetic_element	5000	7000	.	+	.	ID=MGE_contig_2_1;mge=ce:0.7;mgeR=serine:0.8
"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gff') as tmp_file:
            tmp_file.write(sample_gff)
            tmp_file.flush()
            
            try:
                promge_dict, mgeR = promge_parser(tmp_file.name)
                
                # Test structure
                assert 'contig_1' in promge_dict
                assert 'contig_2' in promge_dict
                
                # Test element data
                contig1_element = promge_dict['contig_1'][0]
                assert contig1_element[0] == 'contig_1'
                assert 'insertion_sequence' in contig1_element[1]  # Should contain converted type
                assert contig1_element[2] == '1000'
                assert contig1_element[3] == '3000'
                
                # Test mgeR data
                assert len(mgeR) == 2
                
            finally:
                os.unlink(tmp_file.name)


class TestMapper:
    """Test the mapping function that finds overlaps"""
    
    def test_mapper_no_overlap(self):
        """Test mapper with non-overlapping elements"""
        momofy_dict = {
            'contig_1': [('contig_1', 'insertion_sequence', '1000', '2000', 'IS_1')]
        }
        promge_dict = {
            'contig_1': [('contig_1', 'phage', '5000', '6000', 'contig_1_1')]
        }
        
        momo_unique, pro_unique, final_overlapped, permge_metadata = mapper(momofy_dict, promge_dict)
        
        # Should have unique elements, no overlaps
        assert len(momo_unique) == 1
        assert len(pro_unique) == 1
        assert len(final_overlapped) == 0
        
        # Check metadata
        assert len(permge_metadata) == 2
    
    def test_mapper_with_overlap(self):
        """Test mapper with overlapping elements"""
        momofy_dict = {
            'contig_1': [('contig_1', 'insertion_sequence', '1000', '2000', 'IS_1')]
        }
        promge_dict = {
            'contig_1': [('contig_1', 'phage', '1500', '2500', 'contig_1_1')]
        }
        
        momo_unique, pro_unique, final_overlapped, permge_metadata = mapper(momofy_dict, promge_dict)
        
        # Should have overlapping elements
        assert len(momo_unique) == 0
        assert len(pro_unique) == 0
        assert len(final_overlapped) == 1
        assert len(final_overlapped[0]) == 2  # Two overlapping elements


class TestToPrint:
    """Test the GFF line formatting function"""
    
    def test_to_print_basic(self):
        """Test basic GFF line formatting"""
        metadata_tuple = ('contig_1', 'insertion_sequence', '1000', '2000', 'IS_1')
        genome = 'test_genome'
        source = 'MAP'
        extra_attributes = ['merged_label=MAP_unique', 'subtype=IS5']
        
        result = to_print(metadata_tuple, genome, source, extra_attributes)
        
        # Check that line is properly formatted
        assert result.startswith('contig_1\tMAP\tinsertion_sequence\t1000\t2000')
        assert 'ID=test_genome|contig_1:1000-2000' in result
        assert 'merged_label=MAP_unique' in result
        assert result.endswith('\n')
    
    def test_to_print_type_replacement(self):
        """Test that pipe characters in types are replaced with slashes"""
        metadata_tuple = ('contig_1', 'insertion_sequence|phage', '1000', '2000', 'IS_1')
        genome = 'test_genome'
        source = 'promge/MAP'
        extra_attributes = ['merged_label=complete_overlap']
        
        result = to_print(metadata_tuple, genome, source, extra_attributes)
        
        # Pipe should be replaced with slash
        assert 'insertion_sequence/phage' in result


class TestMerger:
    """Test the merger function that writes output"""
    
    def test_merger_creates_output_file(self):
        """Test that merger creates the expected output file"""
        # Prepare test data
        momo_unique = ['contig_1:1000-2000']
        pro_unique = ['contig_1_1']
        final_overlapped = []
        permge_metadata = {
            'contig_1:1000-2000': ('contig_1', 'insertion_sequence', '1000', '2000', 'IS_1'),
            'contig_1_1': ('contig_1', 'phage', '5000', '6000', 'contig_1_1')
        }
        irdr_dict = {}
        momo_subtypes = {'contig_1:1000-2000': 'IS5'}
        mgeR = {'contig_1_1': 'tyrosine'}
        genome_name = 'test_genome'
        
        # Create temporary directory
        with tempfile.TemporaryDirectory() as tmp_dir:
            original_cwd = os.getcwd()
            os.chdir(tmp_dir)
            
            try:
                merger(
                    momo_unique, pro_unique, final_overlapped, permge_metadata,
                    irdr_dict, momo_subtypes, mgeR, genome_name
                )
                
                # Check that output file was created
                output_file = f"{genome_name}_merged.gff"
                assert os.path.exists(output_file)
                
                # Check file content
                with open(output_file, 'r') as f:
                    content = f.read()
                    assert content.startswith('##gff-version 3\n')
                    assert 'MAP_unique' in content
                    assert 'promge_unique' in content
                    
            finally:
                os.chdir(original_cwd)


class TestIntegration:
    """Integration tests for the complete workflow"""
    
    def test_full_workflow_integration(self):
        """Test the complete workflow with sample data"""
        # Create sample MAP GFF
        map_gff = """##gff-version 3
contig_1	MAP	insertion_sequence	1000	2000	.	+	.	ID=IS_contig_1_1;mobile_element_type=IS5
contig_2	MAP	integron	3000	4000	.	+	.	ID=INT_contig_2_1;mobile_element_type=class1
"""
        
        # Create sample proMGE GFF
        promge_gff = """##gff-version 3
contig_1	proMGE	mobile_genetic_element	1500	2500	.	+	.	ID=MGE_contig_1_1;mge=is_tn:0.8;mgeR=tyrosine:0.9
contig_3	proMGE	mobile_genetic_element	5000	7000	.	+	.	ID=MGE_contig_3_1;mge=phage:0.7;mgeR=serine:0.8
"""
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Write sample files
            map_file = os.path.join(tmp_dir, 'map.gff')
            promge_file = os.path.join(tmp_dir, 'promge.gff')
            
            with open(map_file, 'w') as f:
                f.write(map_gff)
            with open(promge_file, 'w') as f:
                f.write(promge_gff)
            
            # Change to temp directory for output
            original_cwd = os.getcwd()
            os.chdir(tmp_dir)
            
            try:
                # Run the workflow
                momofy_dict, irdr_dict, momo_subtypes = momo_parser(map_file)
                promge_dict, mgeR = promge_parser(promge_file)
                momo_unique, pro_unique, final_overlapped, permge_metadata = mapper(
                    momofy_dict, promge_dict
                )
                merger(
                    momo_unique, pro_unique, final_overlapped, permge_metadata,
                    irdr_dict, momo_subtypes, mgeR, 'test_genome'
                )
                
                # Verify output
                assert os.path.exists('test_genome_merged.gff')
                
                with open('test_genome_merged.gff', 'r') as f:
                    content = f.read()
                    lines = content.strip().split('\n')
                    
                    # Should have header + at least some data lines
                    assert len(lines) >= 2
                    assert lines[0] == '##gff-version 3'
                    
                    # Should contain both unique and overlapped elements
                    assert any('MAP_unique' in line for line in lines)
                    assert any('promge_unique' in line for line in lines)
                    
            finally:
                os.chdir(original_cwd)


# Fixtures for common test data
@pytest.fixture
def sample_gff_line():
    """Sample GFF line for testing"""
    return "contig_1\tMAP\tinsertion_sequence\t1000\t2000\t.\t+\t.\tID=IS_1;mobile_element_type=IS5"


@pytest.fixture
def sample_momofy_dict():
    """Sample momofy dictionary for testing"""
    return {
        'contig_1': [('contig_1', 'insertion_sequence', '1000', '2000', 'IS_1')],
        'contig_2': [('contig_2', 'integron', '3000', '4000', 'INT_1')]
    }


@pytest.fixture
def sample_promge_dict():
    """Sample promge dictionary for testing"""
    return {
        'contig_1': [('contig_1', 'phage', '1500', '2500', 'contig_1_1')],
        'contig_3': [('contig_3', 'conjugative_element', '5000', '7000', 'contig_3_1')]
    }


# Test runner configuration
if __name__ == '__main__':
    pytest.main([__file__, '-v'])
