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
# Instead of direct import, use relative import from the package
import sys
from pathlib import Path

# Add the project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Now import the functions
from postprocessing.mob_merger_v3 import (
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
        gff_line = "contig_1\tISEScan\tinsertion_sequence\t950\t2040\t.\t.\t.\tID=contig_1|insertion_sequence-950:2040;mobile_element_type=IS3_with_TIR"
        result = gff_parser(gff_line)
        expected = ["contig_1", "insertion_sequence", "950", "2040", "contig_1|insertion_sequence-950:2040"]
        assert result == expected
    
    def test_gff_parser_with_complex_attributes(self):
        """Test GFF parsing with complex attributes"""
        gff_line = "contig_2\tproMGE\tmobile_genetic_element\t2000\t6500\t4500\t.\t.\tID=MGE_genome1_contig_2:3500-4500;mge=mi:1;genome_type=COR;mge_type=non-nested;size=4500;n_genes=4;mgeR=int_ctndot:1;name=ISLAND"
        result = gff_parser(gff_line)
        expected = ["contig_2", "mobile_genetic_element", "2000", "6500", "MGE_genome1_contig_2:3500-4500"]
        assert result == expected


class TestMomoParser:
    """Test the MAP/momo GFF parser"""
    
    def test_momo_parser_with_sample_data(self):
        """Test momo parser with sample GFF data"""
        sample_gff = """##gff-version 3
contig_1	ISEScan	insertion_sequence	950	2040	.	.	.	ID=contig_1|insertion_sequence-950:2040;mobile_element_type=IS3_with_TIR
contig_1	ISEScan	terminal_inverted_repeat_element	950	970	.	.	.	ID=contig_1|terminal_inverted_repeat_element-950:970;flanking_site=TIR_1
contig_1	ISEScan	terminal_inverted_repeat_element	2020	2040	.	.	.	ID=contig_1|terminal_inverted_repeat_element-2020:2040;flanking_site=TIR_2
contig_2	ICEfinder	integron	3000	14000	.	.	.	ID=contig_2|integron-3000:14000;mobile_element_type=IME
contig_3	geNomad	plasmid	1	5000	.	.	.	ID=contig_3|plasmid-1:5000;mobile_element_type=plasmid
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
                assert contig1_element[2] == '950'  # start
                assert contig1_element[3] == '2040'  # end
                
                # Test subtypes
                assert len(momo_subtypes) > 0
                
            finally:
                os.unlink(tmp_file.name)
    

class TestPromgeParser:
    """Test the proMGE GFF parser"""
    
    def test_promge_parser_basic(self):
        """Test basic proMGE parsing"""
        sample_gff = """##gff-version 3
contig_1	proMGE	mobile_genetic_element	1500	2500	1000	.	.	ID=MGE_genome1_contig_1:1500-2500;mge=is_tn:1;genome_type=ACC;mge_type=non-nested;size=1000;n_genes=1;mgeR=dde_tnp_is66:1;name=ISLAND
contig_1	.	gene	1501	2100	599	+	.	ID=contig_1_2;Parent=MGE_genome1_contig_1:1501-2100;cluster=contig_1_2;size=599;genome_type=ACC
contig_2	proMGE	mobile_genetic_element	2000	6500	4500	+	.	ID=MGE_genome1_contig_2:3500-4500;mge=mi:1;genome_type=COR;mge_type=non-nested;size=4500;n_genes=4;mgeR=int_ctndot:1;name=ISLAND
contig_4	proMGE	mobile_genetic_element	10000	20000	10000	+	.	ID=MGE_genome1_contig_4:10000-20000;mge=phage_like:1;genome_type=ACC;mge_type=non-nested;size=10000;n_genes=8;mgeR=phage_integrase:1;name=ISLAND
"""
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gff') as tmp_file:
            tmp_file.write(sample_gff)
            tmp_file.flush()
            
            try:
                promge_dict, mgeR = promge_parser(tmp_file.name)
                
                # Test structure
                assert 'contig_1' in promge_dict
                assert 'contig_2' in promge_dict
                assert 'contig_4' in promge_dict
                
                # Test element data
                contig1_element = promge_dict['contig_1'][0]
                assert contig1_element[0] == 'contig_1'
                assert 'insertion_sequence' in contig1_element[1]  # Should contain converted type
                assert contig1_element[2] == '1500'
                assert contig1_element[3] == '2500'
                
                # Test mgeR data
                assert len(mgeR) == 3
                
            finally:
                os.unlink(tmp_file.name)


class TestMapper:
    """Test the mapping function that finds overlaps"""
    
    def test_mapper_no_overlap(self):
        """Test mapper with non-overlapping elements"""
        momofy_dict = {
            'contig_3': [('contig_3', 'plasmid', '1', '5000', 'contig_3|plasmid-1:5000')]
        }
        promge_dict = {
            'contig_4': [('contig_4', 'phage', '10000', '20000', 'contig_4:10000-20000')]
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
            'contig_1': [('contig_1', 'insertion_sequence', '950', '2040', 'contig_1|insertion_sequence-950:2040')]
        }
        promge_dict = {
            'contig_1': [('contig_1', 'insertion_sequence', '1500', '2500', 'contig_1:1500-2500')]
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
        metadata_tuple = ('contig_1', 'insertion_sequence', '950', '2040', 'contig_1:950-2040')
        genome = 'test_genome'
        source = 'MAP'
        extra_attributes = ['merged_label=MAP_unique', 'subtype=IS3_with_TIR']
        
        result = to_print(metadata_tuple, genome, source, extra_attributes)
        
        # Check that line is properly formatted
        assert result.startswith('contig_1\tMAP\tinsertion_sequence\t950\t2040')
        assert 'ID=test_genome|contig_1:950-2040' in result
        assert 'merged_label=MAP_unique' in result
        assert result.endswith('\n')
    

class TestMerger:
    """Test the merger function that writes output"""
    
    def test_merger_creates_output_file(self):
        """Test that merger creates the expected output file"""
        # Prepare test data
        momo_unique = ['contig_3:1-5000']
        pro_unique = ['contig_4:10000-20000']
        final_overlapped = []
        permge_metadata = {
            'contig_3:1-5000': ('contig_3', 'plasmid', '1', '5000', 'contig_3|plasmid-1:5000'),
            'contig_4:10000-20000': ('contig_4', 'phage', '10000', '20000', 'contig_4:10000-20000')
        }
        irdr_dict = {}
        momo_subtypes = {'contig_3:1-5000': 'plasmid'}
        mgeR = {'contig_4:10000-20000': 'phage_integrase'}
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
contig_1	ISEScan	insertion_sequence	950	2040	.	.	.	ID=contig_1|insertion_sequence-950:2040;mobile_element_type=IS3_with_TIR
contig_2	ICEfinder	integron	3000	14000	.	.	.	ID=contig_2|integron-3000:14000;mobile_element_type=IME
"""
        
        # Create sample proMGE GFF
        promge_gff = """##gff-version 3

contig_1	proMGE	mobile_genetic_element	1500	2500	1000	.	.	ID=MGE_genome1_contig_1:1500-2500;mge=is_tn:1;genome_type=ACC;mge_type=non-nested;size=1000;n_genes=1;mgeR=dde_tnp_is66:1;name=ISLAND
contig_4	proMGE	mobile_genetic_element	10000	20000	10000	+	.	ID=MGE_genome1_contig_4:10000-20000;mge=phage_like:1;genome_type=ACC;mge_type=non-nested;size=10000;n_genes=8;mgeR=phage_integrase:1;name=ISLAND
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
    return "contig_1\tISEScan\tinsertion_sequence\t950\t2040\t.\t.\t.\tID=contig_1|insertion_sequence-950:2040;mobile_element_type=IS3_with_TIR"


@pytest.fixture
def sample_momofy_dict():
    """Sample momofy dictionary for testing"""
    return {
        'contig_1': [('contig_1', 'insertion_sequence', '950', '2040', 'contig_1|insertion_sequence-950:2040')],
        'contig_2': [('contig_2', 'integron', '3000', '14000', 'contig_2|integron-3000:14000')]
    }


@pytest.fixture
def sample_promge_dict():
    """Sample promge dictionary for testing"""
    return {
        'contig_1': [('contig_1', 'insertion_sequence', '1500', '2500', 'contig_1:1500-2500')],
        'contig_4': [('contig_4', 'phage', '10000', '20000', 'contig_4:10000-20000')]
    }


# Test runner configuration
if __name__ == '__main__':
    pytest.main([__file__, '-v'])
