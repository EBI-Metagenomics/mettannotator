[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# mettannotator

<img align="right" width="162" height="149" src="media/mettannotator-logo.png">

- [ Introduction ](#intro)
- [ Workflow and tools](#wf)
- [ Installation and dependencies ](#install)
  - [Reference databases](#reference-databases)
- [ Usage ](#usage)
- [ Test ](#test)
- [ Outputs ](#out)
- [Preparing annotations for ENA or GenBank submission](#submission)
- [ Mobilome annotation ](#mobilome)
- [ Credits ](#credit)
- [ Contributions and Support ](#contribute)
- [ Citation ](#cite)

<a name="intro"></a>

## Introduction

**mettannotator** is a bioinformatics pipeline that generates an exhaustive annotation of prokaryotic genomes using existing tools. The output is a GFF file that integrates the results of all pipeline components. Results of each individual tool are also provided.

<a name="wf"></a>

## Workflow and tools

<img src="media/mettannotator-schema.png">
<br />
<br />

The workflow uses the following tools and databases:

| Tool/Database                                                                                    | Version                                       | Purpose                                                                                                                |
| ------------------------------------------------------------------------------------------------ | --------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| [Prokka](https://github.com/tseemann/prokka)                                                     | 1.14.6                                        | CDS calling and functional annotation (default)                                                                        |
| [Bakta](https://github.com/oschwengers/bakta)                                                    | 1.9.3                                         | CDS calling and functional annotation (if --bakta flag is used)                                                        |
| [Bakta db](https://zenodo.org/record/10522951/)                                                  | 2024-01-19 with AMRFinderPlus DB 2024-01-31.1 | Bakta DB (when Bakta is used as the gene caller)                                                                       |
| [Pseudofinder](https://github.com/filip-husnik/pseudofinder)                                     | v1.1.0                                        | Identification of possible pseudogenes                                                                                 |
| [Swiss-Prot](https://www.uniprot.org/help/downloads)                                             | 2024_06                                       | Database for Pseudofinder                                                                                              |
| [InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/)                               | 5.62-94.0                                     | Protein annotation (InterPro, Pfam)                                                                                    |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)                                       | 2.1.11                                        | Protein annotation (eggNOG, KEGG, COG, GO-terms)                                                                       |
| [eggNOG DB](http://eggnog6.embl.de/download/)                                                    | 5.0.2                                         | Database for eggNOG-mapper                                                                                             |
| [UniFIRE](https://gitlab.ebi.ac.uk/uniprot-public/unifire)                                       | 2023.4                                        | Protein annotation                                                                                                     |
| [AMRFinderPlus](https://github.com/ncbi/amr)                                                     | 3.12.8                                        | Antimicrobial resistance gene annotation; virulence factors, biocide, heat, acid, and metal resistance gene annotation |
| [AMRFinderPlus DB](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/)              | 3.12 2024-01-31.1                             | Database for AMRFinderPlus                                                                                             |
| [DefenseFinder](https://github.com/mdmparis/defense-finder)                                      | 1.2.0                                         | Annotation of anti-phage systems                                                                                       |
| [DefenseFinder models](https://github.com/mdmparis/defense-finder-models)                        | 1.2.3                                         | Database for DefenseFinder                                                                                             |
| [GECCO](https://github.com/zellerlab/GECCO)                                                      | 0.9.8                                         | Biosynthetic gene cluster annotation                                                                                   |
| [antiSMASH](https://antismash.secondarymetabolites.org/#!/download)                              | 7.1.0                                         | Biosynthetic gene cluster annotation                                                                                   |
| [SanntiS](https://github.com/Finn-Lab/SanntiS)                                                   | 0.9.3.4                                       | Biosynthetic gene cluster annotation                                                                                   |
| [run_dbCAN](https://github.com/linnabrown/run_dbcan)                                             | 4.1.2                                         | PUL prediction                                                                                                         |
| [dbCAN DB](https://bcb.unl.edu/dbCAN2/download/Databases/)                                       | V12                                           | Database for run_dbCAN                                                                                                 |
| [CRISPRCasFinder](https://github.com/dcouvin/CRISPRCasFinder)                                    | 4.3.2                                         | Annotation of CRISPR arrays                                                                                            |
| [cmscan](http://eddylab.org/infernal/)                                                           | 1.1.5                                         | ncRNA predictions                                                                                                      |
| [Rfam](https://rfam.org/)                                                                        | 14.9                                          | Identification of SSU/LSU rRNA and other ncRNAs                                                                        |
| [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE)                                       | 2.0.9                                         | tRNA predictions                                                                                                       |
| [pyCirclize](https://github.com/moshi4/pyCirclize)                                               | 1.4.0                                         | Visualise the merged GFF file                                                                                          |
| [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline)                                 | 2.0.0                                         | Viral sequence annotation (runs separately)                                                                            |
| [Mobilome annotation pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline) | 2.0                                           | Mobilome annotation (runs separately)                                                                                  |

<a name="install"></a>

## Installation and dependencies

This workflow is built using [Nextflow](https://www.nextflow.io/). It uses containers (Docker or Singularity) making installation simple and results highly reproducible.

- Install [Nextflow version >=21.10](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- Install [Singularity](https://github.com/apptainer/singularity/blob/master/INSTALL.md)
- Install [Docker](https://docs.docker.com/engine/install/)

Although it's possible to run the pipeline on a personal computer, due to the compute requirements, we encourage users to run it on HPC clusters. Any HPC scheduler supported by [Nextflow](https://www.nextflow.io/) is compatible; however, our team primarily uses [Slurm](https://slurm.schedmd.com/) and [IBM LSF](https://www.ibm.com/docs/en/spectrum-lsf) for the EBI HPC cluster, so those are the profiles we ship with the pipeline.

<a name="reference-databases"></a>

### Reference databases

The pipeline needs reference databases in order to work, they take roughly 180G.

| Path                | Size |
| ------------------- | ---- |
| amrfinder           | 217M |
| antismash           | 9.4G |
| bakta               | 71G  |
| dbcan               | 7.5G |
| defense_finder      | 242M |
| eggnog              | 48G  |
| interproscan        | 45G  |
| interpro_entry_list | 2.6M |
| rfam_models         | 637M |
| pseudofinder        | 273M |
| total               | 182G |

`mettannotator` has an automated mechanism to download the databases using the `--dbs <db_path>` flag. When this flag is provided, the pipeline inspects the folder to verify if the required databases are already present. If any of the databases are missing, the pipeline will automatically download them.

Users can also provide individual paths to each reference database and its version if needed. For detailed instructions, please refer to the Reference databases section in the `--help` of the pipeline.

It's important to note that users are not allowed to mix the `--dbs` flag with individual database paths and versions; they are mutually exclusive. We recommend users to run the pipeline with the `--dbs` flag for the first time in an appropriate path and to avoid downloading the individual databases separately.

<a name="usage"></a>

## Usage

### Input file

First, prepare an input file in the CSV format that looks as follows:

`assemblies_sheet.csv`:

```csv
prefix,assembly,taxid
BU_ATCC8492VPI0062,/path/to/BU_ATCC8492VPI0062_NT5002.fa,820
EC_ASM584v2,/path/to/GCF_000005845.2.fna,562
...
```

Here,
`prefix` is the prefix and the locus tag that will be assigned to output files and proteins during the annotation process;
maximum length is 24 characters;

`assembly` is the path to where the assembly file in FASTA format is located;

`taxid` is the NCBI TaxId (if the species-level TaxId is not known, a TaxId for a higher taxonomic level can be used). If the taxonomy is known, look up the TaxID [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi).

#### Finding TaxIds

If NCBI taxonomies of input genomes are not known, a tool such as [CAT/BAT](https://github.com/MGXlab/CAT_pack) can be used.
Follow the [instructions](https://github.com/MGXlab/CAT_pack?tab=readme-ov-file#installation) for getting the tool and downloading the NCBI nr database for it.

If using CAT/BAT, here is the suggested process for making the `mettannotator` input file:

```bash
# Run BAT on each input genome, saving all results to the same folder
CAT bins -b ${genome_name}.fna -d ${path_to_CAT_database} -t ${path_to_CAT_tax_folder} -o BAT_results/${genome_name}

# Optional: to check what taxa were assigned, you can add names to them
CAT add_names -i BAT_results/${genome_name}.bin2classification.txt -o BAT_results/${genome_name}.name.txt -t ${path_to_CAT_tax_folder}
```

To generate an input file for `mettannotator`, use [generate_input_file.py](preprocessing/generate_input_file.py):

```
python3 preprocessing/generate_input_file.py -h
usage: generate_input_file.py [-h] -i INFILE -d INPUT_DIR -b BAT_DIR -o OUTFILE [--no-prefix]

The script takes a list of genomes and the taxonomy results generated by BAT and makes a
mettannotator input csv file. The user has the option to either use the genome file name
(minus the extension) as the prefix for mettannotator or leave the prefix off and fill it
out themselves after the script generates an input file with just the FASTA location and
the taxid. It is expected that for all genomes, BAT results are stored in the same folder
and are named as {fasta_base_name}.bin2classification.txt. The script will use the lowest-
level taxid without an asterisk as the taxid for the genome.

optional arguments:
  -h, --help    show this help message and exit
  -i INFILE     A file containing a list of genome files to include (file name only, with file
                extension, unzipped, one file per line).
  -d INPUT_DIR  Full path to the directory where the input FASTA files are located.
  -b BAT_DIR    Folder with BAT results. Results for all genomes should be in the same folder
                and should be named {fasta_base_name}.bin2classification.txt
  -o OUTFILE    Path to the file where the output will be saved to.
  --no-prefix   Skip prefix generation and leave the first column of the output file empty for
                the user to fill out. Default: False
```

For example:

```bash
python3 generate_input_file.py -i list_of_genome_fasta_files.txt -d /path/to/the/fasta/files/folder/ -b BAT_results/ -o mettannotator_input.csv
```

It is always best to check the outputs to ensure the results are as expected. Correct any wrongly detected taxa before starting `mettannotator`.

Note, that by default the script uses FASTA file names as prefixes and truncates them to 24 characters if they exceed the limit.

### Running mettannotator

Running `mettannotator` with the `--help` option will pull the repository and display the help message:

> [!NOTE]
> We use the `-latest` flag with the `nextflow run` command, which ensures that the latest available version of the pipeline is pulled.
> If you encounter any issues with the `nextflow run` command, please refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/reference/cli.html#run).

```angular2html
$ nextflow run -latest ebi-metagenomics/mettannotator/main.nf --help
N E X T F L O W  ~  version 23.04.3
Launching `mettannotator/main.nf` [disturbed_davinci] DSL2 - revision: f2a0e51af6


------------------------------------------------------
  ebi-metagenomics/mettannotator <version>
------------------------------------------------------
Typical pipeline command:

  nextflow run ebi-metagenomics/mettannotator --input assemblies_sheet.csv -profile docker

Input/output options
  --input                            [string]  Path to comma-separated file containing information about the assemblies with the prefix to be used.
  --outdir                           [string]  The output directory where the results will be saved. You have to use absolute paths to storage on Cloud
                                               infrastructure.
  --fast                             [boolean] Run the pipeline in fast mode. In this mode, InterProScan, UniFIRE, and SanntiS won't be executed, saving
                                               resources and speeding up the pipeline.
  --email                            [string]  Email address for completion summary.
  --multiqc_title                    [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.

Reference databases
  --dbs                              [string]  Folder for the tools' reference databases used by the pipeline for downloading. It's important to note that
                                               mixing the --dbs flag with individual database paths and versions is not allowed; they are mutually
                                               exclusive.
  --interproscan_db                  [string]  The InterProScan reference database, ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/
  --interproscan_db_version          [string]  The InterProScan reference database version. [default: 5.62-94.0]
  --interpro_entry_list              [string]  TSV file listing basic InterPro entry information - the accessions, types and names,
                                               ftp://ftp.ebi.ac.uk/pub/databases/interpro/releases/94.0/entry.list
  --interpro_entry_list_version      [string]  InterPro entry list version [default: 94]
  --eggnog_db                        [string]  The EggNOG reference database folder,
                                               https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#requirements
  --eggnog_db_version                [string]  The EggNOG reference database version. [default: 5.0.2]
  --rfam_ncrna_models                [string]  Rfam ncRNA models, ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genomes-pipeline/ncrna/
  --rfam_ncrna_models_rfam_version   [string]  Rfam release version where the models come from. [default: 14.9]
  --amrfinder_plus_db                [string]  AMRFinderPlus reference database,
                                               https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/. Go to the following
                                               documentation for the db setup https://github.com/ncbi/amr/wiki/Upgrading#database-updates.
  --amrfinder_plus_db_version        [string]  The AMRFinderPlus reference database version. [default: 2023-02-23.1]
  --defense_finder_db                [string]  Defense Finder reference models, https://github.com/mdmparis/defense-finder#updating-defensefinder. The
                                               Microbiome Informatics team provides a pre-indexed version of the models for version 1.2.3 on this ftp location:
                                               ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/defense-finder/defense-finder-models_1.2.3.tar.gz.
  --defense_finder_db_version        [string]  The Defense Finder models version. [default: 1.2.3]
  --antismash_db                     [string]  antiSMASH reference database, go to this documentation to do the database setup
                                               https://docs.antismash.secondarymetabolites.org/install/#installing-the-latest-antismash-release.
  --antismash_db_version             [string]  The antiSMASH reference database version. [default: 7.1.0]
  --dbcan_db                         [string]  dbCAN indexed reference database, please go to the documentation for the setup
                                               https://dbcan.readthedocs.io/en/latest/. The Microbiome Informatics team provides a pre-indexed version of the
                                               database for version 4.0 on this ftp location:
                                               ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/dbcan/dbcan_4.0.tar.gz
  --dbcan_db_version                 [string]  The dbCAN reference database version. [default: 4.1.3_V12]
  --pseudofinder_db                  [string]  Pseudofinder reference database. Mettannotator uses SwissProt as the database for Pseudofinder.
  --pseudofinder_db_version          [string]  SwissProt version. [default: 2024_06]

Generic options
  --multiqc_methods_description      [string]  Custom MultiQC yaml file containing HTML including a methods description.

Other parameters
  --bakta                            [boolean] Use Bakta instead of Prokka for CDS annotation. Prokka will still be used for archaeal genomes.

 !! Hiding 17 params, use --validationShowHiddenParams to show them !!
------------------------------------------------------
If you use ebi-metagenomics/mettannotator for your analysis please cite:

* The nf-core framework
  https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
  https://github.com/ebi-metagenomics/mettannotator/blob/master/CITATIONS.md
------------------------------------------------------

```

Now, you can run the pipeline using:

```bash
nextflow run ebi-metagenomics/mettannotator \
   -profile <docker/singularity/...> \
   --input assemblies_sheet.csv \
   --outdir <OUTDIR> \
   --dbs <PATH/TO/WHERE/DBS/WILL/BE/SAVED>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

#### Running the pipeline from the source code

If the Nextflow integration with Git does not work, users can download the tarball from the releases page. After extracting the tarball, the pipeline can be run directly by executing the following command:

```bash
$ nextflow run path-to-source-code/main.nf --help
```

#### Local execution

The pipeline can be run on a desktop or laptop, with the caveat that it will take a few hours to complete depending on the resources. There is a local profile in the Nextflow config that limits the total resources the pipeline can use to 8 cores and 12 GB of RAM. In order to run it (Docker or Singularity are still required):

```bash
nextflow run -latest ebi-metagenomics/mettannotator \
   -profile local,<docker or singulairty> \
   --input assemblies_sheet.csv \
   --outdir <OUTDIR> \
   --dbs <PATH/TO/WHERE/DBS/WILL/BE/SAVED>
```

### Gene caller choice

By default, `mettannotator` uses Prokka to identify protein-coding genes. Users can choose to use Bakta instead by
running `mettannotator` with the `--bakta` flag. `mettannotator` runs Bakta without ncRNA and CRISPR
annotation as these are produced by separate tools in the pipeline. Archaeal genomes will continue to be annotated using
Prokka as Bakta is only intended for annotation of bacterial genomes.

### Fast mode

To reduce the compute time and the amount of resources used, the pipeline can be executed with the `--fast` flag. When
run in the fast mode, `mettannotator` will skip InterProScan, UniFIRE and SanntiS. This could be a suitable option
for a first-pass of annotation or if computational resources are limited, however, we recommend running the full version
of the pipeline whenever possible.

When generating an input file for a fast mode run, it is sufficient to indicate the taxid of the superkingdom (`2` for
bacteria and `2157` for Archaea) in the "taxid" column rather than the taxid of the lowest known taxon.

<a name="test"></a>

## Test

To run the pipeline using a test dataset, execute the following command:

```bash
wget https://raw.githubusercontent.com/EBI-Metagenomics/mettannotator/master/tests/test.csv

nextflow run -latest ebi-metagenomics/mettannotator \
   -profile <docker/singularity/...> \
   --input test.csv \
   --outdir <OUTDIR> \
   --dbs <PATH/TO/WHERE/DBS/WILL/BE/SAVED>
```

<a name="out"></a>

## Outputs

The output folder structure will look as follows:

```
└─<PREFIX>
   ├─antimicrobial_resistance
   │  └─amrfinder_plus
   ├─antiphage_defense
   │  └─defense_finder
   ├─biosynthetic_gene_clusters
   │  ├─antismash
   │  ├─gecco
   │  └─sanntis
   ├─functional_annotation
   │  ├─dbcan
   │  ├─eggnog_mapper
   │  ├─interproscan
   │  ├─merged_gff
   │  ├─prokka
   │  ├─pseudofinder
   │  └─unifire
   ├─mobilome
   │  └─crisprcas_finder
   ├─quast
   │  └─<PREFIX>
   │      ├─basic_stats
   │      └─icarus_viewers
   ├─rnas
   │  ├─ncrna
   │  └─trna
   ├─multiqc
   │  ├─multiqc_data
   │  └─multiqc_plots
   │      ├─pdf
   │      ├─png
   │      └─svg
   ├─pipeline_info
   │  ├─software_versions.yml
   │  ├─execution_report_<timestamp>.txt
   │  ├─execution_report_<timestamp>.html
   │  ├─execution_timeline_<timestamp>.txt
   │  ├─execution_timeline_<timestamp>.html
   │  ├─execution_trace_<timestamp>.txt
   │  ├─execution_trace_<timestamp>.html
   │  └─pipeline_dag_<timestamp>.html

```

### Merged GFF

The two main output files for each genome are located in `<OUTDIR>/<PREFIX>/functional_annotation/merged_gff/`:

- `<PREFIX>_annotations.gff`: annotations produced by all tools merged into a single file

- `<PREFIX>_annotations_with_descriptions.gff`: a version of the GFF file above that includes descriptions of all InterPro terms to make the annotations human-readable. Not generated if `--fast` flag was used.

Both files include the genome sequence in the FASTA format at the bottom of the file.

Additionally, for genomes with no more than 50 annotated contigs, a Circos plot of the `<PREFIX>_annotations.gff` file is generated and included in the same folder. An example of such plot is shown below:

<img src="media/circos-plot-example.png">

#### Data sources

Below is an explanation of how each field in column 3 and 9 of the final GFF file is populated. In most cases, information is taken as is from the reporting tool's output.

| Feature (column 3)    | Attribute Name (column 9)                                               | Reporting Tool  | Description                                                                                                                                                                                                 |
| --------------------- | ----------------------------------------------------------------------- | --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| ncRNA                 | all\*                                                                   | cmscan + Rfam   | ncRNA annotation (excluding tRNA)                                                                                                                                                                           |
| tRNA                  | all\*                                                                   | tRNAscan-SE     | tRNA annotation                                                                                                                                                                                             |
| LeftFLANK, RightFLANK | all\*                                                                   | CRISPRCasFinder | CRISPR array flanking sequence                                                                                                                                                                              |
| CRISPRdr              | all\*                                                                   | CRISPRCasFinder | Direct repeat region of a CRISPR array                                                                                                                                                                      |
| CRISPRspacer          | all\*                                                                   | CRISPRCasFinder | CRISPR spacer                                                                                                                                                                                               |
| CDS                   | `ID`, `eC_number`, `Name`, `Dbxref`, `gene`, `inference`, `locus_tag`   | Prokka/Bakta    | Protein annotation                                                                                                                                                                                          |
| CDS                   | `product`                                                               | mettannotator   | Product assigned as described in [ Determining the product ](#product)                                                                                                                                      |
| CDS                   | `product_source`                                                        | mettannotator   | Tool that reported the product chosen by mettannotator                                                                                                                                                      |
| CDS                   | `eggNOG`                                                                | eggNOG-mapper   | Seed ortholog from eggNOG                                                                                                                                                                                   |
| CDS                   | `cog`                                                                   | eggNOG-mapper   | COG category                                                                                                                                                                                                |
| CDS                   | `kegg`                                                                  | eggNOG-mapper   | KEGG orthology term                                                                                                                                                                                         |
| CDS                   | `Ontology_term`                                                         | eggNOG-mapper   | GO associations                                                                                                                                                                                             |
| CDS                   | `pfam`                                                                  | InterProScan    | Pfam accessions                                                                                                                                                                                             |
| CDS                   | `interpro`                                                              | InterProScan    | InterPro accessions. In `<PREFIX>_annotations_with_descriptions.gff` each accession is followed by its description and entry type: Domain [D], Family [F], Homologous Superfamily [H], Repeat [R], Site [S] |
| CDS                   | `nearest_MiBIG`                                                         | SanntiS         | MiBIG accession of the nearest BGC to the cluster in the MIBIG space                                                                                                                                        |
| CDS                   | `nearest_MiBIG_class`                                                   | SanntiS         | BGC class of nearest_MiBIG                                                                                                                                                                                  |
| CDS                   | `gecco_bgc_type`                                                        | GECCO           | BGC type                                                                                                                                                                                                    |
| CDS                   | `antismash_bgc_function`                                                | antiSMASH       | BGC function                                                                                                                                                                                                |
| CDS                   | `amrfinderplus_gene_symbol`                                             | AMRFinderPlus   | Gene symbol according to AMRFinderPlus                                                                                                                                                                      |
| CDS                   | `amrfinderplus_sequence_name`                                           | AMRFinderPlus   | Product description                                                                                                                                                                                         |
| CDS                   | `amrfinderplus_scope`                                                   | AMRFinderPlus   | AMRFinderPlus database (core or plus)                                                                                                                                                                       |
| CDS                   | `element_type`, `element_subtype`                                       | AMRFinderPlus   | Functional category                                                                                                                                                                                         |
| CDS                   | `drug_class`, `drug_subclass`                                           | AMRFinderPlus   | Class and subclass of drugs that this gene is known to contribute to resistance of                                                                                                                          |
| CDS                   | `dbcan_prot_type`                                                       | run_dbCAN       | Predicted protein function: transporter (TC), transcription factor (TF), signal transduction protein (STP), CAZyme                                                                                          |
| CDS                   | `dbcan_prot_family`                                                     | run_dbCAN       | Predicted protein family                                                                                                                                                                                    |
| CDS                   | `substrate_dbcan-pul`                                                   | run_dbCAN       | Substrate predicted by dbCAN-PUL search                                                                                                                                                                     |
| CDS                   | `substrate_dbcan-sub`                                                   | run_dbCAN       | Substrate predicted by dbCAN-subfam                                                                                                                                                                         |
| CDS                   | `defense_finder_type`, `defense_finder_subtype`                         | DefenseFinder   | Type and subtype of the anti-phage system found                                                                                                                                                             |
| CDS                   | `uf_prot_rec_fullname`, `uf_prot_rec_shortname`, `uf_prot_rec_ecnumber` | UniFIRE         | Protein recommended full name, short name and EC number according to UniFIRE                                                                                                                                |
| CDS                   | `uf_prot_alt_fullname`, `uf_prot_alt_shortname`, `uf_prot_alt_ecnumber` | UniFIRE         | Protein alternative full name, short name and EC number according to UniFIRE                                                                                                                                |
| CDS                   | `uf_chebi`                                                              | UniFIRE         | ChEBI identifiers                                                                                                                                                                                           |
| CDS                   | `uf_ontology_term`                                                      | UniFIRE         | GO associations                                                                                                                                                                                             |
| CDS                   | `uf_keyword`                                                            | UniFIRE         | UniFIRE keywords                                                                                                                                                                                            |
| CDS                   | `uf_gene_name`, `uf_gene_name_synonym`                                  | UniFIRE         | Gene name and gene name synonym according to UniFIRE                                                                                                                                                        |
| CDS                   | `uf_pirsr_cofactor`                                                     | UniFIRE         | Cofactor names from PIRSR                                                                                                                                                                                   |

\*all attributes in column 9 are populated by the tool
<br>
<br>

<a name="product"></a>

#### Determining the product

The following logic is used by `mettannotator` to fill out the `product` field in the 9th column of the GFF:

<img src="media/mettannotator-product.png">

If the pipeline is executed with the `--fast` flag, only the output of eggNOG-mapper is used to determine the product of proteins that were labeled as hypothetical by the gene caller.

#### Detection of pseudogenes and spurious ORFs

`mettannotator` uses several approaches to detect pseudogenes and spurious ORFs:

- If Bakta is used as the initial annotation tool, `mettannotator` will inherit the pseudogene labels assigned by Bakta.
- `mettannotator` runs Pseudofinder and labels genes that Pseudofinder predicts to be pseudogenes by adding `"pseudo=true"` to the 9th column of the final merged GFF file. If there is a disagreement between Pseudofinder and Bakta and one of the tools calls a gene a pseudogene, it will be labeled as a pseudogene.
- AntiFam, which is a part of InterPro, is used to identify potential spurious ORFs. If an ORF has an AntiFam hit, `mettannotator` will remove it from the final merged GFF file. These ORFs will still appear in the raw outputs of Bakta/Prokka and may appear in other tool outputs.

`mettannotator` produces a report file which is located in the `merged_gff` folder and includes a list of CDS with AntiFam hits and pseudogenes. For each pseudogene, the report shows which tool predicted it.

### Contents of the tool output folders

The output folders of each individual tool contain select output files of the third-party tools used by `mettannotator`. For file descriptions, please refer to the tool documentation. For some tools that don't output a GFF, `mettannotator` converts the output into a GFF.

Note: if the pipeline completed without errors but some of the tool-specific output folders are empty, those particular tools did not generate any annotations to output.

<a name="submission"></a>

## Preparing annotations for ENA or GenBank submission

`mettannotator` produces a final annotation file in GFF3 format. To submit the annotations to data archives, it is first necessary to convert the GFF3 file into the required format, using third-party tools available. `mettannotator` outputs a specially formatted GFF3 file, named `<prefix>_submission.gff` to be used with converters.

### ENA

ENA accepts annotations in the EMBL flat-file format.
Please use [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) to perform the conversion; the repository includes detailed instructions. The two files required for conversion are:

- the genome FASTA file
- `<mettannotator_results_folder>/<prefix>/functional_annotation/merged_gff/<prefix>_submission.gff`

Please note that it is necessary to register the project and locus tags in ENA prior to conversion. Follow links in the [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) repository for more details.

### GenBank

To convert annotations for GenBank submission, please use [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/).
Three files are required:

- the genome FASTA file
- `<mettannotator_results_folder>/<prefix>/functional_annotation/merged_gff/<prefix>_submission.gff`
- Submission template file (can be generated [here](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/))

More instructions on running `table2asn` are available via [GenBank](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/).

<a name="mobilome"></a>

## Mobilome annotation

The mobilome annotation workflow is not currently integrated into `mettannotator`. However, the outputs produced by `mettannotator` can be used to run [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline) and the [mobilome annotation pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline) and the outputs of these tools can be integrated back into the GFF file produced by `mettannotator`.

After installing both tools, follow these steps to add the mobilome annotation:

1. Run the [viral annotation pipeline](https://github.com/EBI-Metagenomics/emg-viral-pipeline):

```bash
nextflow run \
    emg-viral-pipeline/virify.nf \
    -profile <profile> \
    --fasta <genome_fasta.fna> \
    --output <prefix>
```

2. Run the [mobilome annotation pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline):

```bash
nextflow run mobilome-annotation-pipeline/main.nf \
    --assembly <genome_fasta.fna> \
    --user_genes true \
    --prot_gff <mettannotator_results_folder/<prefix>/functional_annotation/merged_gff/<prefix>_annotations.gff \
    --virify true # only if the next two VIRify files exist, otherwise skip this line \
    --vir_gff Virify_output_folder/08-final/gff/<prefix>_virify.gff # only if file exists, otherwise skip this line \
    --vir_checkv Virify_output_folder/07-checkv/\*quality_summary.tsv # only if the GFF file above exists, otherwise skip this line \
    --outdir <mobilome_output_folder> \
    --skip_crispr true \
    --skip_amr true \
    -profile <profile>"
```

3. Integrate the output into the `mettannotator` GFF

```bash
# Add mobilome to the merged GFF produced by mettannotator
python3 postprocessing/add_mobilome_to_gff.py \
    -m <mobilome_output_folder>/gff_output_files/mobilome_nogenes.gff \
    -i <mettannotator_results_folder>/<prefix>/functional_annotation/merged_gff/<prefix>_annotations.gff \
    -o <prefix>_annotations_with_mobilome.gff

# Add mobilome to the GFF with descriptions produced by mettannotator
python3 postprocessing/add_mobilome_to_gff.py \
    -m <mobilome_output_folder>/gff_output_files/mobilome_nogenes.gff \
    -i <mettannotator_results_folder>/<prefix>/functional_annotation/merged_gff/<prefix>_annotations_with_descriptions.gff \
    -o <prefix>_annotations_with_descriptions_with_mobilome.gff
```

4. Optional: regenerate the Circos plot with the mobilome track added

```bash
pip install pycirclize
pip install matplotlib

python3 bin/circos_plot.py \
    -i <prefix>_annotations_with_mobilome.gff \
    -o plot.png \
    -p <prefix> \
    --mobilome
```

<a name="credit"></a>

## Credits

ebi-metagenomics/mettannotator was originally written by the Microbiome Informatics Team at [EMBL-EBI](https://www.ebi.ac.uk/about/teams/microbiome-informatics/)

<a name="contribute"></a>

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

<a name="cite"></a>

## Citations

If you use the software, please cite:

Gurbich TA, Beracochea M, De Silva NH, Finn RD. mettannotator: a comprehensive and scalable Nextflow annotation pipeline for prokaryotic assemblies. bioRxiv 2024.07.11.603040; doi: https://doi.org/10.1101/2024.07.11.603040

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
