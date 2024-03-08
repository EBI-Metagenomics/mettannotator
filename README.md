[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# mettannotator

<img align="right" width="200" height="180" src="media/temp-logo.png">

- [ Introduction ](#intro)
- [ Workflow and tools](#wf)
- [ Setup ](#sp)
- [ Install and dependencies ](#install)
- [ Usage ](#usage)
- [ Inputs ](#in)
- [ Outputs ](#out)
- [ Tests ](#test)
- [ Citation ](#cite)

<a name="intro"></a>
## Introduction
**mettannotator** is a bioinformatics pipeline that generates an exhaustive annotation of prokaryotic genomes using existing tools. The output is a GFF file that integrates results of all pipeline components. Results of each individual tool are also provided.

<a name="wf"></a>
## Workflow and tools
<img width="720" height="482" src="media/mettannotator-schema-temp.png">
<br />
<br />

The workflow uses the following tools and databases:


| Tool/Database                                                                                    | Version           | Purpose                                                                                                                |
|--------------------------------------------------------------------------------------------------|-------------------|------------------------------------------------------------------------------------------------------------------------|
| [Prokka](https://github.com/tseemann/prokka)                                                     | 1.14.6            | CDS calling and functional annotation                                                                                  |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)                                       | 2.1.11            | Protein annotation (eggNOG, KEGG, COG, CAZy, GO-terms)                                                                 |
| [eggNOG DB](http://eggnog6.embl.de/download/)                                                    | 5.0.2             | Database for eggNOG-mapper                                                                                             |
| [InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/)                               | 5.62-94.0         | Protein annotation (InterPro, Pfam)                                                                                    |
| [UniFIRE](https://gitlab.ebi.ac.uk/uniprot-public/unifire)                                       | 2023.4            | Protein annotation                                                                                                     |
| [CRISPRCasFinder](https://github.com/dcouvin/CRISPRCasFinder)                                    | 4.3.2             | Annotation of CRISPR arrays                                                                                            |
| [AMRFinderPlus](https://github.com/ncbi/amr)                                                     | 3.11.4            | Antimicrobial resistance gene annotation; virulence factors, biocide, heat, acid, and metal resistance gene annotation |
| [AMRFinderPlus DB](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/)              | 3.11 2023-02-23.1 | Database for AMRFinderPlus                                                                                             |
| [SanntiS](https://github.com/Finn-Lab/SanntiS)                                                   | 0.9.3.2           | Biosynthetic gene cluster annotation                                                                                   |
| Infernal                                                                                         | 1.1.4             | ncRNA predictions                                                                                                      |
| tRNAscan-SE                                                                                      | 2.0.9             | tRNA predictions                                                                                                       |
| Rfam                                                                                             | 14.9              | Identification of SSU/LSU rRNA and other ncRNAs                                                                        |
| [VIRify](https://github.com/EBI-Metagenomics/emg-viral-pipeline)                                 | 2.0.0             | Viral sequence annotation                                                                                              |
| [Mobilome annotation pipeline](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline) | 2.0               | Mobilome annotation                                                                                                    |



<a name="usage"></a>
## Usage

First, prepare a assemblies_sheet with your input data that looks as follows:

`assemblies_sheet.csv`:

```csv
prefix,assembly,taxid
BU_ATCC8492VPI0062_NT5002,BU_ATCC8492VPI0062_NT5002.fa,820
...
```

Now, you can run the pipeline using:

```bash
nextflow run ebi-metagenomics/mettannotator \
   -profile <docker/singularity/...> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

ebi-metagenomics/mettannotator was originally written by the Microbiome Informatics Team at [EMBL-EBI](https://www.ebi.ac.uk/about/teams/microbiome-informatics/)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
