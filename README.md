[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**ebi-metagenomics/mettannotator** is a bioinformatics pipeline use to annotate assemblies using the following tools:

| Tool/Database      | Version | Purpose |
| ----------- | ----------- |----------- |
| Prokka   | 1.14.6        | Protein annotation       |
| eggNOG-mapper  | 2.1.11        | Protein annotation (eggNOG, KEGG, COG,  CAZy)       |
| eggNOG DB  | 5.0       | Database for eggNOG-mapper       |
| Diamond    | 2.0.11       | Protein annotation (eggNOG)       |
| InterProScan   | 5.62-94.0      | Protein annotation (InterPro, Pfam)       |
| CRISPRCasFinder   | 4.3.2        | Annotation of CRISPR arrays       |
| AMRFinderPlus   | 3.11.4        |   Antimicrobial resistance gene annotation; virulence factors, biocide, heat, acid, and metal resistance gene annotation     |
| AMRFinderPlus DB   | 3.11 2023-02-23.1        | Database for AMRFinderPlus      |
| SanntiS   | 0.9.3.2        | Biosynthetic gene cluster annotation       |
| Infernal   | 1.1.4        | RNA predictions       |
| tRNAscan-SE   | 2.0.9       | tRNA predictions       |
| Rfam   | 14.9        | Identification of SSU/LSU rRNA and other ncRNAs       |
| VIRify   | -        | Viral sequence annotation       |
| MoMofy   | 1.0.0        | Mobilome annotation       |

## Usage

First, prepare a assemblies_sheet with your input data that looks as follows:

`assemblies_sheet.csv`:

```csv
prefix,assembly
BU_ATCC8492VPI0062_NT5002,BU_ATCC8492VPI0062_NT5002.fa
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

ebi-metagenomics/mettannotator was originally written by @mberacochea.

We thank the following people for their extensive assistance in the development of this pipeline:
- @KateSakharova
- @tgurbich

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
