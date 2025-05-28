# repository structure:

```
rnaseq-nf/
├── README.md
├── nextflow.config
├── main.nf
├── modules/
│   ├── local/
│   │   ├── index_genome.nf
│   │   ├── split_reference.nf
│   │   ├── lib_process.nf
│   │   ├── gather_counts.nf
│   │   └── dge_analysis.nf
│   └── nf-core/
│       └── (optional nf-core modules)
├── subworkflows/
│   ├── local/
│   │   ├── prepare_references.nf
│   │   └── process_libraries.nf
│   └── nf-core/
├── workflows/
│   └── rnaseq.nf
├── bin/
│   ├── split_annotations.py
│   ├── split_transcripts_5p3p.py
│   ├── concatenate_trinity_transcripts.py
│   ├── calc_mapping_stats.py
│   └── calc_5p3p_bias.py
├── assets/
│   ├── samplesheet_schema.json
│   └── multiqc_config.yml
├── conf/
│   ├── base.config
│   ├── modules.config
│   ├── test.config
│   └── test_full.config
├── docs/
│   ├── usage.md
│   ├── output.md
│   └── images/
├── test_data/
│   ├── samplesheet.csv
│   ├── genome.fa
│   ├── annotations.gtf
│   └── fastq/
├── .github/
│   └── workflows/
│       ├── ci.yml
│       └── linting.yml
├── .gitignore
├── .nf-core.yml
├── CHANGELOG.md
├── CITATIONS.md
└── LICENSE
```
# rnaseq-nf

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**rnaseq-nf** is a bioinformatics pipeline that performs RNA-seq data analysis including:

- Quality control of raw sequencing reads
- Read alignment to reference genome/transcriptome
- Gene expression quantification
- Differential gene expression analysis
- Pathway enrichment analysis

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run rnaseq-nf -profile test,docker --outdir <OUTDIR>