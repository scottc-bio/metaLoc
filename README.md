# metaLoc protein localisation workflow

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.10.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction
metaLoc is a reproducible [NextFlow](https://www.nature.com/articles/nbt.3820) workflow for protein localisation prediction utilising publicly available tools from both protein amino acid (aa) sequences or, in "meta" mode, from metagenomic nucleotide (nt) sequences. 

The main aim of this workflow is to rapidly accelerate investigation of large protein datasets by combining multiple tools quickly and easily for users to obtain bioinformatically isolated subsets of proteins with characteristics of interest, e.g. *in silico* secretomes.

For each protein sequence provided, or derived from contig sequences provided, the results of the user-selected tools for protein characterisation are provided individually and also in merged results tables ready for post-workflow filtering.

The workflow schematic can be seen below:

<img src="docs/workflow.png" width="70%">

## Core workflow
The core workflow accepts protein sequence .fasta files as the input. For both eukaryotic and prokaryotic sequences the workflow utilises [SignalP 6.0](https://www.nature.com/articles/s41587-021-01156-3) in "euk" mode for eukaryotic sequences and in "other" mode for prokaryotic sequences for signal peptide prediction. Next [DeepLoc 2.1](https://academic.oup.com/nar/article/52/W1/W215/7642068) is utilised for eukaryotic protein localisation, and [DeepLocPro 1.0](https://academic.oup.com/bioinformatics/article/40/12/btae677/7900293) for prokaryotic protein localisation. Finally, transmembrane helices prediction is performed utilising either [Phobius 1.01](https://academic.oup.com/nar/article/35/suppl_2/W429/2920784) or [DeepTMHMM 1.0](https://www.biorxiv.org/content/10.1101/2022.04.08.487609v1) depending on the user preference. Despite its improved accuraccy, DeepTMHMM implemented here uploads sequences to the DeepTMHMM server for processing which could be time consuming for large datasets and unsuitable for sensitive data.

## Meta mode
For metagenomic nucleotide sequences (e.g. contigs) of known eukaryotic or prokaryotic origin, protein sequences can be predicted from the nucleotide sequences. [AUGUSTUS 3.5.0](https://academic.oup.com/bioinformatics/article/24/5/637/202844) is used alongside the user selected gene model to predict coding sequences from eukaryotic contigs. [Prodigal 2.6.3](https://link.springer.com/article/10.1186/1471-2105-11-119) is used for ORF prediction from prokaryotic contigs. For sequences of unknown origin (e.g. shotgun metagenomic contigs), [EukRep 0.6.7](https://pmc.ncbi.nlm.nih.gov/articles/PMC5880246/) is utilised for taxonomic classification of contigs into eukaryotic or prokaryotic contig sequences. Both branches of the pipeline are then performed simultaneously.  

# Requirements

- Nextflow ≥ 25.10.2 (see installation instructions at [Nextflow](https://www.nextflow.io/docs/latest/install.html))
  - Tested in 25.10.2 and 25.10.4
- Miniconda, Conda, or Mamba (see installation instructions at [Conda](https://docs.conda.io/en/latest/))

*N.B. This workflow was built and tested on a Linux system with local execution*

# Installation

Clone the repository and move into the workflow directory.

```bash

git clone https://github.com/scottc-bio/metaLoc.git
cd metaLoc

```

# Licensed software setup

[SignalP 6.0](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0i&platform=fast), [DeepLoc 2.1](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.1&packageversion=2.1&platform=All), and [Phobius 1.01](https://software.sbc.su.se/phobius.html) distributed under separate academic licenses and cannot therefore be distributed in this repository. Users must obtain these tools directly from the official websites and agree to licensing terms.

After downloading the relevant compressed archives for each tool, place them in the 'assets/' directory where the Nextflow workflow will automatically configure and use them during execution. No further manual installation is required.

*Expected file names*

- SignalP 6.0 - signalp-6.0i.fast.tar.gz
- DeepLoc 2.1 - deeploc-2.1.All.tar.gz
- Phobius 1.01 - phobius101_linux.tgz

*The workflow automatically expects files with matching names to exist in the 'assets/' directory. However, the workflow can be directed to files with differing names using parameters for each compressed archive (see [Parameters](docs/README.md#Parameters)).





# Parameters
