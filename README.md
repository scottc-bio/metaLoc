# metaLoc: Protein localisation prediction workflow

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.10.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction
metaLoc is a reproducible [NextFlow](https://www.nature.com/articles/nbt.3820) workflow for protein localisation prediction from protein amino acid sequences or, in `--meta` mode, from metagenomic nucleotide sequences. 

This workflow accelerates the investigation of large protein datasets by combining multiple publicly-available tools for the identification of biologically-relevant subsets (e.g. predicted secretomes).

Individual tool outputs and merged results tables are produced to allow simple downstream filtering and analysis.

The workflow schematic can be seen below:

<img src="docs/workflow.png" width="70%">

## Core workflow
Protein sequences are the input for the core workflow. For both eukaryotic and prokaryotic sequences the workflow utilises [SignalP 6.0](https://www.nature.com/articles/s41587-021-01156-3) in `euk` mode for eukaryotic sequences and in `other` mode for prokaryotic sequences for signal peptide prediction. Next [DeepLoc 2.1](https://academic.oup.com/nar/article/52/W1/W215/7642068) is utilised for eukaryotic protein localisation prediction, and [DeepLocPro 1.0](https://academic.oup.com/bioinformatics/article/40/12/btae677/7900293) for prokaryotic protein localisation prediction. Finally, transmembrane helices prediction is performed by [Phobius 1.01](https://academic.oup.com/nar/article/35/suppl_2/W429/2920784) or [DeepTMHMM 1.0](https://www.biorxiv.org/content/10.1101/2022.04.08.487609v1). 

*N.B. DeepTMHMM in this workflow relies on sequence upload for processing, which could be time consuming and unsuitable for sensitive datasets.*

## Meta mode
In `--meta` mode, protein sequences can be predicted from the nucleotide sequences. [AUGUSTUS 3.5.0](https://academic.oup.com/bioinformatics/article/24/5/637/202844) utilised the user-selected model to predict genes from eukaryotic contigs. [Prodigal 2.6.3](https://link.springer.com/article/10.1186/1471-2105-11-119) performs ORF prediction from prokaryotic contigs. For sequences of unknown origin or from mixed community sample, [EukRep 0.6.7](https://pmc.ncbi.nlm.nih.gov/articles/PMC5880246/) can be utilised for kingdom classification of contigs into eukaryotic or prokaryotic sequences. Both branches of the pipeline are then performed simultaneously.

# Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 25.10.2
  - Tested in 25.10.2 and 25.10.4
- [Miniconda, Conda, or Mamba](https://docs.conda.io/en/latest/) 

*N.B. This workflow was built and tested on a Linux system with local execution*

# Installation

Clone the repository and move into the workflow directory.

```bash

git clone https://github.com/scottc-bio/metaLoc.git
cd metaLoc

```

# Licensed software setup

[SignalP 6.0](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0i&platform=fast), [DeepLoc 2.1](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.1&packageversion=2.1&platform=All), and [Phobius 1.01](https://software.sbc.su.se/phobius.html) are distributed under separate academic licenses and cannot therefore be distributed in this repository. Users must obtain these tools directly from the official websites and agree to licensing terms.

After downloading the relevant compressed archives for each tool, place them in the `assets/` directory where the Nextflow workflow will automatically configure and use them during execution. No further manual installation is required.

**Expected file names:**

- SignalP 6.0 - signalp-6.0i.fast.tar.gz
- DeepLoc 2.1 - deeploc-2.1.All.tar.gz
- Phobius 1.01 - phobius101_linux.tgz

# Initial setup

Running the test profile will verify the install and prepare almost all necessary environments for the workflow except for the optional DeepTMHMM tool.

```bash

nextflow run . -profile test

```

This will use `assets/test.fasta` which contains two nucleotide sequences:

- The first 10 Kbp of the first contig of the *Aspergillus nidulans* var. *acristatus* genome assembly [GCA_047715555.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_047715555.1/)
- The first 10 Kbp of the first contig of the *Escherichia coli* metagenomic genome assembly [GCA_977857295.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_977857295.1/)

# Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--fasta` | Yes | – | Input FASTA file containing protein or contig sequences |
| `--organism` | Yes* | – | Organism type (`euk` for eukaryotic sequences or `other` for prokaryotic sequences). Required when `--meta` is false. |
| `--meta` | No | `false` | Enables metagenomic mode (`euk`, `other`, or `mixed`). |
| `--augustus_model` | Conditional | – | Augustus species model. Required when `--meta` is `euk` or `mixed`. |
| `--min_contig_len` | Conditional | – | Minimum contig length for EukRep filtering. Required when `--meta` is `mixed`. |
| `--deeptmhmm` | No | `false` | Use DeepTMHMM instead of Phobius for TM prediction (set to `true`). |
| `--signalptar` | No | `assets/signalp-6.0i.fast.tar.gz` | Path to SignalP archive. |
| `--deeploctar` | No | `assets/deeploc-2.1.All.tar.gz` | Path to DeepLoc archive. |
| `--phobiustar` | No | `assets/phobius101_linux.tgz` | Path to Phobius archive. |

## AUGUSTUS models
AUGUSTUS *ab initio* gene prediction requires a gene model to be selected. To view the list of possible models use:

```bash

cat assets/augustus_model_list.txt

```
*N.B. Gene-models are species-specific and strongly affect gene prediction*

## EukRep filtering

EukRep is unreliable on short contigs. Therefore a minimum contig length must be selected which will be applied to input sequences prior to kingdom classification. Values >= 3000 bp are recommended, and values <= 1000 bp are strongly discouraged.

# Outputs

Results will be written to the `results/` directory:

- **Tool specific outputs**
  Named according to input .fasta file and containing the full outputs from each tool used in subdirectories.

- **Final output**
  Named according to input .fasta file with `_final` suffix and containing:
  - `final_merged.tsv` - All columns from all tool outputs
  - `final_concise.tsv` - Key columns from outputs

In `meta` mode, additional outputs from AUGUSTUS, Prodigal, and EukRep are produced.

The `work/` directory contains per-process intermediate files generated by Nextflow.

In `mixed` mode, separate final results are produced for eukaryotic and prokaryotic sequences.

# License

This project is distributed under the MIT License. See [LICENSE](LICENSE) for details.

External tools used by this workflow are distributed under their own licenses and must be obtained separately.
