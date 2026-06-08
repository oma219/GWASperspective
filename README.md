# :dna: GWASperspective

Analysis code for the manuscript entitled "*Current challenges in GWAS integration and fine-mapping for variant interpretation*". In this manuscript, we document the current obstacles in utilizing GWAS for downstream functional experiments and varying approaches to combat these issues.

This repo contains the scripts for running the analyses present in the manuscript along with the **plotting scripts** for the figures.

## Installation guide
```sh
git clone https://github.com/oma219/GWASperspective.git
cd GWASperspective
```
The typical installation time for the repository is within a minute on normal desktop computer. For our software dependencies (described below), the time is also typically within minutes since they are R packages or Python scripts.

## How to use this code?
This repo consists mostly of scripts that either run third-party tools like `SuSiE`, `coloc`, `AlphaGenome`, etc. or custom scripts in Python/R that analyze publicly available data.

The repo organized by the figures in our paper and each figure has a corresponding wiki page that goes in-depth on how we ran the code in order perform the analysis. The wiki-pages can be [accessed here](https://github.com/oma219/GWASperspective/wiki).

## Software dependencies
Each software was installed/used following the author's recommendations:
```sh
# coding languages
R (v4.4.1)
Python3 (>=3.10)

# statistial genetics tools
susieR (v0.12.35)( https://github.com/stephenslab/susieR)
PolyFun (commit 8b17b4e) (https://github.com/omerwe/polyfun)
Coloc (v5.2.3) (https://github.com/chr1swallace/coloc)
HyPrColoc (v1.0.0) (https://github.com/cnfoley/hyprcoloc)
GWASBrewer (commit 78e7994)  (https://github.com/jean997/GWASBrewer)

# deep-learning methods
AlphaGenome (v0.2.0) (https://github.com/google-deepmind/alphagenome)
Borzoi (commit 5c93582) (https://github.com/calico/borzoi)
Borzoi Prime  (commit 5c93582)(https://github.com/calico/borzoi-paper/tree/main/extensions/prime)
GPN-MSA (hg38 scores downloaded from https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores)
Enformer (v1) (https://github.com/google-deepmind/deepmind-research/tree/master/enformer)
ChromBPNet (v1.0.1) (https://github.com/kundajelab/chrombpnet)
```

## System requirements
Our analysis code has been tested on Linux machines, it has not be tested on Windows/Mac yet. For the deep-learning related sections (Figure 5), GPUs are utilized to speed up computation however the code can be adapted to run on CPUs.