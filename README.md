# Characterizing Cellular Identities and Gene Expression Profiles in Axolotl Regeneration

## Description

We propose to characterize the population heterogeneity of axolotl cells across the span of regeneration. We hypothesize that novel cell populations derived from the blastema express species specific genes and autonomous patterning and regeneration, offering one explanation for why regeneration is not conserved. To identify key cells in regeneration, we propose to profile cellular populations at the site of injury in regenerating limbs as compared to healthy tissue using single cell RNA sequencing, focusing on rare subpopulations that were previously not identified. Our aims for this project are as follows:  

**Aim 1**. Characterize global processes at each step of regeneration

**Aim 2**. Resolve cellular subpopulations and identify cellular drivers of regeneration 


## Overview

This repo contains data and code used by Tyler Dao and Nicholas Hutchins for 20.440 Biological Networks Project, Spring 2021.  

**Repo structure:**  
* *Data* - single cell gene expression matrices
* *Deliverables* - written material on the project
* *Script* - contains R scripts for the analysis

At present, the code is sufficient to perform preliminary analysis on the scRNA-seq dataset present. Dimensionality reduction was performed for visualization using UMAP, and the differential expression of 4 known tissue formation regulation genes (TNMD, ASPN, SPARC, HMGN2) were overlaid on the clusters. 


## Data

This project leverages single-cell datasets from uninjured axolotl (10X Genomics scRNA-seq of 2,375 axolotl uninjured adult upper arm connective tissue). Cells are listed in rows and genes in columns of the csv file. 

*Dataset courtesy of Gerber T, Murawala P, Knapp D, et al.  Science. 2018;362(6413).*



## Installation and Usage

The following analysis uses the [Seurat](https://satijalab.org/seurat/index.html) library in **R**. To install:

```bash
# Enter commands in R (or R studio, if installed)
install.packages('Seurat')
library(Seurat)
```

To run this code, download the .R code file from this repo along with the csv file in "Data." Currently, the script is written under the assumption that the dataset and code is separated in folders "Data" and "Script" respectively. Do change the directory if your organization structure differs.


## References
1. Leigh ND, Dunlap GS, Johnson K, et al. Transcriptomic landscape of the blastema niche in regenerating adult axolotl limbs at single-cell resolution. Nat Commun. 2018;9(1):5153.
2. Gerber T, Murawala P, Knapp D, et al. Single-cell analysis uncovers convergence of cell identities during axolotl limb regeneration. Science. 2018;362(6413).
3. Qin T, Fan CM, Wang TZ, et al. Single-cell RNA-seq reveals novel mitochondria-related musculoskeletal cell populations during adult axolotl limb regeneration process. Cell Death Differ. 2021;28(3):1110-1125.
4. Oviedo NJ, Beane WS. Regeneration: The origin of cancer or a possible cure?. Semin Cell Dev Biol. 2009;20(5):557-564. doi:10.1016/j.semcdb.2009.04.005
5. Sibai M, Parlayan C, Tuğlu P, Öztürk G, Demircan T. Integrative Analysis of Axolotl Gene Expression Data from Regenerative and Wound Healing Limb Tissues. Sci Rep. 2019;9(1):20280. Published 2019 Dec 30. doi:10.1038/s41598-019-56829-6
6. Brbić, M., Zitnik, et al. MARS: discovering novel cell types across heterogeneous single-cell experiments. Nature methods. 2020, 17(12), 1200–1206.
