# MultiOmicLand: An automated pipeline to generate landscape of microbiome-host interactions through integrative analysis of multi-omic data 
MultiOmicLand is an R workflow that performs integrative analysis for microbiome-host multi-omic data, to identify microbiome-metabolite-host interaction links and generate hypotheses for microbial-host interaction. Central to this workflow is the implementation of a sequential mediation analysis, to identify in silico causality between microbiome and host phenotypes, sequentially mediated through metabolome, host transcriptome and host proteome (`microbiome` -> `metabolome` -> `host transcriptome` -> `host proteome` -> `host phenotype`). The overall workflow consists of the following steps:

1. Dimensionality reduction
2. Module-phenotypic association
3. Sequential mediation analysis
4. Multi-omic link identification
5. Multi-omic link prioritization
6. Driver taxa analysis

## System Requirement
MultiOmicLand is built on R 4.1.0. The running of MultiOmicLand depends on the following R packages and their dependencies:
`pacman`, `gtools`, `verification`, `doParallel`, `foreach`, `magrittr`, `tibble`, `data.table`, `dplyr`, `reshape2`, `WGCNA`, `mediation`
These packages should be updated to the latest version before running the workflow to avoid compatibility problem.

## Installation
The required packages can be installed from CRAN by running the following lines in R:
```
   install.packages("pacman")
   install.packages("gtools")
   install.packages("verification")
   install.packages("doParallel")
   install.packages("foreach")
   install.packages("magrittr")
   install.packages("tibble")
   install.packages("data.table")
   install.packages("dplyr")
   install.packages("reshape2")
   install.packages("WGCNA")
   install.packages("mediation")
   install.packages("tidyverse")
``` 
## Data Requirement
The running of MultiOmicLand requires the input of metagenome, metabolome, host transcriptome and proteome (optional) datasets in the form of data matrices (feature by sample), as well as metadata for demographic and clinical data of interest.

1. `Metagenome`: The metagenomic data should be an abundance matrix of KEGG Orthologs (KOs) by samples. This data can be generated through reads-based (i.e. HUMAnN3) or contigs-based (i.e. assembly, gene prediction and KO mapping) approaches.
2. `Metabolome`: The metabolomic data should be an abundance matrix of metabolites (in KEGG or HMDB ID) by samples. This data can be obtained through targeted or non-targeted metabolomic profiling.
3. `Host transcriptome`: The host transcriptome data should be an abundance matrix of genes (in Hugo Gene Symbol) by samples. This data can be obtained by RNA Sequencing followed by reads mapping (Hisat) and gene calling (RSEM, Subread), or by microarray-based approaches.
4. `Host proteome (optional)`: The host proteome data should be an abundance matrix of proteins (in Hugo Gene Symbol or other IDs that can be mapped to Hugo Gene Symbol) by samples. This data can be otained by customized arrays or other metaproteome approaches.

A demo of the above data is provided in the "source.data" folder.

## Workflow Structure and Usage
The overall workflow is presented in the wrapper scripts pipeline.detail.R and pipeline.simpl.R that contained R scrripts needed to run the pipeline with a detailed or simplified description, respectively. 

Under the working directory, there needs to be a "function" directory with all the R function files to be sourced by the wrapper, a "source.data" directory with all the omic datasets, and a "database" directory that contains all the database-related files. The "function" directory can be downloaded as is. The "database" directory can be downloaded as is or customized by the users.

To run MultiOmicLand, implement the package follow the instruction above, and run the scripts step-by-step in pipeline.detail.R or pipeline.simpl.R.

## Note
1. Users need to provide two files, to map the metabolite IDs in their own metabolomic data to the corresponding IDs in MetaCyc and STITCH databases. The metabolite match can be performed by ID conversion (i.e. in MetaboAnalyst) and/or by compound structural search. Two example files (cmpd2metabo.txt, metabo2CIDm.txt) are provided in the "database" folder.

2. It is advised not to use any sample IDs starting with numbers as R will prepend them with a 'X', which could mess up some scripts in the workflow. Rename these sample IDs consistently in all input files at the beginning.
