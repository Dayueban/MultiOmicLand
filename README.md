# MultiOmicLand: An automated pipeline to generate microbiome-host interaction landscape through integrative analysis of multi-omic data
MultiOmicLand is an R workflow that performs integrative analysis for microbiome-host multi-omic data, to identify microbiome-metabolite-host interaction links and generate hypotheses for microbial-host interaction. The overall workflow consists of the following steps:

1. Dimensionality reduction
2. Module-phenotypic association
3. Sequential mediation analysis
4. Multi-omic link identification
5. Multi-omic link prioritization
6. Driver taxa analysis

## System Requirement
MultiOmicLand is built on R 4.1.0. The running of MultiOmicLand depends on the following R packages and their dependencies:
`pacman`, `gtools`, `verification`, `doParallel`, `foreach`, `magrittr`, `tibble`, `data.table`, `dplyr`, `reshape2`, `WGCNA`, `mediation`

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
``` 
## Data Requirement
The running of MultiOmicLand requires the input of metagenome, metabolome and host transcriptome dataset in the form of data matrices (feature by sample), as well as metadata for demographic and clinical data of interest.

1. Metagenome: The metagenomic data should be an abundance matrix of KEGG Orthologs (KOs) by samples. This data can be generated through reads-based (i.e. HUMAnN3) or contigs-based (i.e. assembly, gene prediction and KO mapping) approaches.
2. Metabolome: The metabolomic data should be an abundance matrix of metabolites (in KEGG or HMDB ID) by samples. This data can be obtained through targeted or non-targeted metabolomic profiling.
3. Host transcriptome: The host transcriptome data should be an abundance matrix of genes (in Hugo Gene Symbol) by samples. This data can be obtained by RNA Sequencing followed by reads mapping (Hisat) and gene calling (RSEM, Subread), or by microarray-based approaches.

A small demo of these data is provided in the "source.data" folder with the prefix 'Demo'.

## Workflow Structure
The overall workflow is presented in the wrapper scripts pipeline.detail.R and pipeline.simpl.R that contained R scrripts needed to run the pipeline with a detailed or simplified description, respectively. 

Under the working directory, there needs to be a "function" directory with all the R function files, a "source.data" directory with all the omic data, and a "database" directory that contains all the database-related files. The "function" directory can be downloaded as is. The "database" directory can be downloaded as is or customized by the users.
