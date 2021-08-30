setwd("C:/Users/Wang Zhang/Desktop/COPD_multiomics/combine_analysis/pipeline/pipeline")

'under the working directory, there needs to be a "function" directory with all the R function files;
a "source.data" directory with all the omic quantification data;
a "database" directory which contains all the database-related files'

source("functions/ssGSEA2.0.R")
source("functions/io.R")
source("functions/utils.R")
source("functions/SupportingFunc.R")
source("functions/wgcna.R")
source("functions/glm.sigModules.R")
source("functions/MediationAnalysis-parallel.R")
source("functions/ko2cmpd.R")
source("functions/ko2metabo.R")
source("functions/MetaG.MetaB.link.R")
source("functions/MetaB.HostT.link.R")
source("functions/HostT.HostP.link.R")
source("functions/LOSO.ko.R")
source("functions/LOSO.delta.r.R")
source("functions/rf_MetaG.MetaB.HostT.links.R")


if(!suppressPackageStartupMessages(require("pacman"))){
  install.packages("pacman")
}
library(pacman)
p_load(gtools)
p_load(verification)
p_load(doParallel)
p_load(foreach)
p_load(magrittr)
library(tibble)
library(data.table)
library(dplyr)
library(reshape2)
library(WGCNA)
library(mediation)
library(tidyverse)
library(ranger)


## #########################################################################################
##
## dimention reduction
##
## #########################################################################################

# dimension reduction for metagenomics data ----------------

# explanation of arguments :  ####

ssGSEA2(input.ds,              ## input data file in gct format, first column (Name) must contain gene symbols
        output.prefix,         ## prefix used for output tables
        gene.set.databases,    ## list of genesets (in gmt format) to evaluate enrichment on 
        sample.norm.type    = c("rank", "log", "log.rank", "none"),  ## sample normalization
        weight= 0,             ## when weight==0, all genes have the same weight; if weight>0,
                               ##    actual values matter, and can change the resulting score
        statistic           = c("area.under.RES", "Kolmogorov-Smirnov"), ## test statistic
        output.score.type   = c("NES", "ES"),
        nperm               = 200,    ## number of random permutations for NES case
        combine.mode        = c("combine.off", "combine.replace", "combine.add"),  ## default combine.off
        ## "combine.off" do not combine *_UP and *_DN versions in a single score.
        ## "combine.replace" combine *_UP and *_DN versions in a single score.
        ## "combine.add" combine *_UP and *_DN versions in a single score and
        ##    add it but keeping the individual *_UP and *_DN versions.
        min.overlap         = 10, 
        correl.type         = c("rank", "z.score", "symm.rank"),  ## correlation type: "rank"(default), "z.score", "symm.rank"
        global.fdr          = FALSE,   ## if TRUE calculate global FDR; else calculate FDR sample-by-sample
        extended.output     = TRUE,    ## if TRUE the result GCT files will contain statistics about gene coverage etc.  
        par                 =T,   ## if TRUE enable multi-core              
        spare.cores         =1,   ## if par=T, spare.cores define number of cores spared 
        export.signat.gct   =T,   ## if TRUE gct files with expression values for each signature will be generated
        param.file          =T,   ## if TRUE parameter file will be generated  
        log.file            ='ssGEA2.log',
        outputDir           ="1_DimReduction")


# example of usage  :  ####
ssGSEA2(input.ds = "source.data/metagenome.gct", 
        gene.set.databases = "source.data/KEGG_modules.gmt",
        output.prefix = "metaG", outputDir = "1_DimReduction",
        min.overlap = 2, weight = 0, statistic = 'area.under.RES', output.score.type = "NES", nperm = 50,
        export.signat.gct = F, 
        par = T, spare.cores = 5)




# dimension reduction for host transcriptomics data and metabolomics data ----------------

# explanation of arguments :  ####
wgcna(input.ds,   ## input data in txt format, where rows are features and columns are different samples.  
      min.ModuleSize = 10,  ## minimum number of features in a module 
      cut.Height     = 0.1, ## ????
      outputDir      = "1_DimReduction", 
      output.prefix  = "hostT"  )


# example of usage  :  ####
wgcna(input.ds = "source.data/transcriptome.txt", output.prefix = "hostT", outputDir = "1_DimReduction")
wgcna(input.ds = "source.data/metabolome.txt", output.prefix = "metaB", outputDir = "1_DimReduction")




## #########################################################################################
##
## Identify significant modules: 
## modules with significant association to the disease status (COPD / heathy)
##
## #########################################################################################


# explanation of arguments :  ####
glm.sigModules(input.ds,  ## module quantification in the format of a data frame object (with rownames being features and colnames being samples) 
                          ## or the path to the module quantification file generated by ssGSEA2 (gct file) or by WGCNA (txt file)
               meta.file, ## meta data prepared in txt file, must contain a column named SampleID and one column named Y. 
                          ## this function aim to find significant association between omics modules and Y. 
               glm.family = "binomial",   ## "Guassian" for continuous and unbounded Y; 
                                          ## "binomial" for binary Y (0/1 or proportions of "successes" and "failures"); 
                                          ##"Gamma" or "inverse.gaussian" for continuous non-negative  Y;  
                                          ## "poisson" when Y is counts and when mean is equal to variance;  
                                          ## "quasipoisson" when Y is counts
               glm.p = 0.1)   ## p-value in glm model to identify significant modules, default 0.1



# example of usage:  ####
# identify MetaG significant modules 
MetaG.sigMods <- glm.sigModules(input.ds = "1_DimReduction/metaG-combined.gct", 
                                meta.file="source.data/meta.txt",
                                glm.family = "binomial",
                                glm.p = 0.5) ## use 0.25 for test data, default is 0.1

# organize MetaB data into a dataframe with rownames being modules and colnames being samples
MetaB.Mod.dat <- fread("1_DimReduction/metaB.module_eigengene.txt",data.table = F) 
#  dplyr::filter(!grepl("#",`#NAME`,fixed=T)) 
#MetaB.Mod.dat[-1] <- sapply(MetaB.Mod.dat[-1], as.numeric)
MetaB.Mod.dat <- MetaB.Mod.dat %>% tibble::column_to_rownames("V1") %>% t() %>% data.frame()

# then identify metaG significant modules 
MetaB.sigMods <- glm.sigModules(input.ds = MetaB.Mod.dat,
                                meta.file="source.data/meta.txt", 
                                glm.family = "binomial",
                                glm.p = 0.5)


# organize HostT data into a dataframe with rownames being modules and colnames being samples
HostT.Mod.dat <- fread("1_DimReduction/hostT.module_eigengene.txt",data.table = F)
#  dplyr::filter(!grepl("#",`#NAME`,fixed=T)) 
#HostT.Mod.dat[-1] <- sapply(HostT.Mod.dat[-1], as.numeric)
HostT.Mod.dat <- HostT.Mod.dat %>% tibble::column_to_rownames("V1") %>% t() %>% data.frame()

# then identify HostT significant modules ------
HostT.sigMods <- glm.sigModules(input.ds = HostT.Mod.dat,
                                meta.file="source.data/meta.txt",
                                glm.family = "binomial",
                                glm.p = 0.5)


# identify HostP significant modules -------
HostP.data <- data.frame(fread("source.data/sputum_cyto.txt"), row.names = 1 )
HostP.sigFeatures <- glm.sigModules(input.ds = HostP.data,
                                    meta.file="source.data/meta.txt",
                                glm.family = "binomial",
                                glm.p = 0.25)





## #######################################################################################################
##
## mediation analysis 
## Mediation analysis takes long for large amount of module pairs - suggest to run on a multi core server
##
## #######################################################################################################


# This function analyze the impact of Treat omic module on Y (clinical factor), mediated by Mediator omic module. This function generates a 
# resulting file in outputDir with name formated as  "Treator_affects_Y_through_Mediator_parallel.txt"
# explanation of arguments :  #### 
MediationAnalysis_parallel(
    Treat.omic,   ##  "MetaG", "MetaB" or "HostT"  to define output file name and column name in output files
    Mediator.omic, ##   "MetaG", "MetaB" or "HostT"  to define output file name and column name in output files
    Treat.omic.input,##   Treat omic module quantity file as a dataframe with rownames being features and colnames being samples; or as path to the file generated from dimention reduction functions
    Mediator.omic.input, ##   Mediator omic module quantity file as a dataframe with rownames being features and colnames being samples; or as path to the file generated from dimention reduction functions
    meta.mediate, ##   metadata file containing clinical variables
    Y,            ##   the clinical variable to which the impact by Treat omic modules and Mediator omic modules
    Treat.omic.sigModules,   ##   Treat omic modules generated from the glm.sigModules function 
    Mediator.omic.sigModules,    ##    Mediator omic modules generated from the glm.sigModules function
    outputDir, ##   output directory, default "mediation.out"
    log.file, ##   log file, default "mediation.log"
    threads ##   number of threads to perform parallel calculation
)


# example of usage: ####
# MetaG modules affects NUE through MetaB modules  -----------

MediationAnalysis_parallel(Treat.omic = "MetaG", Mediator.omic = "MetaB", 
                           Treat.omic.input = "metaG_DR/metaG-combined.gct",
                           Mediator.omic.input = MetaB.Mod.dat,
                           meta.mediate = "source.data/meta.mediation.txt", Y = "NEU",
                           Treat.omic.sigModules = MetaG.sigMods,
                           Mediator.omic.sigModules = MetaB.sigMods,
                           log.file = "mediation.parallel.log",
			   outputDir = "2_Mediation",
                           threads = 25)
MetaG.MetaB.NEU_medres <- fread("2_Mediation/MetaG_affects_NEU_through_MetaB_parallel.txt")


# MetaB modules affects NUE through HostT modules -----------

MediationAnalysis_parallel(Treat.omic = "MetaB", Mediator.omic = "HostT", 
                           Treat.omic.input = MetaB.Mod.dat,
                           Mediator.omic.input = HostT.Mod.dat,
                           meta.mediate = "source.data/meta.mediation.txt", Y = "NEU",
                           Treat.omic.sigModules = MetaB.sigMods,
                           Mediator.omic.sigModules = HostT.sigMods,
                           log.file = "mediation.parallel.log",
			   outputDir = "2_Mediation",
                           threads = 25)
MetaB.HostT.NEU_medres <- fread("2_Mediation/MetaB_affects_NEU_through_HostT_parallel.txt")


# HostT modules affects NUE through HostP features ----------- 

MediationAnalysis_parallel(Treat.omic = "HostT", Mediator.omic = "HostP", 
                           Treat.omic.input = HostT.Mod.dat,
                           Mediator.omic.input = HostP.data,
                           meta.mediate = "source.data/meta.mediation.NEU.txt", Y = "NEU",
                           Treat.omic.sigModules = HostT.sigMods,
                           Mediator.omic.sigModules = HostP.sigFeatures,
			   outputDir = "2_Mediation",
                           threads = 30)
HostT.HostP.NEU_medres <- fread("2_Mediation/HostT_affects_NEU_through_HostP_parallel.txt")




## #######################################################################################################
##
## Biological link analysis  
## 
##
## #######################################################################################################

# MetaG - MetaB links ------------
# 1). ko2cmpd and ko2metabo  
ko2cmpd(dbDir = "database")

#' The ko2cmpd function generates information on substrates, products and substrate&products (in the format of metacyc compound names) for each KO number. 
#' Users should provide the path to a database directory, which must contain the belowing files:
#' ko01000.keg; metacyc_reactions.txt. 
#' This function should be run when the user goes through the pipeline for the first time and whenever the user updated the "ko01000.keg" and "metacyc_reactions.txt" files. 
#' Intermediate and resulting files in the format of RData will be stored in the same database directory.


ko2metabo(dbDir = "database")
#' The ko2metabo function generates information on substrates, products and substrate&products (in the format of metabolomic ids) for each KO number. 
#' Users should provide the path to a database directory, which must contain the cmpd2metabo.txt file and KO2CMPD.lists.RData  file.
#' This function should be run after the ko2cmpd() function
#' Resulting files in the format of RData will be stored in the same database directory.


# 2) metaG-metaB links
# The MetaG.MetaB.link function allows you to identify MetaG - MetaB module pairs with biological links. 
# user should provide a metabo.KEGGmodule.match_file
# explaination of arguments : ####
MetaG.MetaB.link(
  mediation.res,         ##  resulting data frame generated from MediationAnalysis() function
  MetaG_module.feature_file,    ##   KEGG_modules.tab file, source ???
  MetaG_quantity_file,   ##   metagenomic quantification file on the feature level
  KO2METABO_file,        ##   KO2METABO.lists.RData generated from ko2metabo() function
  MetaB_module.feature_file,    ##  module assignment file of the metabolomic featrues, generated by the wgcna function
  MetaB_quantity_file,   ##    metabolomic quantification file on the feature level
  metabo.KEGGmodule.match_file,  ##  metabo.KEGGmodule.match_file
  ACME.p.co,             ##   cutoff of ACME p_value to identify potentially linked MetaG - MetaB modules
  output.dir,            ##   default "biological.links"
  output.prefix,         ##   default ""
  log.file               ##   default "MetaG.MetaB.link.log"
)


# example of usage: ####
MetaG.MetaB.links <- MetaG.MetaB.link( MetaG.MetaB.NEU_medres, 
                                       MetaG_module.feature_file =  "database/KEGG_modules.tab", 
                                       KO2METABO_file = "database/KO2METABO.lists.RData", ## this is an output
                                       MetaB_module.feature_file = "DR_wgcna/metaB.module_assign.txt",
                                       MetaB_quantity_file = "source.data/metabolome.txt",
                                       MetaG_quantity_file = "source.data/metagenome.gct",
                                       metabo.KEGGmodule.match_file = "database/Metabo.KEGGModule.match.txt",
				       output.dir = "3_Biological_Links",
                                       ACME.p.co = 0.25)  





# metaB-hostT links ------------
# The  MetaB.HostT.links  function allows you to identify MetaB - HostT module pairs with biological links.  
# User prepares a  metabo2CIDm file using the perl scripts provided on github.
# Suggest to run the perl scripts on a linux server given the large data size. 
# explaination of arguments : ####

MetaB.HostT.links(
  mediation.res,                     ##   resulting data frame generated from MediationAnalysis() function
  MetaB_quantity_file,               ##   metabolomic quantification file on the feature level
  MetaB_module.feature_file,         ##   module assignment file of the metabolomic featrues, generated by the wgcna function
  HostT_quantity_file,               ##   host transcritomic quantificaiton file on the feature level
  HostT_module.feature_file,         ##   module assignment file of the host transcritomic featrues, generated by the wgcna function
  METABO2CIDm_file,        ##    match file between metabo ID and CIDm ID, generated by the perl scripts provided
  CIDm.receptor_file,      ##    link file between CIDm ID and receptor, available on the githup repo ????
  ACME.p.co,               ##    cutoff of ACME p value to define MetaB-HostT module pairs with potential biological links, default 0.25
  CIDm.receptor.score.co,  ##    default 700
  output.dir,              ##    default "biological.links"
  output.prefix,           ##    default is ""
  log.file                 ##    default  "MetaB.HostT.link.log"
)
  




# example of usage: ####
MetaB.HostT.links <- MetaB.HostT.link(mediation.res = MetaB.HostT.NEU_medres,
                                      MetaB_quantity_file = "source.data/metabolome.txt",  
                                      MetaB_module.feature_file = "DR_wgcna/metaB.module_assign.txt", 
                                      HostT_quantity_file = "source.data/transcriptome.txt", 
                                      HostT_module.feature_file = "DR_wgcna/hostT.module_assign.txt", 
                                      METABO2CIDm_file =  "database/metabo2CIDm.txt",  
                                      CIDm.receptor_file = "database/all_cidm_receptor.txt",  
				      output.dir = "3_Biological_Links",
                                      ACME.p.co = 0.25)  


# hostT-hostP links ------------
# The HostT.HostP.link function allows you to identify biological linkes between HostT modules and HostP features.
# User should provide a protein.gene_file and a Pthway2Gene_file . 
# explaination of arguments : ####

HostT.HostP.link(
  mediation.res,              ##     resulting data frame generated from MediationAnalysis() function
  HostT_quantity_input,       ##     HostT quantification on the feature level as a file path, or as a dataframe with rownames being features and colnames being samples
  HostP_quantity_input,       ##     HostP quantification on the feature level as a file path, or as a dataframe with rownames being features and colnames being samples
  HostT_module.feature_file,  ##     module assignment file of the host transcritomic featrues, generated by the wgcna function
  HostP_protein.gene_file,    ##     protein names and corresponding gene names, sourced from  ????????????????
  Pthway2Gene_file,           ##     Pathway information file (.gmt file) sourced from ????????????
  ACME.p.co,                  ##     cutoff of ACME p value to define MetaB-HostT module pairs with potential biological links, default 0.25
  output.dir,                 ##     default "biological.links"
  output.prefix,              ##     default ""
  log.file                    ##     default "HostT.HostP.link.log"
)
 

# example of usage: ####
HostT.HostP.links <- HostT.HostP.link(mediation.res = HostT.HostP.NEU_medres,
                                      HostT_quantity_input = "source.data/transcriptome.txt",
                                      HostP_quantity_input = HostP.data,
                                      HostT_module.feature_file = "DR_wgcna/hostT.module_assign.txt",
                                      HostP_protein.gene_file = "database/protein_info.txt",
                                      Pthway2Gene_file = "database/pathway.gmt",
				      output.dir = "3_Biological_Links",
                                      ACME.p.co = 0.25)    


### #######################################################################################################
##
##  Random forest
##
## #######################################################################################################

# first identify the MetaG-MetaB-HostT links ---------

MetaG.MetaB.links <- fread("3_Biological_Links/MetaG.MetaB.modules.linked.txt",select = 1, col.names = "V1") %>% unique() %>%
  mutate(MetaG.module = sapply(strsplit(V1, "_", fixed = T), "[[", 1),
         MetaB.module = sapply(strsplit(V1, "_", fixed = T), "[[", 2))


MetaB.HostT.links <- fread("3_Biological_Links/MetaB.HostT.modules.linked.txt",select = 1, col.names = "V1")  %>% 
  unique()  %>% # not unique !! need to check link script
  mutate(MetaB.module = sapply(strsplit(V1, "_", fixed = T), "[[", 1),
         HostT.module = sapply(strsplit(V1, "_", fixed = T), "[[", 2))

MetaG.MetaB.HostT.links <- NULL
for(i in c(1:nrow(MetaG.MetaB.links)) ){
  
  gb.pair = MetaG.MetaB.links$V1[i]
  
  bm = strsplit(gb.pair, "_", fixed = T)[[1]][2]
  
  bt.pairs <- MetaB.HostT.links$V1[which(MetaB.HostT.links$MetaB.module == bm)] 
  
  tmp <- expand.grid(gb.pair, bt.pairs)
  
  MetaG.MetaB.HostT.links <- bind_rows(MetaG.MetaB.HostT.links, tmp)
  
} 

MetaG.MetaB.HostT.links <- MetaG.MetaB.HostT.links %>% mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  mutate(MetaG = sapply( strsplit(Var1,"_",fixed = T), "[[", 1) ,
         MetaB = sapply( strsplit(Var1,"_",fixed = T), "[[", 2) , 
         HostT = sapply( strsplit(Var2,"_",fixed = T), "[[", 2))


# then create a predicted variable data frame ----- 
meta <- fread("source.data/meta.mediation.txt") %>% select(SampleID, NEU) 
#meta <- fread("source.data/meta.txt") %>% select(SampleID, Disease) #%>% mutate(Y = as.factor(as.character(Disease)))

GZ.sp <-meta$SampleID[!grepl("^Z", meta$SampleID)] 


# last, perform random forest analysis -----------
#' The rf_MetaG.MetaB.HostT.Links function allows you to evaluate the accuracy of predicting e.g. a clinical variable by the MetaG.MetaB.HostT links through random forest modeling.  
#' This function generates a performance table recording prediction performance of each MetaG.MetaB.HostT link, as rsme and rsq for "regression", 
#' and as accuracy and roc_auc for "classification". 

# explaination of arguments : ####
rf_MetaG.MetaB.HostT.Links(
  MetaG.Mod.input,          ##    the gct file generated by ssGSEA2 or a data frame with rownames being the modules and colnames being the samples
  MetaB.Mod.input,          ##    the txt file generated by wgcna or a data frame with rownames being the modules and colnames being the samples
  HostT.Mod.input,          ##    the txt file generated by wgcna or a data frame with rownames being the modules and colnames being the samples
  MetaG.MetaB.HostT.link.df,       ##    a data frame with columns named "MetaG", "MetaB, and "HostT" 
  PredictedVar.input,       ##    a data frame with one column named "SampleID" and another column giving the predicted variable (should be factor if the prediction is classification)
  PredictionType,           ##    the type of prediction, "regression" or "classification". 
  Training.samples,         ##    samples used as training set. If defined, all other samples will be used for validation. If not defined, the total samples will be randomly split into training set (3/4) and validation set (1/4) 
  log.file,                 ##    default "rf_by.links.log"
  output.dir,               ##    default "rf_performance"
  output.prefix             ##    default ""
)
  

# example of usage: ####

rf.Performance_cls <- rf_MetaG.MetaB.HostT.Links(MetaG.Mod.input = "1_DimReduction/metaG-combined.gct",
                                                 MetaB.Mod.input = MetaB.Mod.dat,  
                                                 HostT.Mod.input = HostT.Mod.dat,  
                                                 MetaG.MetaB.HostT.link.df = MetaG.MetaB.HostT.links,  
                                                 PredictedVar.input = meta,   
                                                 PredictionType = "regression",  
						 output.dir = "4_RandomForest",
                                                 Training.samples = GZ.sp) 

## #######################################################################################################
##
## LOSO analysis   
##
## #######################################################################################################

# 1) generates KO abundance files for each species excluded


#' The LOSO.ko function allows you to perform LOSO (leave one species out) analysis. 
#' Users should provide 1) gene quantification data, 2) information about Species annotaion of bins, 3) connection between bins and genes through scaffold ids, and 4) KO annotaion of genes. 
#' This function generates one KO abundance file each time after excluding one Species from the gene quantification data.  

# explaination of arguments : ####

LOSO.ko(
  geneDepth_df,     ##  input data frame of gene quantification, with rownames being gene ids and columns being samples; gene ids should be formated as "scaffoldID_number"
  gene.ko_file,     ##    path to file giving informaiton on gene ids and KO number. should contain columns named "Gene" and  "KOnumber". gene ids should be formated in "scaffold_number"
  bin.scaffold_file,  ##  path to file giving information on bins and scaffold id, should contain columns named "Scaffold" and "Bin"
  bin.species_file,   ##   path to file giving information on bins and Species annotaion, should contain columns named "Bin" and "Species"
  log.file            ##  default 'LOSO.ko.log'
)


# example of usage: ####
LOSO.ko(geneDepth_file = "source.data/geneDepth.txt",  
        gene.ko_file = "source.data/ko_noeuk.txt",
        bin.scaffold_file =  "source.data/all_membership.txt",
        bin.species_file = "source.data/bin_species.txt") # memory issue for large geneDepth file, can run on linux server




# 2) perform ssGSEA for each ko abundance file
if(!dir.exists("5_LOSO")) dir.create("5_LOSO")

koFiles <- list.files("5_LOSO/", full.names = T)
for(kof in koFiles){
  specs <-sub("\\.gct$", "", sub("^ko\\.abund_rm\\.", "", basename(kof)) ) 
  
  ssGSEA2(input.ds = kof, 
          gene.set.databases = "source.data/KEGG_modules.gmt",
          output.prefix = specs,
          min.overlap = 2, weight = 0, statistic = 'area.under.RES', output.score.type = "NES",nperm = 100,
          par = T,export.signat.gct = F,
          outputDir = "5_LOSO")
  
}




# 3) calculate delta spearman.r by each species
#' The LOSO.delta.r function calculate delta spearman r in the LOSO analysis. 
#' Users provide an original quantification input for MetaG modules without leaving any species out (MetaG.mod.before);   
#' a directory path containing all module quantification files after removing one species each time (MetaG.mod.after.dir);  
#' a MetaB module quantification input (MetaB.mod.input); 
#' pairs of MetaG-MetaB modules on which the impacted by the species are of interest (module.pairs).    
#' The script calculates the spearman r values between the interested MetaG-MetaB module pairs before and after removing one Species from the metagenomic data each time. 
#' A delta spearman r is calculated for each species. The delta spearman r values are z-score standardized across each MetaG-MetaB module pair. 

# explaination of arguments : ####

LOSO.delta.r(
  module.pairs,           ##  a data frame with columns named "MetaG.module" and "MetaB.module"
  MetaG.mod.before,       ##  gct file produced by ssGSEA2
  MetaG.mod.after.dir,    ##  a directory containing gct files produced by ssGSEA2, one for each species removed. 
                          ##  File names must be formated as "Species name-combined.gct"
  MetaB.mod.input,        ##  MetaB quantification input as a data frame with rownames being modules and colnames being samples, or as a txt file produced by wgcna 
  log.file,               ##  default "LOSO.deltaR.log",
  error.file,             ##  to record error in reading gct files, default "LOSO.deltaR.error.txt"
  output.dir,             ##  default "LOSO.deltaR.out",
  output.prefix,          ##  default  "" 
)
 

# example of usage: ####

# create a test module pair data frame:
test.ModulePairs = fread("biological.links/MetaG.MetaB.modules.linked.txt", select = 1, data.table = F) %>% 
  mutate(MetaG.module = sapply(strsplit(MetaG.MetaB_modulePair,"_", fixed = T),"[[", 1)) %>%
  mutate(MetaB.module = sapply(strsplit(MetaG.MetaB_modulePair,"_", fixed = T),"[[", 2))


LOSO.results <- LOSO.delta.r(module.pairs = test.ModulePairs,
                             MetaG.mod.before = "metaG_DR/metaG-combined.gct",
                             MetaG.mod.after.dir = "LOSO_metaG_DR",
                             MetaB.mod.input = MetaB.Mod.dat,
		             output.dir = "5_LOSO")  

