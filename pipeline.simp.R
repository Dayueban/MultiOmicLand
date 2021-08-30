setwd("C:/Users/Wang Zhang/Desktop/COPD_multiomics/combine_analysis/pipeline/pipeline")

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

ssGSEA2(input.ds = "source.data/metagenome.gct", 
        gene.set.databases = "source.data/KEGG_modules.gmt",
        output.prefix = "metaG", outputDir = "1_DimReduction",
        min.overlap = 2, weight = 0, statistic = 'area.under.RES', output.score.type = "NES", nperm = 50,
        export.signat.gct = F, 
        par = T, spare.cores = 5)

wgcna(input.ds = "source.data/transcriptome.txt", output.prefix = "hostT", outputDir = "1_DimReduction")
wgcna(input.ds = "source.data/metabolome.txt", output.prefix = "metaB", outputDir = "1_DimReduction")




## #########################################################################################
##
## Identify significant modules: 
## modules with significant association to the disease status (COPD / heathy)
##
## #########################################################################################

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
ko2metabo(dbDir = "database")

MetaG.MetaB.links <- MetaG.MetaB.link( MetaG.MetaB.NEU_medres, 
                                       MetaG_module.feature_file =  "database/KEGG_modules.tab", 
                                       KO2METABO_file = "database/KO2METABO.lists.RData", ## this is an output
                                       MetaB_module.feature_file = "DR_wgcna/metaB.module_assign.txt",
                                       MetaB_quantity_file = "source.data/metabolome.txt",
                                       MetaG_quantity_file = "source.data/metagenome.gct",
                                       metabo.KEGGmodule.match_file = "database/Metabo.KEGGModule.match.txt",
				       output.dir = "3_Biological_Links",
                                       ACME.p.co = 0.25)  

MetaB.HostT.links <- MetaB.HostT.link(mediation.res = MetaB.HostT.NEU_medres,
                                      MetaB_quantity_file = "source.data/metabolome.txt",  
                                      MetaB_module.feature_file = "DR_wgcna/metaB.module_assign.txt", 
                                      HostT_quantity_file = "source.data/transcriptome.txt", 
                                      HostT_module.feature_file = "DR_wgcna/hostT.module_assign.txt", 
                                      METABO2CIDm_file =  "database/metabo2CIDm.txt",  
                                      CIDm.receptor_file = "database/all_cidm_receptor.txt",  
				      output.dir = "3_Biological_Links",
                                      ACME.p.co = 0.25)  

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

LOSO.ko(geneDepth_file = "source.data/geneDepth.txt",  
        gene.ko_file = "source.data/ko_noeuk.txt",
        bin.scaffold_file = "source.data/all_membership.txt",
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

# create a test module pair data frame:
test.ModulePairs = fread("biological.links/MetaG.MetaB.modules.linked.txt", select = 1, data.table = F) %>% 
  mutate(MetaG.module = sapply(strsplit(MetaG.MetaB_modulePair,"_", fixed = T),"[[", 1)) %>%
  mutate(MetaB.module = sapply(strsplit(MetaG.MetaB_modulePair,"_", fixed = T),"[[", 2))


LOSO.results <- LOSO.delta.r(module.pairs = test.ModulePairs,
                             MetaG.mod.before = "1_DimReduction/metaG-combined.gct",
                             MetaG.mod.after.dir = "5_LOSO",
                             MetaB.mod.input = MetaB.Mod.dat,
			     output.dir = "5_LOSO")  

