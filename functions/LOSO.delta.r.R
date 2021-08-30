#' LOSO.delta.r function
#'
#' 
#' This function calculate delta spearman r in the LOSO analysis. 
#' Users provide an original quantification input for MetaG modules without leaving any species out (MetaG.mod.before);   
#' a directory path containing all module quantification files after removing one species each time (MetaG.mod.after.dir);  
#' a MetaB module quantification input (MetaB.mod.input); 
#' pairs of MetaG-MetaB modules on which the impacted by the species are of interest (module.pairs).    
#' The script calculates the spearman r values between the interested MetaG-MetaB module pairs before and after removing one Species from the metagenomic data each time. 
#' A delta spearman r is calculated for each species. The delta spearman r values are z-score standardized across each MetaG-MetaB module pair. 

#' 
#' 
#' @param  module.pairs    a data frame with columns named "MetaG.module" and "MetaB.module"
#' @param  MetaG.mod.before    gct file produced by ssGSEA2
#' @param  MetaG.mod.after.dir    a directory containing gct files produced by ssGSEA2, one for each species removed. File names must be formated as "Species name-combined.gct"
#' @param  MetaB.mod.input  MetaB quantification input as a data frame with rownames being features and colnames being samples, or as a txt file produced by wgcna 
#' @param  log.file   default "LOSO.deltaR.log"
#' @param  error.file   to record error in reading gct files, default "LOSO.deltaR.error.txt"
#' @param  output.dir  default "LOSO.deltaR.out"
#' @param  output.prefix  
#' 
#' @examples
#' LOSO.delta.r(module.pairs =  df, MetaG.mod.before = "metaG_DR/metaG-combined.gct", MetaG.mod.after.dir = "LOSO_metaG_DR", 
#'              MetaB.mod.input = MetaB.Mod.dat )



LOSO.delta.r <- function(
  module.pairs, 
  MetaG.mod.before = "metaG_DR/metaG-combined.gct",
  MetaG.mod.after.dir = "LOSO_metaG_DR",
  MetaB.mod.input = MetaB.Mod.dat,
  log.file = "LOSO.deltaR.log",
  error.file = "LOSO.deltaR.error.txt",
  output.dir = "LOSO.deltaR.out",
  output.prefix = ""
  
){
  
  
  ## ##########################################################################################
  ## 
  ##  libraries
  ## 
  ## ##########################################################################################
  
  
 # if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
 # if (!requireNamespace("cmapR", quietly = TRUE)) BiocManager::install("cmapR")
  
 # library(cmapR)
  library(doParallel)
  library(dplyr)
  library(data.table)
  
  
  ## ##########################################################################################
  ## 
  ## import data 
  ## 
  ## ##########################################################################################
  
  cat(paste(as.character(Sys.time()), '\n'),  file=log.file, append=T)
  cat('Importing data: \n',  file=log.file, append=T)
  
  # MetaG original module data
  input.ds <- MetaG.mod.before
  gct.unique <- NULL
  dataset <- try(parse.gctx(input.ds), silent = T)
  if(class(dataset) != 'try-error' ){
    
    m <- dataset@mat
    gene.names <- dataset@rid
    gene.descs <- dataset@rdesc
    sample.names <- dataset@cid
    sample.descs <- dataset@cdesc
    
  } else {
    
    ## - cmapR functions stop if ids are not unique
    ## - import gct using readLines and make ids unique
    if(length(grep('rid must be unique', dataset) ) > 0) {
      gct.tmp <- readLines(input.ds)
      #first column
      rid <- gct.tmp %>% sub('\t.*','', .)
      #gct version
      ver <- rid[1]
      #data and meta data columns
      meta <- strsplit(gct.tmp[2], '\t') %>% unlist() %>% as.numeric()
      if(ver=='#1.3')
        rid.idx <- (meta[4]+3) : length(rid)
      else
        rid.idx <- 4:length(rid)
      
      #check whether ids are unique
      if(length(rid[rid.idx]) > length(unique(rid[rid.idx]))){
        warning('rids not unique! Making ids unique and exporting new GCT file...\n\n')
        #make unique
        rid[rid.idx] <- make.unique(rid[rid.idx], sep='_')
        #other columns
        rest <- gct.tmp %>% sub('.*?\t','', .)
        rest[1] <- ''
        gct.tmp2 <- paste(rid, rest, sep='\t') 
        gct.tmp2[1] <-  sub('\t.*','',gct.tmp2[1])
        
        #export
        gct.unique <- sub('\\.gct', '_unique.gct', input.ds)
        writeLines(gct.tmp2, con=gct.unique)
        
        #import using cmapR functions
        dataset <- parse.gctx(fname = gct.unique)
        
        ## extract data 
        m <- dataset@mat
        gene.names <- sub('_[0-9]{1,5}$', '', dataset@rid)
        gene.descs <- dataset@rdesc
        sample.names <- dataset@cid
        sample.descs <- dataset@cdesc
      }
      
    } else { #end if 'rid not unique'
      
      ########################################################
      ## display a more detailed error message if the import 
      ## failed due to other reasons than redundant 'rid'
      stop("\n\nError importing GCT file using 'cmapR::parse.gctx()'. Possible reasons:\n\n1) Please check whether you have the latest version of the 'cmapR' installed. Due to submission to Bioconductor the cmap team changed some naming conventions, e.g 'parse.gctx()' has been renamed to 'parse.gctx()'.\n2) The GCT file doesn't seem to be in the correct format! Please see take a look at https://clue.io/connectopedia/gct_format for details about GCT format.
\nError message returned by 'cmapR::parse.gctx()':\n\n", dataset, '\n\n')
    } 
  } #end if try-error
  MetaG.mod.bf <- m %>% t() %>% data.frame()
  #head(MetaG.mod.bf[,1:6])
  
  
  
  # MetaB module data
  if(class(MetaB.mod.input) == "character"){
    if(grepl("\\.gct$", MetaB.mod.input)) m1 <- parse.gctx(MetaB.mod.input)@mat %>% t() %>% data.frame() else m1 <- data.frame(fread(MetaB.mod.input), row.names=1)
  }else{
    m1 = MetaB.mod.input
    feature.abb_df1 <- cbind.data.frame(feature = rownames(m1),
                                        abb = paste("feature",seq(1,nrow(m1),1),sep = ""),
                                        stringsAsFactors = F)
    rownames(m1) <- sapply(rownames(m1), function(x) feature.abb_df1$abb[which(feature.abb_df1$feature == x)])
    m1 <- t(m1) %>% as.data.frame(stringsAsFactors=F)
    colnames(m1) <- sapply(colnames(m1), function(x) feature.abb_df1$feature[which(feature.abb_df1$abb == x)])
  }
  MetaB.mod <- m1
  #head(MetaB.mod[,1:6])
  
  #all(rownames(MetaB.mod) %in% rownames(MetaG.mod.bf))
  
  MetaG.B.bf <- merge(MetaG.mod.bf, MetaB.mod, by=0) 
  MetaG.B.bf <- MetaG.B.bf[complete.cases(MetaG.B.bf),] # remove NA because otherwise spearman r might be NA
  
  ## ##########################################################################################
  ## 
  ## LOSO analysis 
  ## 
  ## ##########################################################################################
  cat('Starting LOSO analysis: \n',  file=log.file, append=T)
  
  AllSpecies = sub("\\-combined\\.gct","",  basename(list.files(MetaG.mod.after.dir,full.names = T, pattern = "combined.gct"))) 
  if(length(AllSpecies) == 0) stop("\n\nError importing GCT file, GCT file names must contain '-combined.gct' as exported by ssGSEA2 function\n\n")
  
  
  # combination of the module pairs and species 
  allCombs <- expand.grid(paste(module.pairs$MetaG.module, module.pairs$MetaB.module, sep = "_"),
                          AllSpecies)
  
  
  
  moduelPair.res <- NULL
  for( i in c(1:nrow(allCombs))){
    mdp = strsplit(as.character(allCombs$Var1[i]), "_", fixed = T)[[1]]
    
    metaG.md = sub("MetaG\\.", "", mdp[1]  )
    metaB.md = sub("MetaB\\.", "",  mdp[2])
    
    r.bf = cor(MetaG.B.bf[,metaG.md], MetaG.B.bf[,metaB.md], method = "spearman")
    
    
    spc = as.character(allCombs$Var2[i])
    
    spc.mod.files = list.files(MetaG.mod.after.dir, full.names = T, pattern = spc)
    spc.mod.file = spc.mod.files[grepl("-combined.gct", spc.mod.files, fixed = T)]
    
    # cat(paste('module pairs: ', paste(mdp,collapse = " | "), " \nSpecies: ",spc,"\n", sep=""), file=log.file, append=T)
    if(i %% 100 ==  0) cat(paste("-----------progress: ",i," out of ", nrow(allCombs), " combinations.", sep = ""), file=log.file, append=T)
    
    m <- NA
    readGctx_success <- F
    if(T){
      input.ds <- spc.mod.file
      gct.unique <- NULL
      dataset <- try(parse.gctx(input.ds), silent = T)
      if(class(dataset) != 'try-error' ){
        
        m <- dataset@mat
        gene.names <- dataset@rid
        gene.descs <- dataset@rdesc
        sample.names <- dataset@cid
        sample.descs <- dataset@cdesc
        readGctx_success <- T
      } else {
        
        ## - cmapR functions stop if ids are not unique
        ## - import gct using readLines and make ids unique
        if(length(grep('rid must be unique', dataset) ) > 0) {
          gct.tmp <- readLines(input.ds)
          #first column
          #rid <- gct.tmp %>% sub('\t.*','', .)
          rid <-  sub('\t.*','', gct.tmp)
          #gct version
          ver <- rid[1]
          #data and meta data columns
          meta <- as.numeric( unlist( strsplit(gct.tmp[2], '\t'))) 
          if(ver=='#1.3'){rid.idx <- (meta[4]+3) : length(rid)} else {rid.idx <- 4:length(rid)}
          
          #check whether ids are unique
          if(length(rid[rid.idx]) > length(unique(rid[rid.idx]))){
            warning('rids not unique! Making ids unique and exporting new GCT file...\n\n')
            #make unique
            rid[rid.idx] <- make.unique(rid[rid.idx], sep='_')
            #other columns
            rest <- sub('.*?\t','', gct.tmp)
            rest[1] <- ''
            gct.tmp2 <- paste(rid, rest, sep='\t') 
            gct.tmp2[1] <-  sub('\t.*','',gct.tmp2[1])
            
            #export
            gct.unique <- sub('\\.gct', '_unique.gct', input.ds)
            writeLines(gct.tmp2, con=gct.unique)
            
            #import using cmapR functions
            dataset <- parse.gctx(fname = gct.unique)
            
            ## extract data 
            m <- dataset@mat
            gene.names <- sub('_[0-9]{1,5}$', '', dataset@rid)
            gene.descs <- dataset@rdesc
            sample.names <- dataset@cid
            sample.descs <- dataset@cdesc
            
            readGctx_success <- T
            
            
          }
          
        } else { #end if 'rid not unique'
          
          ########################################################
          ## display a more detailed error message if the import 
          ## failed due to other reasons than redundant 'rid'
          cat(paste("\n\nError importing GCT file using 'cmapR::parse.gctx()' for species: ", spc,
                    "\nPossible reasons:\n\n1) Please check whether you have the latest version of the 'cmapR' installed. Due to submission to Bioconductor the cmap team changed some naming conventions, e.g 'parse.gctx()' has been renamed to 'parse.gctx()'.\n2) The GCT file doesn't seem to be in the correct format! Please see take a look at https://clue.io/connectopedia/gct_format for details about GCT format.
\nError message returned by 'cmapR::parse.gctx()':\n\n", sep=""),  file = error.file, append=T)
          
          # next
          
          
        } 
      } 
      
      
    }# wrap the script of importing gct file  
    
    if(!readGctx_success ){
      delta.r_vec = c(metaG.md, metaB.md, spc, NA)
      names(delta.r_vec) <- c("MetaG.module","MetaB.module","Species","delta.spearman.r")
    }else{
      
      spc.metaG.mod.aft <- data.frame(t(m))
      MetaG.B.aft <- merge(spc.metaG.mod.aft, MetaB.mod, by=0) 
      MetaG.B.aft <- MetaG.B.aft[complete.cases(MetaG.B.aft),] # remove NA because otherwise spearman r might be NA
      r.aft = cor(MetaG.B.aft[,metaG.md], MetaG.B.aft[,metaB.md], method = "spearman")
      
      delta.r = r.aft - r.bf
      
      delta.r_vec = c(metaG.md, metaB.md, spc, delta.r)
      names(delta.r_vec) <- c("MetaG.module","MetaB.module","Species","delta.spearman.r")
    }
      
    #delta.r_vec
    
    moduelPair.res <- bind_rows(moduelPair.res, delta.r_vec)
  }
  
  
  final_res <- moduelPair.res %>% mutate(delta.spearman.r = as.numeric(delta.spearman.r)) %>%
    group_by(MetaG.module, MetaB.module) %>% mutate(zscore = (delta.spearman.r - mean(delta.spearman.r, na.rm=T))/sd(delta.spearman.r, na.rm=T))
  
  
  cat('Saving results as file: \n',  file=log.file, append=T)
  if(!dir.exists(output.dir)) dir.create(output.dir)
  write.table(final_res, file = paste(output.dir,"/", output.prefix, "_LOSO_delta.spearman.r.txt", sep = ""), sep = "\t", quote = F, row.names = F)
  
  return(final_res)
  
  
}

