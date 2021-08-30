#' MediationAnalysis function
#'
#' This function analyze the impact of Treat omic module on Y (clinical factor), mediated by Mediator omic module
#' Dependent packages includes "mediation","data.table" and "dplyr" 
#' Calculations takes long time for large module amount.  
#' @param  Treat.omic    "MetaG", "MetaB" or "HostT"  to define output file name and column name in output files
#' @param  Mediator.omic   "MetaG", "MetaB" or "HostT"  to define output file name and column name in output files
#' @param  Treat.omic.input   Treat omic module quantity file as a dataframe with rownames being features and colnames being samples; or as path to the file generated from dimention reduction functions
#' @param  Mediator.omic.input   Mediator omic module quantity file as a dataframe with rownames being features and colnames being samples; or as path to the file generated from dimention reduction functions
#' @param  meta.mediate   metadata file containing clinical variables
#' @param  Y   the clinical variable to which the impact by Treat omic modules and Mediator omic modules
#' @param  Treat.omic.sigModules   Treat omic modules generated from the glm.sigModules function 
#' @param  Mediator.omic.sigModules    Mediator omic modules generated from the glm.sigModules function
#' @param  outputDir   output directory, default "mediation.out"
#' @param  log.file   log file, default "mediation.log"
#' 
#' 
#' @examples
#' MediationAnalysis(Treat.omic = "MetaG", Mediator.omic = "MetaB", Treat.omic.input = "path/MetaG/quantity/file", 
#'                   Mediator.omic.input = "path/MetaB/quantity/file", meta.mediate = "path/meta/file", Y = "diseaseSeverity",
#'                   Treat.omic.sigModules = MetaG.sigMods, Mediator.omic.sigModules = MetaB.sigMods)




MediationAnalysis <- 
  function(
    Treat.omic,
    Mediator.omic,
    Treat.omic.input,
    Mediator.omic.input,
    meta.mediate,
    Y,
    Treat.omic.sigModules,
    Mediator.omic.sigModules,
    outputDir = "mediation.out",
    log.file = "mediation.log"
  ){
    
    library(mediation)
    if(class(Treat.omic.input) == "character"){
      if(grepl("\\.gct$", Treat.omic.input)) m1 <- parse.gctx(Treat.omic.input)@mat %>% t() %>% data.frame() else m1 <- data.frame(fread(Treat.omic.input), row.names=1)
    }else{
      m1 = Treat.omic.input
      feature.abb_df1 <- cbind.data.frame(feature = rownames(m1),
                                         abb = paste("feature",seq(1,nrow(m1),1),sep = ""),
                                         stringsAsFactors = F)
      rownames(m1) <- sapply(rownames(m1), function(x) feature.abb_df1$abb[which(feature.abb_df1$feature == x)])
      m1 <- t(m1) %>% as.data.frame(stringsAsFactors=F)
      colnames(m1) <- sapply(colnames(m1), function(x) feature.abb_df1$feature[which(feature.abb_df1$abb == x)])
    }
    
    if(class(Mediator.omic.input) == "character"){
      if(grepl("\\.gct$", Mediator.omic.input)) m2 <- parse.gctx(Mediator.omic.input)@mat %>% t() %>% data.frame() else m2 <- data.frame(fread(Mediator.omic.input), row.names=1)
    }else{
      m2 = Mediator.omic.input
      feature.abb_df2 <- cbind.data.frame(feature = rownames(m2),
                                          abb = paste("feature",seq(1,nrow(m2),1),sep = ""),
                                          stringsAsFactors = F)
      rownames(m2) <- sapply(rownames(m2), function(x) feature.abb_df2$abb[which(feature.abb_df2$feature == x)])
      m2 <- t(m2) %>% as.data.frame(stringsAsFactors=F)
      colnames(m2) <- sapply(colnames(m2), function(x) feature.abb_df2$feature[which(feature.abb_df2$abb == x)])
      
    }
   
    meta_df <- fread(meta.mediate, data.table = F)
    
    if(!dir.exists(outputDir)) dir.create(outputDir)
    
    cat(paste(as.character(Sys.time()), '\n'),  file=log.file, append=T)
    cat('Performing mediation analysis: \n',  file=log.file, append=T)
    
    Mediation.results <- NULL
    for(md1 in Treat.omic.sigModules){
      
      for(md2 in Mediator.omic.sigModules){
        
        id = paste(Treat.omic, ".", md1, "_", Mediator.omic, ".", md2, "_", Y,sep = "" )
        cat(paste(id, "\n", sep = ""),  file=log.file, append=T)
        
        dat <- base::merge(base::merge(meta_df, m1 %>% dplyr::select(all_of(md1)), by.x = "SampleID", by.y = 0),
                           m2 %>% dplyr::select(all_of(md2)), by.x = "SampleID", by.y=0)
        
        colnames(dat)[(ncol(dat)-1):ncol(dat)] <- c("Treat", "Mediator")
        
        # write.table(dat, file = paste(outputDir, "/",Treat.omic, "_", md1, "_", Mediator.omic, "_", md2,".txt",sep = ""),
        #             quote = F, row.names = F, sep = "\t")
        
        colnames(dat)[which(colnames(dat) == Y)] <- "Y"
        
        # remove na
        
        dat <- dat %>% filter(!is.na(Y)) %>% filter(!is.na(Treat)) %>% filter(!is.na(Mediator))
        
        #model.m=lm(Mediator ~ Treat+Age+Gender+CurrentSmoking+ICS,dat)
        model.m = lm(as.formula( paste("Mediator ~ Treat + ", 
                                       paste(colnames(dat)[!(colnames(dat) %in% c("SampleID","Y","Mediator","Treat"))], collapse = " + "),
                                       sep = "") ), data = dat)
        
        #model.y=lm(Y~Treat+Mediator+Age+Gender+CurrentSmoking+ICS,dat)
        model.y = lm(as.formula( paste("Y ~ Treat + Mediator + ", 
                                       paste(colnames(dat)[!(colnames(dat) %in% c("SampleID","Y","Mediator","Treat"))], collapse = " + "),
                                       sep = "") ), data = dat)
        
        summary = summary(mediate(model.m,model.y,treat="Treat",mediator="Mediator",boot=F,sims=1000))
        #capture.output(summary,file="mediator_out.txt",append=FALSE)
        res <- capture.output(summary,append=FALSE)
        
        #sub( "^()\\s", "\\1", res[7])
        tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
        tmp <- tmp[tmp != "" & tmp!="."]
        tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
        ACME.p <- tmp[length(tmp)]
        
        tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
        tmp <- tmp[tmp != "" & tmp!="."]
        tmp <- tmp[!grepl("*",tmp,fixed = T) ]
        ADE.p <- tmp[length(tmp)]
        
        tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
        tmp <- tmp[tmp != "" ]
        i_str = which(grepl("Mediated", tmp))
        prop.mediated <- tmp[(i_str + 1)]
        
        spearman.r = cor(dat$Treat, dat$Mediator, method = "spearman")
        
        
        
        vec = c(id, ACME.p, ADE.p, prop.mediated, spearman.r)
        names(vec) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated", "spearman.r")
        
        Mediation.results <- bind_rows(Mediation.results, vec)
        
      }
    }
    
    cat('Mediation analysis generates a result file named Treat_affects_Y_through_Mediator.txt. \n ',  file=log.file, append=T)
    write.table(Mediation.results, file = paste(outputDir, "/", Treat.omic,"_affects_", Y, "_through_", Mediator.omic,".txt",sep = "" ),
                quote = F, row.names = F, sep = "\t")
    
    return(Mediation.results)
  }