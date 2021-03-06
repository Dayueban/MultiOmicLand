#' MetaB.HostT.link function
#'
#' This function allows you to identify MetaB - HostT module pairs with biological links.  
#' Dependent packages includes "dplyr" and "data.table". 
#' 
#' @param  mediation.res   resulting data frame generated from MediationAnalysis() function
#' @param  MetaB_quantity_file   metabolomic quantification file on the feature level
#' @param  MetaB_module.feature_file   module assignment file of the metabolomic featrues, generated by the wgcna function
#' @param  HostT_quantity_file   host transcritomic quantificaiton file on the feature level
#' @param  HostT_module.feature_file  module assignment file of the host transcritomic featrues, generated by the wgcna function
#' @param  METABO2CIDm_file    match file between metabo ID and CIDm ID, generated by the perl scripts provided
#' @param  CIDm.receptor_file  link file between CIDm ID and receptor, available on the githup repo ????
#' @param  ACME.p.co    cutoff of ACME p value to define MetaB-HostT module pairs with potential biological links
#' @param  CIDm.receptor.score.co  ????
#' @param  output.dir
#' @param  output.prefix
#' @param  log.file
#' 
#' @examples
#' MetaB.HostT.link(mediation.res = MetaB.HostT.NUE_medres, 
#'                  MetaB_quantity_file = "source.data/metabolome_sub.txt", 
#'                  MetaB_module.feature_file = "DR_wgcna/metaB.module_assign.txt",
#'                  HostT_quantity_file = "source.data/transcriptome.txt",
#'                  HostT_module.feature_file = "DR_wgcna/hostT.module_assign.txt" ,
#'                  METABO2CIDm_file =  "database/metabo2CIDm.txt",
#'                  CIDm.receptor_file = "database/all_cidm_receptor.txt",
#'                  ACME.p.co = 0.25,
#'                  CIDm.receptor.score.co = 700 )
#'                  
  

MetaB.HostT.link <- function(
  mediation.res = MetaB.HostT.NUE_medres,  
  MetaB_quantity_file = "source.data/metabolome_sub.txt",  
  MetaB_module.feature_file = "DR_wgcna/metaB.module_assign.txt", 
  HostT_quantity_file = "source.data/transcriptome.txt", 
  HostT_module.feature_file = "DR_wgcna/hostT.module_assign.txt", 
  METABO2CIDm_file =  "database/metabo2CIDm.txt",  
  CIDm.receptor_file = "database/all_cidm_receptor.txt", 
  ACME.p.co = 0.25,
  CIDm.receptor.score.co = 700,
  output.dir = "biological.links",
  output.prefix = "",
  log.file = "MetaB.HostT.link.log"
){
  
  
  cat(paste("\n\n",as.character(Sys.time()), '\n'),  file=log.file, append=T)
  cat("Importing data : \n", file=log.file, append=T)
  
  
  ## ###############################################
  ##
  ##  import data 
  ##
  ## ###############################################
  library(data.table)
  library(dplyr)
  
  MetaB.dat <- data.frame(fread(MetaB_quantity_file), row.names = 1) %>% t() %>% as.data.frame()
  MetaB_module.feature <- fread(MetaB_module.feature_file, data.table = F, col.names =c("Feature","Module"))
  HostT.dat <- data.frame(fread(HostT_quantity_file), row.names = 1) %>% t() %>% as.data.frame()
  HostT_module.feature <- fread(HostT_module.feature_file, data.table = F, col.names =c("Feature","Module"))
  metabo2CIDm <- fread(METABO2CIDm_file, data.table = F, header = F, col.names = c("Metabo","CIDm"))
  
  # load and organize CIDm.receptor data 
  conn <- file(CIDm.receptor_file, open="r")
  linn <-readLines(conn)
  close(conn)
  
  
  numElements <- sapply(linn, function(x) length(strsplit(x, "\t", fixed = T)[[1]]))
  table(numElements) # 4 levels of lengths
  
  tmp <- linn[numElements == 4]
  subdf1 <- unname(sapply(tmp,function(x) strsplit(x,"\t")[[1]] )) %>% t() %>% data.frame(stringsAsFactors = F)
  colnames(subdf1) <- c("CIDm","ENSP","geneName","linkType")
  
  
  tmp <- linn[numElements == 5]
  subdf2 <- unname(sapply(tmp,function(x) strsplit(x,"\t")[[1]] )) %>% t() %>% data.frame(stringsAsFactors = F)
  colnames(subdf2) <- c("CIDm","ENSP","geneName","description", "linkType")
  
  tmp <- linn[numElements == 9]
  subdf3 <- unname(sapply(tmp,function(x) strsplit(x,"\t")[[1]] )) %>% t() %>% data.frame(stringsAsFactors = F)
  colnames(subdf3) <- c("CIDm","ENSP","geneName","linkType","X1","X2","X3","X4", "score") 
  
  tmp <- linn[numElements == 10]
  subdf4 <- unname(sapply(tmp,function(x) strsplit(x,"\t")[[1]] )) %>% t() %>% data.frame(stringsAsFactors = F)
  colnames(subdf4) <- c("CIDm","ENSP","geneName","description", "linkType","X1","X2","X3","X4","score") 
  
  CIDm.receptor <- bind_rows(subdf1, subdf2, subdf3, subdf4) %>% filter(score >= CIDm.receptor.score.co)

  
  # merge data
  MetaB.HostT.dat <- merge(MetaB.dat, HostT.dat, by=0)
  
  
  ## ###############################################
  ##
  ##  perform linke analysis 
  ##
  ## ###############################################
  cat("Performing link analysis : \n", file=log.file, append=T)
  # identify MetaB-HostT module pairs 
  #MetaB.HostT.modPairs <- strsplit((mediation.res %>% dplyr::filter(ACME.p <= ACME.p.co))$Treat_Mediator_Y,"_",fixed = T)
  MetaB.HostT.modPairs <- (mediation.res %>% dplyr::filter(ACME.p <= ACME.p.co))$Treat_Mediator_Y
  
  
  # identify links
  links <- NULL
  for(mdp in MetaB.HostT.modPairs){
    # mdp = MetaB.HostT.modPairs[[1]]
    
    ConfirmedLink = FALSE
    #metab.md = sub("ME(.*)$","\\1",sub("MetaB\\.(.*)$","\\1",mdp[1]) )
    #hostt.md = sub("ME(.*)$","\\1",sub("HostT\\.(.*)$","\\1",mdp[2]) )
    
    parts = strsplit(mdp,"_", fixed = T)[[1]]
    metab.md = sub("ME(.*)$","\\1",sub("MetaB\\.(.*)$","\\1",parts[1]) )
    hostt.md = sub("ME(.*)$","\\1",sub("HostT\\.(.*)$","\\1",parts[2]) )
    
    #cat(paste("\n\nAnalyzing module pairs: ", parts[1], "_",parts[2],"\n", sep = ""), file=log.file, append=T)
    i_mdp <- which(MetaB.HostT.modPairs == mdp)
    if(i_mdp %% 1000 == 0) cat(paste(" -------------------------- progress: ", i_mdp, " out of ", length(MetaB.HostT.modPairs)," module pairs -------------------------- \n", sep = ""), file=log.file, append=T)
    
    
    
    
    metab.ftr = MetaB_module.feature$Feature[MetaB_module.feature$Module == metab.md]
    hostt.ftr = HostT_module.feature$Feature[HostT_module.feature$Module == hostt.md]
    
    
    for(bf in metab.ftr){
      # bf=metab.ftr[2]
      # writeLines(paste("the ", which(metab.ftr == bf), " bf", sep = "") )
      
      if(ConfirmedLink ) break
      cidm = metabo2CIDm$CIDm[which(metabo2CIDm$Metabo == bf)]
      
      if(length(cidm) == 0) {
        writeLines(paste("MetaB feature ",bf, " (belonged to the ", parts[1]," module) doesn't have matching CIDm", sep=""))
        next
      }
      
      if(length(cidm) > 1) cidm <- as.character((data.frame(table(cidm)) %>% dplyr::arrange(desc(Freq)))$cidm[1])
      
      receptors_df = CIDm.receptor %>% filter(CIDm == cidm)
      receptors <- receptors_df$geneName
      
      receptors <- receptors[receptors != ""]
      
      if( length(receptors) == 0 ) next 
      if(!any(receptors %in% hostt.ftr)) next
      
      # verify the links
      recptr_check <- receptors[receptors %in% hostt.ftr]
      receptors_df <- receptors_df %>% filter(geneName %in%  recptr_check)
      
      for(i in c(1:nrow(receptors_df))){
        if(ConfirmedLink ) break
        rcpt <- receptors_df$geneName[i]
        lkType <- receptors_df$linkType[i]
        
        if(lkType %in% c("catalysis","reaction", "expression")) next
        
        if(lkType == "binding"){
          link = c(paste(parts[1], parts[2], sep = "_" ), bf, rcpt )
          names(link) <- c("MetaB.HostT_modulePair", "MetaB.feature","HostT.feature")
        }else {
          
          if(all(c(bf, rcpt) %in% colnames(MetaB.HostT.dat))) cor.dat <- MetaB.HostT.dat %>% dplyr::select(all_of(c(bf, rcpt))) else next
          test.r = cor(cor.dat[,1], cor.dat[,2], method = "spearman")
          
          if(lkType == "activation" & test.r > 0){
            link = c(paste(parts[1], parts[2], sep = "_" ), bf, rcpt )
            names(link) <- c("MetaB.HostT_modulePair", "MetaB.feature","HostT.feature")
            
          }else if(lkType == "inhibition" & test.r < 0){
            link = c(paste(parts[1], parts[2], sep = "_" ), bf, rcpt )
            names(link) <- c("MetaB.HostT_modulePair", "MetaB.feature","HostT.feature")
          }else next  #all other scenarios don't count as a biological link
          
        }
        
        
        links <- bind_rows(links,link)
        remove(link)
        ConfirmedLink = T
        
      } # loop through all the receptors to be checked 
      
    }# loop through metab.features
    
    
    
  } # loop through module pairs
  
  if(!dir.exists(output.dir)) dir.create(output.dir)
  write.table(links, file = paste(output.dir,"/",output.prefix,"MetaB.HostT.modules.linked.txt",sep = ""), sep = "\t", quote = F, row.names = F)
  
  
  return(links)
}

