#' LOSO.ko function
#'
#' This function allows you to perform LOSO (leave one species out ) analysis. 
#' Users should provide 1) gene quantification data, 2) information about Species annotaion of bins, 3) connection between bins and genes through scaffold ids, and 4) KO annotaion of genes. 
#' This function generates one KO abundance file each time after excluding one Species from the gene quantification data.  
#' Resulting gct files are stored in the "LOSO_ko.abund" directory.
#' 
#' @param  geneDepth_df  input file path of gene quantification, with the first column being gene ids which should be formated as "scaffold_number"
#' @param  gene.ko_file   path to file giving informaiton on gene ids and KO number. should contain columns named "Gene" and  "KOnumber". gene ids should be formated in "scaffold_number"
#' @param  bin.scaffold_file path to file giving information on bins and scaffold id, should contain columns named "Scaffold" and "Bin"
#' @param  bin.species_file  path to file giving information on bins and Species annotaion, should contain columns named "Bin" and "Species"
#' 
#' @examples
#' LOSO.ko(geneDepth_df = gene.Depth_df, 
#'         gene.ko_file="source.data/ko_noeuk.txt", 
#'         bin.scaffold_file =  "source.data/all_membership.txt",
#'         bin.species_file = "source.data/bin_species.txt")
#' 
#' 
#' 


LOSO.ko <- function(
  geneDepth_file = "F:/temp.data/geneDepth.txt",  
  gene.ko_file = "source.data/ko_noeuk.txt",
  bin.scaffold_file =  "source.data/all_membership.txt",
  bin.species_file = "source.data/bin_species.txt",
  log.file='LOSO.ko.log'
){
  
  ## #####################################################
  ##
  ## load data and delete unused data to save memory
  ##
  ## #####################################################
  cat(paste(as.character(Sys.time()), '\n'),  file=log.file, append=T)
  cat('Performing LOSO.ko analysis: \n',  file=log.file, append=T)
  cat('Importing data: \n',  file=log.file, append=T)
  
  library(data.table)
  library(tibble)
  library(reshape2)
  library(dplyr)
  
  gene.KO_df <- fread(gene.ko_file,  header = T) %>% dplyr::select(Gene, KOnumber)
  
  geneDepth_df <- fread(geneDepth_file, data.table = F) 
  colnames(geneDepth_df)[1] <- "GeneID"
  
  geneDepth_df <- geneDepth_df %>%   # ,nrows = 300 to test
    dplyr::filter(GeneID %in%  gene.KO_df$Gene) %>% 
    tibble::column_to_rownames("GeneID") 
  
  gc()
  
  geneDepth_df <- merge(geneDepth_df,
                        gene.KO_df, by.x = 0, by.y="Gene") 
  colnames(geneDepth_df)[which(colnames(geneDepth_df) == "Row.names")] <- "Gene"
  remove(gene.KO_df)
  gc() 
  
  
  
  ## #####################################################
  ##
  ## create species - scaffold file 
  ##
  ## #####################################################
  bin.scaf_df <- fread(bin.scaffold_file) %>% dplyr::select(Scaffold, Bin)
  bin.species_df <- fread(bin.species_file) %>% dplyr::select(Bin, Species)
  
  species.scaffold_df <- merge(bin.species_df, bin.scaf_df, by="Bin")
  remove(bin.scaf_df, bin.species_df)
  gc() 
  
  if(!dir.exists("5_LOSO") ) dir.create("5_LOSO")
  
  uniqSpecies = unique(species.scaffold_df$Species)
  rows.scaf = sub("(_\\d+$)","", geneDepth_df$Gene)
  
  cat(paste('Excluding', length(uniqSpecies),' species one by one: \n'),  file=log.file, append=T)
  for(species in uniqSpecies){
    i=which(uniqSpecies == species)
    cat(paste(i,".  ", species, ' \n', sep = ""),  file=log.file, append=T)
    
    #  species = uniqSpecies[3]
    
    scafs = species.scaffold_df$Scaffold[species.scaffold_df$Species == species]
    
    
    geneDepth_df_sub <- geneDepth_df[!(rows.scaf %in% scafs),]
    
    res <- geneDepth_df_sub %>%
      dplyr::select(-Gene) %>%
      # reshape2::melt(id.var = "KO", )
      dplyr::group_by(KOnumber) %>%
      dplyr::summarise_if(is.numeric, sum, na.rm = TRUE) %>%
      tibble::column_to_rownames("KOnumber")
    
    
    ## export gct
    gct.tmp <- new('GCT')
    gct.tmp@mat <- data.matrix( res )
    gct.tmp@rid <- rownames(res)
    gct.tmp@cid <- colnames(res)
    gct.tmp@cdesc <- data.frame(id=colnames(res))
    gct.tmp@rdesc <- data.frame(row.names = rownames(res), Description = rownames(res))
    fnn = paste("5_LOSO/ko.abund_rm.", species, ".gct", sep = "")
    gct.tmp@src <- fnn
    
    write.gct(gct.tmp, ofile=fnn, appenddim = F)  
  }
  
}
