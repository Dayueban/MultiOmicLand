#' ko2metabo function 
#'
#' This function generates information on substrates, products and substrate&products (in the format of metabolomic ids) for each KO number. 
#' Dependent packages includes "data.table" and "dplyr". 
#' Users should provide the path to a database directory, which must contain the cmpd2metabo.txt file and KO2CMPD.lists.RData  file.
#' This function should be run after the ko2cmpd() function
#' Resulting files in the format of RData will be stored in the same database directory.
#' 
#' @param  dbDir  database directory containing the "cmpd2metabo.txt" and "KO2CMPD.lists.RData" files
#' 
#'  
#' @examples
#' ko2metabo(dbDir = "path/of/database/dir" )
#' 



ko2metabo <- function(
  dbDir = "database"
){
  
  
 ##  generate KO.Substrates.metab, KO.Products.metab and KO.SubsProd.metab lists
  
  cpmd2metab <- fread(paste(dbDir,"/cmpd2metabo_idmatch_100_99.txt",sep = ""),data.table = F )
  load(paste(dbDir,"/KO2CMPD.lists.RData",sep = ""))
  
  KO.Substrates.metab_list <-  vector("list", length = length(KO.Substrates_list))
  KO.Products.metab_list <-  vector("list", length = length(KO.Substrates_list))
  KO.SubsProd.metab_list <-  vector("list", length = length(KO.Substrates_list))
  
  for(i in c(1:length(KO.Substrates_list))){
    
    ko <- names(KO.Substrates_list)[i]
    
    substrates <- KO.Substrates_list[[ko]]
    products <- KO.Products_list[[ko]]
    subs.prod <- KO.SubsProd_list[[ko]]
    
    substrates.metab <- unlist(sapply(substrates, function(x) cpmd2metab$Metabolome[cpmd2metab$MetaCyc_compounds == x]))
    products.metab <- unlist(sapply(products, function(x) cpmd2metab$Metabolome[cpmd2metab$MetaCyc_compounds == x]))
    subs.prod.metab <- unlist(sapply(subs.prod, function(x) cpmd2metab$Metabolome[cpmd2metab$MetaCyc_compounds == x]))
    
    if(!is.null(substrates.metab) & length(substrates.metab) > 0 )  KO.Substrates.metab_list[[i]] <- substrates.metab; names(KO.Substrates.metab_list)[i] <- ko
    if(!is.null(products.metab) & length(products.metab) > 0 )  KO.Products.metab_list[[i]] <- products.metab; names(KO.Products.metab_list)[i] <- ko
    if(!is.null(subs.prod.metab) & length(subs.prod.metab) > 0 )  KO.SubsProd.metab_list[[i]] <- subs.prod.metab; names(KO.SubsProd.metab_list)[i] <- ko
  }
  
  save(KO.Products.metab_list, KO.SubsProd.metab_list, KO.Substrates.metab_list,
       file = paste(dbDir, "/KO2METABO.lists.RData",sep = ""))
  
}