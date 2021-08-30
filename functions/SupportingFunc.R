

#' Read.GeneSets.db2 function
#'
#' This function allows you to read gene set database file in gmt format 
#' @param  gs.db  input data file in gct format, first column (Name) must contain gene symbols
#' @param  thres.min prefix used for output tables
#' @param  thres.max list of genesets (in gmt format) to evaluate enrichment on
#' @examples 
#' Read.GeneSets.db2("db/file/path")





#############################################################
##
##  import gene sets into R workspace
##
## 20161013 modified by kk
#############################################################
Read.GeneSets.db2 <- function (gs.db, thres.min = 2, thres.max = 2000) {
  ## read gmt files
  temp <- readLines(gs.db)
  temp <- strsplit(temp, '\t')
  temp.size.G <- sapply(temp, function(x) length(x)-2)
  
  ## filter gene sets according to size
  rm.idx <- which(temp.size.G < thres.min | temp.size.G > thres.max)
  if(length(rm.idx) > 0){
    temp <- temp[-rm.idx]
    temp.size.G <- temp.size.G[-rm.idx]
  }
  
  max.Ng <- length(temp)         ## number of signature sets
  temp.size.G <- sapply(temp, function(x) length(x)-2)
  max.size.G <- max(temp.size.G) ## maximal size
  
  gs <- lapply(temp, function(x)x[3:length(x)])
  gs.names <- sapply(temp, function(x)x[1])
  gs.desc <- sapply(temp, function(x)x[2])
  
  ## check whether gene sets are unique
  gs.unique <- lapply(gs, unique)
  gs.unique.size.G <- sapply(gs.unique, length)
  gs.not.unique.idx <- which(gs.unique.size.G < temp.size.G)
  if( length(gs.not.unique.idx) > 0 ){
    warning("\n\nDuplicated gene set members detected. Removing redundant members from:\n\n", paste(gs.names[gs.not.unique.idx], collapse='\n'))
    gs <- gs.unique
    temp.size.G <- gs.unique.size.G 
  }
  size.G <- temp.size.G
  names(gs) <- names(gs.names) <- names(gs.desc) <- names(size.G) <- gs.names
  
  
  return(list(N.gs = max.Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc,
              size.G = size.G, max.N.gs = max.Ng))
}


#' Read.GeneSets.db2 function
#'
#' 
#' @param  string  
#' @param  nChar 
#' @param  add.dots 
#' @examples 
#' chopString("dskdf")






#################################################
##   Given a string and a number of characters
##   the function chops the string to the
##   specified number of characters and adds
##   '...' to the end.
## parameter
##   string     - character
##   nChar      - numeric
## value
##   string of 'nChar' characters followed
##     by '...'
##################################################
chopString <- function(string, nChar=10, add.dots=T)
{
  string.trim <- strtrim(string, nChar)
  
  if(add.dots)
    string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ], '...')
  if(!add.dots)
    string.trim[ which(nchar(string) > nChar) ] <-  paste(string.trim[which(nchar(string) > nChar) ])
  
  return(string.trim)
}


#' Read.GeneSets.db2 function
#' translate a color name into rgb space
#' 
#' @param  color  
#' @param  alpha 
#' @param  maxColorValue 
#' @examples 
#' my.col2rgb("blue")


##########################################################################################################
## translate a color name into rgb space
##
## changelog:  20100929 implementation
##########################################################################################################
my.col2rgb <- function(color, alpha=80, maxColorValue=255){
  
  out <- vector( "character", length(color) )
  
  for(col in 1:length(color)){
    
    col.rgb <- col2rgb(color[col])
    
    out[col] <- rgb(col.rgb[1], col.rgb[2], col.rgb[3], alpha=alpha, maxColorValue=maxColorValue)
    
  }
  return(out)
}

