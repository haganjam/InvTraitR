
# Test for packages and provide warnings for the versions

warning_messages <- function(x) {
  
  # get installed packages
  in.packs <- installed.packages()[,c(1, 3)]
  in.packs <- as.data.frame(in.packs)
  row.names(in.packs) <- NULL
  
  # set the packages used
  packs.used <- c("taxize", "dplyr", "igraph", "here")
  
  # generate warnings about the packages that need to be installed
  eval( if (any( !(packs.used %in% in.packs$Package) )) {
    
    # error if packages are not installed
    error.text <- paste("These functions require the following packages to be installed:", 
                        paste(packs.used, collapse = ", "), sep = " " )
    stop(error.text)
    
  } else {
    
    # warning if packages are installed
    warn.packs <- c("dplyr: 1.0.7", "igraph: 1.2.11", "taxize: 0.9.99", "here: 1.0.1")
    
    warn.text <- paste("These functions were written using the version of the following packages:", 
                       paste(warn.packs, collapse = ", "), sep = " " )
    
    # warning message
    message(warn.text)
    
  } )
  
  # check package versions currently installed
  pack.ver <- in.packs[(in.packs$Package %in% packs.used), ]
  pack.ver <- paste(pack.ver$Package, pack.ver$Version, sep = ": ")
  pack.cond <- sort(pack.ver) == sort(warn.packs)
  
  packs.df <- data.frame(installed.packages = sort(pack.ver),
                         original.packages = sort(warn.packs),
                         pack.match = pack.cond)
  
  eval( if( any( !pack.cond ) ) {
    
    message("Currently installed packages do not match original version")
    print(packs.df)
    
  } else {
    
    message("All package versions match with original packages")
    print(packs.df)
    
  } )
  
  # warn about the R-version used to write these 
  message(paste("These functions were written using:", "R version 4.1.2 (2021-11-01)", sep = " "))
  
}

warning_messages()

### END
