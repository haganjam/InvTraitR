
# Test for packages and provide warnings for the versions

# get installed packages
in.packs <- installed.packages()[,c(1, 3)]
in.packs <- as.data.frame(in.packs)
row.names(in.packs) <- NULL

# set the packages used
packs.used <- c("taxize", "dplyr", "igraph", "here")

# generate warnings about the packages that need to be installed
if (any( !(packs.used %in% in.packs$Package) )) {
  
  # error if packages are not installed
  error.text <- paste("These functions require the following packages to be installed:", 
                      paste(packs.used, collapse = ", "), sep = " " )
  stop(error.text)
  
} else {
  
  # warning if packages are installed
  warn.text <- c("dplyr: 1.0.7", "igraph: 1.2.11", "taxize: 0.9.99", "here: 1.0.1")
  
  warn.text <- paste("These functions were written using the version of the following packages:", 
                     paste(warn.text, collapse = ", "), sep = " " )
  
  # warning message
  warning(warn.text)
  
}

# warn about the R-version used to write these 
message(paste("these functions were written using:", "R version 4.1.2 (2021-11-01)", sep = " "))

### END
