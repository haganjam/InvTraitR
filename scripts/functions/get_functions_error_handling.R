
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# Modified get functions for better error handling

# args:

# database_function: name of the function - "itis" (get_tsn), "bold" - (get_boldid), "gbif" - (get_gbifid)
# taxon_name: name that you want to query
# ask_or_not: accept imperfect matches (default is FALSE)
# tries: how many times to try the function if it keeps getting an error (default = 5)

get_taxon_id <- function(database_function = "itis", taxon_name, ask_or_not = FALSE, tries = 5) {
  
  # make sure the correct packages are loaded
  if (any( !(c("taxize") %in% installed.packages()[,1]) )) {
    
    stop("error, this functions requires taxize to be installed")
    
    warning("the function was written using taxize_0.9.99")
    
  } else {
    
    warning("you have taxize installed but, as a warning, this function was written using taxize_0.9.99")
    
  }
  
  # choose the correct function based on the database
  if(database_function == "itis") {
    
    func_string <- "get_tsn"
    
  } else if (database_function == "bold") {
    
    func_string <- "get_boldid"
    
  } else if (database_function == "gbif") {
    
    func_string <- "get_gbifid"
    
  } else {
    
    stop("error, choose a supported database: itis, bold, gbif")
    
  }
  
  x <- try(stop("!"), silent = TRUE)
  i <- 1
  while( class(x) == "try-error" ) {
    
    if (i > tries){
      break
    }
    
    i <- i + 1
    x <- try( do.call(func_string, list(sci = taxon_name, ask = ask_or_not)) )
    
  }
  
  if ( class(x) == "try-error" ) {
    
    attr(x, "match") <- "could not access database"
    warning("could not access database, use more tries or run function at a later stage")
    
  }
  
  return(x)
  
}

### END
