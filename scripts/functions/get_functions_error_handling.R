
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# Modified get functions for better error handling

# args:

# database_function: name of the function - itis (get_tsn), bold - (get_boldid), gbig - (get_gbifid)
# taxon_name: name that you want to query
# ask_or_not: accept imperfect matches (default is FALSE)

get_taxon_id <- function(database_function = "get_tsn", taxon_name, ask_or_not = FALSE) {
  
  x <- "try_error"
  i <- 1
  while( class(x) == "try-error" ) {
    
    if (i > tries){
      break
    }
    
    i <- i + 1
    x <- try( do.call("get_boldid", list(sci = taxon_name, ask = ask_or_not)) )
    
  }
  
  if ( class(x) == "try-error" ) {
    
    attr(x, "match") <- "could not access database"
    warning("could not access database, use more tries or run function at a later stage")
    
  }
  
  return(x)
  
}

### END
