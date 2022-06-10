
#' @title get_downstream_taxa()
#' 
#' @description Function to get downstream taxa using the taxize package
#' 
#' @details This is a modified version of the downstream() function from the 
#' taxize package that includes error handling when a name is not found immediately.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param sci_id - taxon id for the database being used. Taxon ids are obtained from
#' the set of get_() functions in taxize. Or from the modified version written for this
#' analysis: get_taxon_id()
#' @param db - database to draw from ("itis", "bold", "gbif")
#' @param downto - taxonomic rank to get downstream taxa until
#' @param intermediate - whether to also get intermediate taxa or not

get_downstream_taxa <- function(sci_id, db,
                                downto, intermediate = FALSE, tries = 5) {
  
  if (is.null(sci_id)) {
    
    x <- NULL
    
  } else {
    
    library(taxize)
    x <- try(stop("!"), silent = TRUE)
    i <- 1
    while( class(x) == "try-error" ) {
      
      if (i > tries){
        break
      }
      
      i <- i + 1
      x <- try( do.call("downstream", list(sci_id = sci_id, db = db,
                                           downto = downto, intermediate = intermediate)),
                silent = TRUE)
      
    }
    
    if ( class(x) == "try-error" ) {
      
      x <- NULL
      warning("could not access downstream taxa, use more tries or run function at a later stage")
      
    }
    
  }
  
  return(x)
  
}

### END
