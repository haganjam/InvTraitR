
# itis downstream function: downstream_itis

# args:
# ord.id - taxon ID of the order of the equation name
# ord.name - name of the order of the equation name

downstream_itis <- function(ord.id, ord.name) { 
  
  # test for correct packages
  if (any( !(c("taxize", "dplyr", "igraph") %in% installed.packages()[,1]) )) {
    stop("error, this functions requires taxize dplyr and igraph to be installed")
    warning("the function was written using taxize_0.9.99, dplyr_1.0.7 and igraph_1.2.11")
    
  } else {
    warning("you have taxize and dplyr installed but, as a warning, this function was written using taxize_0.9.99 and dplyr_1.0.2")
  }
  message("this function was written using R version 4.0.2 (2020-06-22)")
  
  # load relevant packages
  library(dplyr)
  library(taxize)
  
  # get downstream taxa and process into a usable data.frame
  downtax.top <- downstream(sci_id = ord.id, downto = "genus", db = "itis", intermediate = TRUE)
  downtax.top <- bind_rows(downtax.top[[1]]$intermediate)
  
  # get taxonomic distance matrix
  source(here("scripts/functions/taxon_matrix_itis_function.R"))
  
  # extract the missing ranks
  missing.ranks <- rbind(downtax.top[downtax.top$rankname != "genus", c("rankname", "taxonname")], c("order", ord.name) )
  names(missing.ranks) <- c("parentrank", "parentname")
  
  # join the missing ranks
  downtax.top <- left_join(downtax.top, missing.ranks, by = "parentname")
  
  # apply taxonomic weights
  weights <- mapply(function(x, y) tax.d[which(row.names(tax.d) == x), which(colnames(tax.d) == y) ],
                    x = downtax.top$rankname,
                    downtax.top$parentrank)
  downtax.top$weights <- unlist(weights, use.names = FALSE)
  
  # create the distance matrix
  d.mat <- 
    downtax.top %>%
    select(from = parentname, to = taxonname, weights)
  
  # use igraph to create a graph from the matrix
  d.g <- graph_from_data_frame(d = d.mat, directed=FALSE)
  
  # produce a distance matrix using the taxonomic weights
  d.g.dist <- distances(
    d.g,
    v = V(d.g),
    to = V(d.g),
    weights = d.mat$weights,
    mode = c("all"),
    algorithm = c("bellman-ford")
  )
  
  # get a data.frame of unique taxonomic names and their 
  par.dat <- downtax.top[, c("parentname", "parentrank")]
  names(par.dat) <- c("taxonname", "rank")
  tax.dat <- downtax.top[, c("taxonname", "rankname")]
  names(tax.dat) <- c("taxonname", "rank")
  
  tax.names <- unique( rbind(par.dat, tax.dat) )
  
  # make a list of the distance matrix and the taxonomic names
  return(list(tax_distance = d.g.dist, tax_names = tax.names))
  
}

### END
