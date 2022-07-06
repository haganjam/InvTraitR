#'
#' @title downstream_itis()
#' 
#' @description Get downstream taxa from order from the itis database
#' 
#' @details This function is used to get all downstream taxa from an order and process it
#' into a matrix of taxonomic distances and a vector of taxon names
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param ord.id - taxon ID of the order
#' @param ord.name - name of the order
#' 
#' @return list with two elements: 
#' 1. tax_distance - taxonomic distance matrix (saved as a sparse matrix); 
#' 2. tax_names - vector of all the names in the taxonomic distance matrix
#' 

downstream_itis <- function(ord.id, ord.name) { 
  
  # load relevant packages
  library(dplyr)
  library(taxize)
  
  # get downstream with error handling and multiple tries
  source(here("scripts/create_database/02_get_downstream_taxa_function.R"))
  
  # get downstream taxa and process into a usable data.frame
  downtax.top <- get_downstream_taxa(sci_id = ord.id, downto = "genus", db = "itis", intermediate = TRUE)
  
  # if the function cannot access downstream taxa then return NAs
  if (identical(downtax.top, NULL)) {
    return(list(tax_distance = NA, tax_names = NA))
  }
  
  # bind the downstream taxa into a data.frame
  downtax.top <- bind_rows(downtax.top[[1]]$intermediate)
  
  # get taxonomic distance matrix
  source(here("scripts/create_database/05_taxon_matrix_itis_function.R"))
  
  # extract the missing ranks
  missing.ranks <- rbind(downtax.top[downtax.top$rankname != "genus", c("rankname", "taxonname")], c("order", ord.name) )
  names(missing.ranks) <- c("parentrank", "parentname")
  
  # join the missing ranks
  downtax.top <- left_join(downtax.top, missing.ranks, by = "parentname")
  
  # remove species names as they should not be searched
  downtax.top <- 
    downtax.top %>%
    filter(rankname != "species")
  
  # apply taxonomic weights
  weights <- mapply(function(x, y) tax.itis[which(row.names(tax.itis) == x), which(colnames(tax.itis) == y) ],
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
  
  # convert symmetrical values in upper matrix to zeros
  d.g.dist [upper.tri(d.g.dist , diag = FALSE)] <- 0
  
  # convert the distance matrix into a sparse matrix
  library(Matrix)
  d.g.dist  <- Matrix(d.g.dist, sparse = TRUE)
  
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
