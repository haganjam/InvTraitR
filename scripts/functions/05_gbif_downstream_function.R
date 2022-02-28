
# Get downstream taxa gbif function: downstream_gbif

# args:
# ord.id - taxon ID of the order of the equation name
# ord.name - name of the order of the equation name

downstream_gbif <- function(ord.id, ord.name) {
  
  # load relevant packages
  library(dplyr)
  library(taxize)
  
  # get downstream with error handling and multiple tries
  source(here("scripts/functions/04_get_downstream_taxa_function.R"))
  
  # get downstream taxa and process into a usable data.frame
  downtax.top <- get_downstream_taxa(sci_id = ord.id, downto = "family", db = "gbif", intermediate = FALSE)
  downtax.top <- downtax.top[[1]]
  downtax.top$parentname <- ord.name
  downtax.top$parentrank <- "order"
  
  # loop over all families and get their downstream taxa
  downtax.int <- vector("list", length = nrow(downtax.top))
  for (i in 1:nrow(downtax.top)) {
    
    x <- get_downstream_taxa(sci_id = downtax.top$key[i], downto = "genus", db = "gbif", intermediate = FALSE)
    if (length(x[[1]]) == 0 ) {
      x <- NULL
    } else {
      x <- x[[1]]
      x$parentname <- downtax.top$name[i]
      x$parentrank <- "family"
      downtax.int[[i]] <- x
    }
  }
  
  # bind the intermediate taxa into a single data.frame
  downtax.top <- bind_rows(bind_rows(downtax.top), (downtax.int))
  
  # create the distance matrix
  d.mat <- 
    downtax.top %>%
    select(from = parentname, to = taxonname)
  
  # use igraph to create a graph from the matrix
  d.g <- graph_from_data_frame(d = d.mat, directed=FALSE)
  
  # produce a distance matrix using the taxonomic weights
  d.g.dist <- distances(
    d.g,
    v = V(d.g),
    to = V(d.g),
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
