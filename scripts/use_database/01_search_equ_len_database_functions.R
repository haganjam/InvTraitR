
# Function to replace target.name with accepted name

# args
# target.name - name to test for synonyms
# data.base - taxonomic database to use (itis, gbif)

syn_correct <- function(target.name, data.base = "itis") {
  
  # load the relevant functions
  source(here("scripts/functions/01_get_taxon_id_function.R"))
  
  # if the input name is a species then extract the genus
  search.name <- extract_genus(binomial = target.name)
  
  # extract target.name ID
  syn.id <- get_taxon_id(database_function = data.base, 
                         taxon_name = search.name, ask_or_not = FALSE, tries = 5)
  
  # if the search.name is not in the database, then throw a warning
  if (is.na(syn.id)) {
    message( paste("Target taxon name is not in the ", data.base, " database", sep = "" ) )
  }
  
  # extract synonyms associated with the target.name ID
  syn.df <- synonyms(syn.id)[[1]]
  
  # if the target.name is not the accepted name then replace it with the accepted name
  if ( any(is.na(syn.df)) ) { 
    
    print(paste("No synonymns for this taxon name in the ", data.base, " database", sep = "" ))
    
  }  else if ( any(syn.df$sub_tsn != syn.df$acc_tsn) ) { 
    
    warning(paste(nrow(syn.df), " synonyms: ", paste(syn.df$syn_name, collapse = ", "), " | Using accepted name: ", syn.df$acc_name[1], sep = ""  ))
    target.name <- syn.df$acc_name[1] 
    search.name <- extract_genus(binomial = target.name)
    
  } else { 
    
    print(paste("Synonymns present but not accepted in the ", data.base, " database", sep = "" )) 
    
  }
  
  return(search.name)
  
}

# Function to get best id from either the length or the equation data

# args
# target.name - name of the taxon to get the equation and length for
# id_info - e_ti/l_ti
# equ_len - equation or length data
# length_only - TRUE
# data.base - taxonomic database to use (itis, gbif)
# max_tax_dist - 6
# d.dist - taxonomic distance matrices of the orders

get_tax_distance <- function(target.name, id_info, equ_len, length_only = TRUE,
                             data.base = "itis", max_tax_dist = 6,
                             d.dist) {
  
  # use the synonymn function here to output an object called "search.name"
  search.name <- syn_correct(target.name = target.name, data.base = data.base)
  
  # if length.only = TRUE then subset equations with only length data
  if (length_only & (equ_len == "equation") ) {
    
    id_info1 <- id_info[ sapply(id_info, function(x) x$id %in% equ_id$id_only_equ_ID) ]
    
  } else { id_info1 <- id_info }
  
  # extract the orders present in the equation database
  orders <- sapply( id_info1, function(x) x$order)
  
  # subset the distance matrices for the equations
  dist <- d.dist[ sapply(d.dist, function(x) x$order %in% orders) ]
  
  # select the taxonomic distance matrix consistent with the taxonname
  dist <- dist[ sapply(dist, function(x) search.name %in% x$tax_names$taxonname) ]
  
  # select the correct taxonomic distance matrix
  dist.m <- dist[[1]]$tax_distance
  print(paste("dimension: ", paste(dim(dist.m), collapse = " x ") ))
  
  # get the list of taxonomic names
  tax.names <- dist[[1]]$tax_names$taxonname
  
  # extract the order from the list
  id.order <- dist[[1]]$order
  
  # get the equations with the correct order
  y <- sapply(id_info1, function(y) y$order == id.order)
  if (sum(y) == 0) {
    
    print(NA)
    stop(paste("No suitable ", equ_len, " in the data for the target taxon name"))
    
  }
  id_info1 <- id_info1[y]
  
  # get the distance from the suitable equations or lengths
  tax.dist <- 
    
    sapply(id_info1, function(x) { 
      
      if (x$name == target.name) {
        
        d <- 0
        
      } else if ( (search.name %in% tax.names) ) {
        
        z <- extract_genus( x$name )
        d1 <- dist.m[which(row.names(dist.m) == z), which(colnames(dist.m) == search.name) ]
        d2 <- dist.m[which(row.names(dist.m) == search.name), which(colnames(dist.m) == z ) ]
        d3 <- max(c(d1, d2))
        
        # add extra distances for the species level
        if (x$rank == "species") { d3 <- d3+0.25 } # add species level distance
        if (length(unlist( strsplit(x = target.name, split = " ", fixed = TRUE) )) > 1) { d3 <- d3+0.25 }
        
      } 
      
      else {d3 <- NA}
      
      return(d3)
      
    } )
  
  # remove the dist.m object to save memory
  rm(dist.m)
  
  # get the id's of the lowest taxonomic distances
  id.min <- which( near(min(tax.dist), tax.dist) )
  print(tax.dist[id.min])
  
  # check if the minimum taxonomic is within the chosen threshold: max_tax_dist
  if ( any(id.min > max_tax_dist ) ) {
    
    print(NA)
    stop(paste( paste("No suitable ", equ_len, " in the data for the target taxon name"), 
                "as the maximum taxonomic distance is greater than ", max_tax_dist, dep = "") )
    
  }
  
  # add equation range data i.e. min length individual and max length individual...
  
  # choose output type i.e. single mass value or a spreadsheet with options?
  
  # get the taxonomic rank table
  tax.hier <- c("order", "suborder", "infraorder", "section", "subsection", "superfamily",
                "family", "subfamily", "tribe", "subtribe", "genus")
  
  # if there are more than one suitable equation with the same minimum taxonomic distance
  # choose the higher taxonomic level
  if (length(id.min) > 1 ) {
    
    ranks <- sapply(id_info1[id.min], function(x) { x$rank } )
    x <- which(tax.hier %in% ranks )
    y <- which(x == min(x))
    best.id <- id_info1[id.min[y]][[1]]$id 
    
  } else {
    
    best.id <- id_info1[id.min][[1]]$id
    
  }
  
  # print the best.equation ID
  print(paste("Best fitting id: ", best.id) )
  
  return(best.id)
  
  }

# test the function
# id.out <- get_tax_distance(target.name = "Toxomerini", 
                           # id_info = l_ti, 
                           # equ_len = "length", 
                           # length_only = FALSE,
                           # data.base = "itis", max_tax_dist = 3,
                           # d.dist = d.dist
                           # )

### END

