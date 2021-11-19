
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# Function to search the taxonomic database and pull out the correct equation ID and taxonomic distance

# args

# focal_taxa_name - name of the taxon that the database will search for
# life_stage - name of the life stage for the taxon name of interest
# rank.difference - acceptable diffference in rank between focal taxon and equation in the database
# length_only - whether the function should only search for equations based on length
# taxon_database - the taxon database that should be searched
# var_input_data - database associated with the taxon_database that has details about the equations

taxon_equation_searcher <- function(focal_taxa_name, life_stage = NA, rank.difference = 1, 
                                    length_only = TRUE, 
                                    taxon_database,
                                    var_input_data) {
  
  # if only length equations are requested, then we subset them out first
  if (length_only == TRUE) {
    
    v <- lapply(split(var_input_data, var_input_data$equation_id), function(y) {
      
      ifelse(!any(y$size_measurement != "body_length"), unique(y$equation_id), NA)
      
    })
    v <- unlist(v, use.names = FALSE)
    
    tx.data.ids <- unlist(lapply(taxon_database, function(x) { x[["equation_id"]] } ), use.names = FALSE)
    tx.data.ids[tx.data.ids %in% v[!is.na(v)]]
    
    taxon_database <- taxon_database[tx.data.ids %in% v[!is.na(v)]]
    
  }
  
  # set-up the while loop starting values
  equ_id <- NA
  pm <- NA
  rank_distance <- NA
  i <- 1
  
  # randomise the order of database to improve sampling efficiency
  l.db <- length(taxon_database)
  n <- sample(x = 1:l.db, l.db, replace = FALSE)
  x <- taxon_database[n]
  
  while ( is.na(equ_id) & ( class(z) != "numeric" ) ) {
    
    if (i == l.db){
      break
    }
    
    # extract information from list
    eq <- x[[i]][["equation_id"]]
    sy <- x[[i]][["synonymns"]]
    ls <- x[[i]][["life_stage"]]
    
    ti <- x[[i]][["taxonomic_information"]]
    
    foc_tax_search <- focal_taxa_name
    
    # if it is a synonymn then replace it with the actual focal name
    if ( foc_tax_search %in% sy ) {
      
      foc_tax_search <- ti[ti$focal_taxa == 1,]$name
      
    }
    
    # search the taxonomic database
    if (foc_tax_search %in% ti[["name"]]) {
      
      rank_distance <- abs( ti[ti[["focal_taxa"]] == 1, ]$rank_number - ti[ti[["name"]] == foc_tax_search, ]$rank_number )
      pm <- ifelse(ti[ti[["focal_taxa"]] == 1, ]$name == foc_tax_search, TRUE, FALSE)
      equ_id <- ifelse(z <= rank.difference, eq, NA)
      
    } else {
      
      rank_distance <- NA
      pm <- NA
      equ_id <- NA
      
    }
    
    # check if it is the correct life-stage
    if ( !is.na(life_stage) & (life_stage == ls)  ) {
      
      equ_id <- equ_id
      
    } else if ( is.na(life_stage) & is.na(ls) ) {
      
      equ_id <- equ_id
      
    } else {
      
      equ_id <- NA
      
    }
    
    i <- i + 1
    
  }
  
  # package the output
  output <- 
    list("equation_id" = equ_id,
         "perfect_match" = pm,
         "taxonomic_difference" = rank_distance)
  
  return(output)
  
}




