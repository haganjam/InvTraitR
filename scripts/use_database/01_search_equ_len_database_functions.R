
### Functions to search the equation and length databases

# load relevant libraries
library(here)
library(dplyr)

## Load the different databases

## Equation data

# load the equation database
if (!exists("equ_id")) {
  equ_id <- readRDS(file = here("database/equation_vars_database.rds"))
}

# load the taxon information database from the equations
if (!exists("e_ti")) {
  e_ti <- readRDS(file = here("database/itis_taxon_identifiers_equation.rds"))
}

## Default length data
if (!exists("len_id")) {
  len_id <- readRDS(file = here("database/default_length_database.rds"))
}

# load the taxon information database from the equations
if (!exists("l_ti")) {
  l_ti <- readRDS(file = here("database/itis_taxon_identifiers_length.rds"))
}

## Taxonomic distance data 

# read in the taxonomic distance database
if (!exists("d.dist")) {
  d.dist <- readRDS(file = here("database/itis_order_taxon_information.rds"))
}


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
    
    message(paste("Synonymns present but not accepted in the ", data.base, " database", sep = "" )) 
    
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
  u <- sapply(dist, function(x) search.name %in% x$tax_names$taxonname)
  dist <- dist[ u ]
  
  # check if the name is found in the order taxonomic matrices
  if (all(u == FALSE)) {
    message(paste("No suitable ", equ_len, " in the data for the target taxon name"))
    message("Returning an NA")
    return(NA)
  }
  
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
    
    message(paste("No suitable ", equ_len, " in the data for the target taxon name"))
    message("Returning an NA")
    return(NA)
    
  }
  id_info1 <- id_info1[y]
  
  # get the distance from the suitable equations or lengths
  tax.dist <- 
    
    sapply(id_info1, function(x) { 
      
      if (x$name == target.name) {
        
        d3 <- 0
        
      } else if ( (search.name %in% tax.names) ) {
        
        z <- extract_genus( x$name )
        d1 <- dist.m[which(row.names(dist.m) == z), which(colnames(dist.m) == search.name) ]
        d2 <- dist.m[which(row.names(dist.m) == search.name), which(colnames(dist.m) == z ) ]
        d3 <- max(c(d1, d2))
        
        # add extra distances for the species level
        if (x$rank == "species") { d3 <- d3+0.25 } # add species level distance
        if (length(unlist( strsplit(x = target.name, split = " ", fixed = TRUE) )) > 1) { d3 <- d3+0.25 }
        
      } 
      
      else {
      
      message(paste("No suitable", equ_len) )
      message("Returning an NA")
      return(NA)
      
      }
      
      return(d3)
      
    } )
  
  # remove the dist.m object to save memory
  rm(dist.m)
  
  # get the id's of the lowest taxonomic distances
  id.min <- which( near(min(tax.dist), tax.dist) )
  print(paste("taxonomic distances: ", paste(tax.dist[id.min], collapse = " ") ) )
  
  # check if the minimum taxonomic is within the chosen threshold: max_tax_dist
  if ( any(id.min > max_tax_dist ) ) {
    
    message(paste( paste("No suitable ", equ_len, " in the data for the target taxon name"), 
                   "as the maximum taxonomic distance is greater than ", max_tax_dist, dep = "") )
    message("Returning an NA")
    return(NA)
    
  }
  
  # create an output data.frame
  dist.df <- data.frame(id = sapply(id_info1[id.min], function(x) x$id),
                        rank = sapply(id_info1[id.min], function(x) { x$rank } ),
                        dist_to_target = tax.dist[id.min])
  
  return(dist.df)

  }


# Function to select the best equation from the set of suitable equations

# args
# target.name - name of the target.taxon
# life.stage - life stage
# id_info - database of id information (l_ti/e_ti)
# equ_len - "equation" or "length"
# data.base - "itis" or "gbif"
# max_tax_dist - maximum acceptable taxonomic difference
# d.dist - list of orders with distance matrices

get_matching_len_equ <- function(target.name, life.stage, id_info, equ_len,
                                 data.base = "itis", max_tax_dist, d.dist,
                                 len_equ_data, length_only = FALSE) {
  
  # get a set of equations within the suitable maximum taxonomic distance
  df <- get_tax_distance(target.name = target.name, 
                         id_info = id_info, 
                         equ_len = equ_len, 
                         length_only = length_only,
                         data.base = data.base, max_tax_dist = max_tax_dist,
                         d.dist = d.dist)
  
  # if there is no suitable equation within the maximum taxonomic distance, then we return an NA
  if (is.na( unlist(df)[1] ) & (length(df) == 1 )) {
    return(NA)
  }
  
  # join these id's to the length data
  df <- left_join(as_tibble(df), len_equ_data, by = "id")
  message("Suitable taxonomic distances:")
  print(df)
  
  # determine if the life-stages match the chosen equations
  
  # test if the life stages are explicitly different
  if ( all(!is.na(df$life_stage)) & is.na(life.stage) | # if all are not NA and life.stage is NA
       !is.na(life.stage) & all(is.na(df$life_stage)) | # if life.stage is not NA and all life.stages are NA
        !is.na(life.stage) & all(life.stage != df$life_stage[!is.na(df$life_stage)] )  ) { # if life.stage is not NA and not equal to any not NA life stages
    
    message(paste("No suitable", equ_len,  ": Life stages differ") )
    message("Returning an NA")
    return(NA)
    
  }
  
  # check if the life stages match and select that equation
  x <- (df$life_stage == life.stage) & !is.na(df$life_stage)
  if ( any( !is.na(x) &  (x == TRUE ) ) ) {
    
    df.out <- df[x,]
    
  # if no direct matching life stages but all are NA, then choose the highest ranking taxon  
  } else if (any(is.na(df$life_stage)) & is.na(life.stage)) {
    
    # get the entries where the life.stage is NA
    df <- df[is.na(df$life_stage),]
    
    # set the taxonomic hierarchy
    tax.hier <- c("order", "suborder", "infraorder", "section", "subsection", "superfamily",
                  "family", "subfamily", "tribe", "subtribe", "genus", "species")
    
    # get the highest ranking equation
    rank.in <- tax.hier[ min( which(tax.hier %in% df$rank ) ) ]
    
    # output the full table and the length
    df.out <- df[df$rank == rank.in, ]
    
    # if the taxonomic ranks are equal, then we use them all 
  } else {
    
    df.out <- df
    
  }
  
  message("Suitable taxonomic distances and life stages:")
  print(df.out)
  return(df.out)
  
}

## Production function

# write a function to get mass output from lengthd and a given taxon name
get_mass_from_length <- function(target.name, 
                                 target.length,
                                 life.stage,
                                 data.base = "itis", 
                                 max_tax_dist = 6, 
                                 length_only = TRUE,
                                 default_length = FALSE,
                                 output = "algorithmic") {
  
  # set the seed
  set.seed(54728587)
  
  # test the if the length_only argument has been specified
  if (is.na(length_only)) {
    stop("Specify if length_only = TRUE/FALSE")
  }
  
  # get the suitable equations
  equ.out <- get_matching_len_equ(target.name = target.name, 
                                  life.stage = life.stage, 
                                  id_info = e_ti, equ_len = "equation",
                                  data.base = data.base, max_tax_dist = max_tax_dist, 
                                  d.dist = d.dist,
                                  len_equ_data = equ_id$equation_data, 
                                  length_only = length_only)
  
  # if there is no suitable equation within the maximum taxonomic distance, then we return an NA
  if (is.na( unlist(equ.out)[1] ) & (length(equ.out) == 1 )) {
    return(list("mass_data" = NA,
                "equation_data" = NA))
  }
  
  # if length_only = FALSE is specified and output is FALSE
  if (length_only == FALSE) {
    
    # join the suitable equations with the variable input data
    equ_join <- left_join(equ.out, equ_id$variable_input_data[, -c(7,8)], 
                          by = c("id", "db_taxon", "life_stage"))
    
    # test if all the equations are body_length equations anyway
    if( all(equ_join$size_measurement == "body_length")) {
      message("All suitable equations are length_only = TRUE, default lengths can be used or target.lengths can be specified")
    }
    
    # return the equ_join data.frame as an ouput
    return(equ_join)
    
  }
  
  # decide whether we want the function to get default length data
  if (default_length) {
    
    len.out <- get_matching_len_equ(target.name = target.name, life.stage = life.stage, 
                                    id_info = l_ti, equ_len = "length",
                                    data.base = data.base, max_tax_dist = 4, 
                                    d.dist = d.dist,
                                    len_equ_data = len_id)
    
    # if the len.out output is NA, then we stop the function
    if(is.na( unlist(len.out)[1] ) & (length(len.out) == 1 )) {
      return(NA)
    }
    
    # calculate the mean, sd and n of all individual entries left over
    len.in <- c("mean" = mean(len.out$length_mid_mm),
                "sd" = sd(len.out$length_mid_mm),
                "n" = length(len.out$length_mid_mm))
    
    len.x <- list(len.out, len.in)
    
    # use the mean of these multiple suitable length values
    message(paste("Using the mean length: ", round(len.in[1], 2), " mm ", 
                  "(sd: ", round(len.in[2], 2), ", n: ", len.in[3], " )", " to calculate mass", sep = "" ))
    len.use <- len.in[1]
    names(len.use) <- NULL
    
  } else if ( (default_length == FALSE ) & all(!is.na(target.length)) ) { 
    
    len.use <- target.length
    
  } else if( (default_length == FALSE ) & any(is.na(target.length)) ) {
    
    stop("default_length = FALSE and target.length is NA: Specify a target.length(s) or set default_length = TRUE")
    
  }
  
  # get the mass from the length data for each equation in the database
  mass_list <- 
    lapply(equ.out$id, function(x) {
      
      # get the suitable equation data
      equ.dat <- equ.out[equ.out$id == x,]  
      
      # match the variable input data to the equation(s)
      var.dat <- equ_id$variable_input_data[equ_id$variable_input_data$id == x, ]
      
      # calculate the mass from the length
      
      # assign the length(s) to the variable
      assign(x = var.dat[["variable"]], value = len.use)
      
      # implement the equation
      parsed_eq <- parse(text = equ.out[equ.out$id == x,][["equation"]] )
      
      # combine mass data with useful other data
      mass_df <- tibble(target_name = target.name,
                        target_life_stage = life.stage,
                        id = x,
                        size = len.use,
                        size_measurement = var.dat$size_measurement,
                        mass = eval(parsed_eq) )
      
      # generate a clean output data.frame
      full_join(equ.dat,
                mass_df, 
                by = "id" ) %>%
        select(target_name, target_life_stage, id, db_taxon, rank, life_stage, dist_to_target,
               size, size_measurement, mass, dw_ww, unit)
      
    } )
  
  # make sure that only one output is specified
  if ( !(output %in% c("algorithmic", "full")) )  {
    stop("Choose appropriate output i.e. algorithmic or full")
  }
  
  # if we choose the length algorithmically
  if (output == "algorithmic") {
    mass_out <- mass_list[sample(1:length(mass_list), 1)][[1]]
  }
  
  # if we choose the full output, then output this
  if (output == "full") {
    mass_out <- bind_rows(mass_list)
    mass_out <- list("equation_data" = equ.out,
                     "mass_data" = mass_out)
  }
  
  # if we used the default length algorithm then we add this to the output
  if (default_length) {
    mass_out <- list("mass_data" = mass_out,
                     "default_length_data" = len.x[[1]])
  }
  
  return(mass_out)
  
}

# test the function
x <- get_mass_from_length(target.name = "Dugesia aenigma", 
                     target.length = rnorm(n = 5, mean = 10, sd = 1),
                     life.stage = NA,
                     data.base = "itis",
                     max_tax_dist = 10, 
                     length_only = TRUE,
                     default_length = FALSE,
                     output = "full")
x
View(x$equation_data)
View(x$mass_data)

### END
