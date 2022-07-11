

# save the igraphs and search them directly
v.x <- V(d.g) 
v.x[which(attr(v.x, "names") == "Silvanidae")]

distances(d.g, 
          v.x[which(attr(v.x, "names") == "Silvanidae")],
          v.x[which(attr(v.x, "names") == "Anexantha")],
          mode = c("all"),
          algorithm = c("bellman-ford"))




# load relevant libraries
library(here)
library(dplyr)
library(igraph)
library(Matrix)

# check for the relevant libraries
source(here("scripts/create_database/01_version_package_warnings.R"))

# Load the different databases

# Equation data

# load the equation database
if (!exists("equ_id")) {
  equ_id <- readRDS(file = here("database/equation_vars_database.rds"))
}

# load the taxon information database from the equations
if (!exists("e_ti")) {
  e_ti <- readRDS(file = here("database/itis_taxon_identifiers_equation.rds"))
}

# Default length data
if (!exists("len_id")) {
  len_id <- readRDS(file = here("database/default_length_database.rds"))
}

# load the taxon information database from the equations
if (!exists("l_ti")) {
  l_ti <- readRDS(file = here("database/itis_taxon_identifiers_length.rds"))
}

# Taxonomic distance data 

# read in the taxonomic distance database
if (!exists("d.dist")) {
  d.dist <- readRDS(file = here("database/itis_order_taxon_information.rds"))
  
  # remove any missing values
  d.dist <- d.dist[sapply(d.dist, function(x) !is.null(x$order) )]
}

# Supplementary databases

# load the supplementary equation database for difficult to identify taxa
if (!exists("equ_id2")) {
  equ_id2 <- readRDS(file = here("database/equation_vars_supp_database.rds"))
}

# load the supplementary default length database for difficult to identify taxa
if (!exists("len_id2")) {
  len_id2 <- readRDS(file = here("database/default_length_supp_database.rds"))
}


#'
#' @title syn_correct()
#' 
#' @description Replace taxon name with accepted name if it is a synonymn
#' 
#' @details Sometimes, the taxon name that you input to get an equation or length for
#' is not an accepted name but rather a synonymn of some accepted name. This function
#' will determine whether the inputted name is a synonymn and then change it to the accepted
#' name in the "itis" database. Only the "itis" database is supported at the moment.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param target.name - name to test for synonymns
#' @param data.base - taxonomic database to use (only "itis" is currently supported)
#' 

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

#'
#' @title get_tax_distance()
#' 
#' @description Function to get equation/length id with the lowest taxonomic distance
#' 
#' @details This function takes a taxon name that the user is interested in and searches
#' either the length or equation database for the equation or length with the lowest taxonomic
#' distance from the target taxon name.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param target.name - name of the taxon to get the equation and length data for
#' @param id_info - database to search e_ti/l_ti
#' @param equ_len - "equation" or "length" data
#' @param length_only - if TRUE, then we only consider equations that are based on body length
#' @param data.base - taxonomic database (only "itis" is currently supported)
#' @param max_tax_dist - maximum acceptable taxonomic distance between target.name and equation/length data
#' @param d.dist - taxonomic distance matrices for the orders
#' 

get_tax_distance <- function(target.name, id_info, equ_len, length_only = TRUE,
                             data.base = "itis", max_tax_dist = 6,
                             d.dist) {
  
  # create an output data.frame
  dist.df <- data.frame(id = NA,
                        rank = NA,
                        dist_to_target = NA)
  
  # use the synonymn function here to output an object called "search.name"
  search.name <- syn_correct(target.name = target.name, data.base = data.base)
  
  # if length.only = TRUE then subset equations with only length data
  if (length_only & (equ_len == "equation") ) {
    
    id_info <- id_info[ sapply(id_info, function(x) x$id %in% equ_id$id_only_equ_ID) ]
    
  }
  
  # extract the orders present in the equation database
  orders <- sapply( id_info, function(x) x$order)
  
  # subset the distance matrices for the equations
  dist <- d.dist[ sapply(d.dist, function(x) x$order %in% unique(orders) ) ]
  
  # select the taxonomic distance matrix consistent with the taxonname
  u <- sapply(dist, function(x) search.name %in% x$tax_names$taxonname)
  dist <- dist[ u ]
  
  # check if the name is found in the order taxonomic matrices
  if (all(u == FALSE)) {
    message(paste("No suitable ", equ_len, " in the data for the target taxon name"))
    message("Returning an NA")
    return(dist.df)
  }
  
  # select the correct taxonomic distance matrix: there can only be one correct one
  dist.m <- dist[[1]]$tax_distance
  
  # get the list of taxonomic names
  tax.names <- dist[[1]]$tax_names$taxonname
  
  # extract the order from the list
  id.order <- dist[[1]]$order
  
  # get the equations with the correct order
  y <- sapply(id_info, function(y) y$order == id.order)
  if (sum(y) == 0) {
    
    message(paste("No suitable ", equ_len, " in the data for the target taxon name"))
    message("Returning an NA")
    return(dist.df)
    
  }
  id_info <- id_info[y]
  
  # get the distance from the suitable equations or lengths
  tax.dist <- 
    
    sapply(id_info, function(x) { 
      
      # check if the same and rank are present in the database
      if (any( is.na(c(x$name, x$rank)) ) | identical(x$rank, character(0)) ) {
        return(dist.df)
      }
      
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
        return(dist.df)
        
      }
      
      return(d3)
      
    } )
  
  # remove the dist.m object to save memory
  rm(dist.m)
  
  # get the id's of the lowest taxonomic distances
  id.min <- which( near(min(tax.dist, na.rm = TRUE), tax.dist) )
  
  # check if the minimum taxonomic is within the chosen threshold: max_tax_dist
  if ( any(tax.dist[id.min] > max_tax_dist ) ) {
    
    message(paste( paste("No suitable ", equ_len, " in the data for the target taxon name"), 
                   "as the maximum taxonomic distance is greater than ", max_tax_dist, sep = "") )
    message("Returning an NA")
    return(dist.df)
    
  }
  
  # update the output data.frame
  # create an output data.frame
  dist.df <- data.frame(id = sapply(id_info[id.min], function(x) x$id),
                        rank = sapply(id_info[id.min], function(x) { x$rank } ),
                        dist_to_target = tax.dist[id.min])
  
  return(dist.df)

  }

#'
#' @title get_matching_len_equ()
#' 
#' @description Function to select the best equation from the set of suitable equations
#' 
#' @details This function takes a taxon name that the user is interested in and searches
#' either the length or equation database for the equation or length that best matches the 
#' target taxon name. This is first based on extracting the database entries with the lowest
#' taxonomic entries from the target taxon name. Then the method chooses the equation based on whether
#' the life-stages match. Matching life-stages are directly chosen. If no life-stage information
#' is provided for the target or in the database and there are multiple database entries with
#' equal taxonomic distance, then the highest ranking database entry is chosen.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param target.name - name of the taxon to get the equation and length data for
#' @param life.stage - life stage of the target.name
#' @param id_info - database to search e_ti/l_ti
#' @param equ_len - equation or length data
#' @param length_only - if TRUE, then we only consider equations that are based on body length
#' @param data.base - taxonomic database (only "itis" is currently supported)
#' @param max_tax_dist - maximum acceptable taxonomic distance between target.name and equation/length data
#' @param d.dist - taxonomic distance matrices for the orders
#' @param len_equ_data - either equ_id$equation_data for equations or len_id for length data
#' @param life_stage_ignore - if TRUE, then all equations with suitable taxonomic distance are outputted
#' 

get_matching_len_equ <- function(target.name, life.stage, id_info, equ_len,
                                 data.base = "itis", max_tax_dist, d.dist,
                                 len_equ_data, length_only = FALSE,
                                 life_stage_ignore = FALSE) {
  
  # get a set of equations within the suitable maximum taxonomic distance
  df <- get_tax_distance(target.name = target.name, 
                         id_info = id_info, 
                         equ_len = equ_len, 
                         length_only = length_only,
                         data.base = data.base, max_tax_dist = max_tax_dist,
                         d.dist = d.dist)
  
  
  # if there is no suitable equation within the maximum taxonomic distance, then we return an NA
  if (any(is.na(unlist(df)))) {
    return(NA)
  }
  
  # join these id's to the length data
  df <- left_join(as_tibble(df), len_equ_data, by = "id")
  
  # do we care about life-stages?
  if (life_stage_ignore) {
    return(df)
  }
  
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
  
  # add search for the ecoregion
  
  return(df.out)
  
}

#'
#' @title get_mass_from_length()
#' 
#' @description Function to get mass from length for a given taxon name
#' 
#' @details This function takes a taxon name that the user is interested in and searches
#' either the length or equation database for the equation or length that best matches the 
#' target taxon name. This is first based on extracting the database entries with the lowest
#' taxonomic entries from the target taxon name. Then the method chooses the equation based on whether
#' the life-stages match. Matching life-stages are directly chosen. If no life-stage information
#' is provided for the target or in the database and there are multiple database entries with
#' equal taxonomic distance, then the highest ranking database entry is chosen.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param target.name - name of the taxon to get the equation and length data for
#' @param target.length - length (mm) of the target.name
#' @param life.stage - life stage of the target.name
#' @param life_stage_ignore - if TRUE, then all equations with suitable taxonomic distance are outputted
#' @param data.base - taxonomic database (only "itis" is currently supported)
#' @param max_tax_dist - maximum acceptable taxonomic distance between target.name and equation/length data
#' @param length_only - if TRUE, then we only consider equations that are based on body length
#' @param output - three options:
#' "algorithmic" - if there are multiple equivalent equations, choose one randomly
#' "full" - if there are multiple equivalent equations options, all are used and reported
#' 

# write a function to get mass output from length and a given taxon name
get_mass_from_length <- function(target.name, 
                                 target.length,
                                 life.stage,
                                 life_stage_ignore = FALSE,
                                 data.base = "itis", 
                                 max_tax_dist = 6, 
                                 length_only = TRUE,
                                 output = "algorithmic") {
  
  # set the seed
  set.seed(54587)
  
  # test the if the length_only argument has been specified
  if (is.na(length_only)) {
    stop("Specify if length_only = TRUE/FALSE")
  }
  
  if (any(is.na(target.length))) {
    default_length <- TRUE
  } else {
    default_length <- FALSE
  }
  
  # get the suitable equations
  equ.out <- get_matching_len_equ(target.name = target.name, 
                                  life.stage = life.stage, 
                                  id_info = e_ti, equ_len = "equation",
                                  data.base = data.base, max_tax_dist = max_tax_dist, 
                                  d.dist = d.dist,
                                  len_equ_data = equ_id$equation_data, 
                                  length_only = length_only,
                                  life_stage_ignore = FALSE)
  
  # if there is no suitable equation within the maximum taxonomic distance, then we return an NA
  if (identical(equ.out, NA)) {
    return(NA)
  }
  
  # if length_only = FALSE is specified and output is FALSE
  if (length_only == FALSE) {
    
    # join the suitable equations with the variable input data
    equ_join <- left_join(equ.out, equ_id$variable_input_data[, -c(7,8)], 
                          by = c("id", "db_taxon", "life_stage"))
    
    # return the equ_join data.frame as an ouput
    return(equ_join)
    
  }
  
  # decide whether we want the function to get default length data
  if (default_length) {
    
    len.out <- get_matching_len_equ(target.name = target.name, life.stage = life.stage, 
                                    id_info = l_ti, equ_len = "length",
                                    data.base = data.base, max_tax_dist = max_tax_dist, 
                                    d.dist = d.dist,
                                    len_equ_data = len_id)
    
    # if the len.out output is NA, then we stop the function
    if( identical(len.out, NA) ) {
      return(NA)
    }
    
    # take the mean if there are multiple measurements to be used to calculate length
    len.use <- mean(len.out$length_mid_mm, na.rm = TRUE)
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
                        length = len.use,
                        length_measurement = var.dat$size_measurement,
                        default_length = default_length,
                        mass = eval(parsed_eq) )
      
      # generate a clean output data.frame
      full_join(equ.dat,
                mass_df, 
                by = "id" ) %>%
        select(target_name, target_life_stage, id, db_taxon, rank, life_stage, dist_to_target,
               length, length_measurement, default_length, mass, dw_ww, unit)
      
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
    mass_out <- c("equation_data" = list(mass_out),
                  "default_length_data" = list(len.out) )
  }
  
  return(mass_out)
  
}

#'
#' @title get_taxa_info()
#' 
#' @description Function to get all relevant information from the database for a vector of names
#' 
#' @details This function takes a vector of taxa names and associated life-stages and returns
#' all information in the database that is within an appropriate taxonomic distance and if the
#' life-stage matches. Information in both the equation
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param target.name - a vector of taxon names to get the equation and length data for
#' @param life.stage - vector of life stages of the target.name
#' @param data.base - taxonomic database (only "itis" is currently supported)
#' @param max_tax_dist - maximum acceptable taxonomic distance between target.name and equation/length data
#' 
#' @return list with two elements:
#' 1. equation_info: all available equation information for the target.name
#' 2. length_info: all available length information for the target.name
#' 
#' @warnings warnings (> 1 result : no direct match found)
#' from taxize's get_ function when no name directly matches the input
#' 

get_taxa_info <- function(target.name,
                          life.stage = NA,
                          data.base = "itis",
                          max_tax_dist = 5) {
  
  # perform some checks on the 
  if (length(life.stage) != length(target.name)) {
    life.stage <- rep(NA, length(target.name))
    message("Life stages not provided for each target.name, using NA for all")
  }
  
  # check that the target.name vector is a character vector
  if(!is.character(target.name)) {
    stop("target.name argument must be a character vector")
  }
  
  # check the target.names don't contain weird characters
  if(any(grepl(pattern = "/.|,|_", x = target.name))) {
    stop("target.name argument contains characters other than spaces and letters")
  }
  
  # check that the length of the string is two words
  z <- sapply(target.name, function(x) length(unlist( strsplit(x = x, split = " ", fixed = TRUE) )) )
  if(any(z > 2)) {
    stop("target.name argument contains characters with more than two words")
  }
  
  # remove these broad taxa from the data.frame
  target.name <- target.name[!(target.name %in% equ_id2$equation_data$db_taxon)]
  
  equ.df <- vector("list", length = length(target.name))
  len.df <- vector("list", length = length(target.name))
  
  for(i in seq_along(target.name)) {
    
    # get all the equation data when not considering only body length data
    x <- get_mass_from_length(target.name = target.name[i], 
                              target.length = NA,
                              life.stage = life.stage[i],
                              data.base = data.base,
                              max_tax_dist = max_tax_dist, 
                              length_only = FALSE,
                              output = "full"
    )
    
    # get all default length data
    y <- get_mass_from_length(target.name = target.name[i], 
                              target.length = NA,
                              life.stage = life.stage[i],
                              data.base = data.base,
                              max_tax_dist = max_tax_dist, 
                              length_only = TRUE,
                              output = "full"
    )
    
    # if data were reported, then output into the list
    if(identical(x, NA)) {
      equ.df[[i]] <- tibble("target.name" = target.name[i])
    } else {
      equ.df[[i]] <- as_tibble(bind_cols(data.frame("target.name" = target.name[i] ),
                                         x))
    }
    
    # if data were reported, then output into the list
    if( !("default_length_data" %in% names(y)) ) {
      len.df[[i]] <- tibble("target.name" = target.name[i])
    } else {
      len.df[[i]] <- as_tibble(bind_cols(data.frame("target.name" = target.name[i] ),
                                         y$default_length_data))
    }
    
  }
  
  list.out <- list("equation_info" = bind_rows(equ.df),
                   "length_info" = bind_rows(len.df),
                   "supp_info" = c(equ_id2, list(len_id2)))
  
  return(list.out)
  
}

#'
#' @title get_taxa_mass()
#' 
#' @description Function to get mass from length data algorithmically
#' 
#' @details This function takes a data.frame with data on the taxa name, taxa life-stage
#' and length data and outputs relevant mass data based on algorithmically choosing the best
#' equation present in the database by matching taxonomic distance and life-stage. If length
#' data is not provided, we use default values that also match the target taxa by taxonomic
#' distance and life-stage. 
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param data.base - taxonomic database (only "itis" is currently supported)
#' @param max_tax_dist - maximum acceptable taxonomic distance between target.name and equation/length data
#' @param data - data.frame with taxa names, life-stages and lengthd ata
#' @param target.name.col - name of the column in 'data' containing the taxa names
#' @param life.stage.col - name of the column in 'data' containing the life-stages
#' @param length.col - name of the column in 'data' containing the length data
#' 
#' @return data.frame with mass data for each length and additional meta.data relating
#' to the equation used to determine the length and, if necessary, the length data
#' 
 
get_taxa_mass <- function(data.base = "itis",
                          max_tax_dist = 5,
                          data,
                          target.name.col,
                          life.stage.col,
                          length.col) {
  
  # copy the data.frame into a df object
  df <- data
  
  # extract a vector of names
  target.name <- df[[target.name.col]]
  
  # check that the target.name vector is a character vector
  if(!is.character(target.name)) {
    stop("target.name argument must be a character vector")
  }
  
  # check the target.names don't contain weird characters
  if(any(grepl(pattern = "/.|,|_", x = target.name))) {
    stop("target.name argument contains characters other than spaces and letters")
  }
  
  # check that the length of the string is two words
  z <- sapply(target.name, function(x) length(unlist( strsplit(x = x, split = " ", fixed = TRUE) )) )
  if(any(z > 2)) {
    stop("target.name argument contains characters with more than two words")
  }
  
  # make a separate data.frame for these taxa
  df2 <- df[(target.name %in% equ_id2$equation_data$db_taxon), ]
  
  # remove these broad taxa from the data.frame
  df <- df[!(target.name %in% equ_id2$equation_data$db_taxon), ]
  
  # split the unique taxa and lifestage combinations and if length data are present or not
  length_na <- sapply(df[[length.col]], function(x) ifelse(is.na(x), 1, 0) )
  dfx <- split(df, paste(df[[target.name.col]], df[[life.stage.col]], length_na, sep = "_"))
  
  # if a taxa has length NA then make it a single group
  dfx <- 
    lapply(dfx, function(x) { 
      z <- is.na(x[[length.col]] )
      if(sum(z) > 1) {
        x[1,] 
      } else {
        x
      } }
    )
  
  m.df <- vector("list", length = length(dfx))
  for(i in seq_along(dfx)) {
    
    x <- get_mass_from_length(target.name = dfx[[i]][[target.name.col]][1], 
                              target.length = dfx[[i]][[length.col]],
                              life.stage = dfx[[i]][[life.stage.col]][1],
                              data.base = data.base,
                              max_tax_dist = max_tax_dist, 
                              length_only = TRUE,
                              output = "algorithmic"
    )
    
    if (identical(x, NA)) {
      
      m.df[[i]] <- tibble(target_name = dfx[[i]][[target.name.col]][1],
                          target_life_stage = dfx[[i]][[life.stage.col]][1],
                          size = dfx[[i]][[length.col]][1])
      
    } else {
      
      if ( c("equation_data") %in% names(x) ) {
        y <- x$equation_data
        y$length_id <- paste(x$default_length_data$id, collapse = "_")
      } else {
        y <- x
        y$length_id <- NA
      }
      
      # reorganise the columns
      y <- 
        y %>%
        select(target_name, target_life_stage, id, db_taxon, rank, life_stage, dist_to_target,
               length_id, length_measurement, default_length, length, mass, dw_ww, unit)
      
      m.df[[i]] <- y
      
    }
    
  }
  
  if (nrow(df2) > 0) {
    message(paste(df2$taxon, collapse = ", "))
    message("The rank of these taxa is too high, see supplementary database for suggestions")
    message("The supplementary database is exported as a list element in the output: $supp_dat")
    
    return(list( mass_data = bind_rows(m.df),
                 supp_dat = c(equ_id2, list(len_id2))) )
    
  }
  
  return(bind_rows(m.df))
  
}

### END
