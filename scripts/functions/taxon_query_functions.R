
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# Functions to extract taxonomic information for each equation in the database

# test for correct packages
if (any( !(c("taxize", "dplyr") %in% installed.packages()[,1]) )) {
  
  stop("error, this functions requires taxize and dplyr to be installed")
  
  warning("the function was written using taxize_0.9.99 and dplyr_1.0.2")
  
} else {
  
  warning("you have taxize and dplyr installed but, as a warning, this function was written using taxize_0.9.99 and dplyr_1.0.2")
  
}

message("this function was written using R version 4.0.2 (2020-06-22)")


# args:

# database_function: name of the function - "itis" (get_tsn), "bold" - (get_boldid), "gbif" - (get_gbifid)
# taxon_name: name that you want to query
# ask_or_not: accept imperfect matches (default is FALSE)
# tries: how many times to try the function if it keeps getting an error (default = 5)

get_taxon_id <- function(database_function = "itis", taxon_name, ask_or_not = FALSE, tries = 5) {
  
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
  
  library(taxize)
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


# Function that takes two vectors and determines if they are nearly equal etc.

# args

# x: vector of values
# y: vector of values
# tol: tolerance level for equality
# mode: "ae" almost equal, "gt" greater than etc.

# https://stackoverflow.com/questions/25154930/less-or-equal-for-floats-in-r
near_equal <- function( x , y , tol = 1.5e-8 , mode = "ae" ){
  
  ae <- mapply( function(x,y) isTRUE( all.equal( x , y , tolerance = tol ) ) , x , y )    
  gt <- x > y
  lt <- x < y
  if( mode == "ae" )
    return( ae )
  if( mode == "gt" )
    return( gt )
  if( mode == "lt" )
    return( lt )
  if( mode == "ne.gt" )
    return( ae | gt )
  if( mode == "ne.lt" )
    return( ae | lt )
}



# Function to extract upstream and downstream taxonomic information from a focal taxonomic name

# args

# x.name: taxon name to be queried
# data.base: database to query from ("itis", "bold", "gbif" are supported)
# rank.diff: acceptable difference in taxonomic level for a valid equation (e.g. species to genus = 1)
# tries: how many tries to query the online databases
# ask_or_not: accept direct matches only (i.e. FALSE which is the default)
# equ.dat: database containing the equation data
# taxon_var: variable name in equ.dat that contains the taxon name associated with the equation
# equation_id: variable name in equ.dat that contains the taxon name

Get_taxonomic_info <- function(x.name, 
                               equ.id = NA,
                               life_stage = NA,
                               data.base = "itis", 
                               rank.diff = 1, tries = 5, ask_or_not = FALSE,
                               equ.dat,
                               taxon_var = "equation_target_taxon",
                               equation_id = "equation_id") {
  
  if (is.na(equ_id)) {
    
    stop("this function requires an equation ID to match with the database")
    
  }
  
  # set important data for the function
  
  # set-up the taxonomic ranks that will be considered
  tax.ranks <- c('class','subclass', 
                 'superorder','order','suborder',
                 'superfamily','family', 'subfamily',
                 'genus', 'species') 
  
  # write this into a data.frame
  df.tax <- data.frame(rank = tax.ranks)
  
  # add a quantitative hierarchical taxonomic information
  rank_number <- c(1, (1+(1/3)), (1+(2/3)), 2, (2+(1/3)), (2+(2/3)), 3, (3+(1/3)), (3+(2/3)), (4+(2/3)))
  df.tax$rank_number <- rank_number
  
  # set the error message
  error_NA <- "NA due to ask=FALSE & no direct match found"
  
  
  # get the taxon ID from the relevant database
  
  # use the get_taxon_id function to obtain the taxon_id from the correct database
  # in addition, this function implements error handling when the databases are 
  # not accessed successfully
  x.taxa.a <- get_taxon_id(database_function = data.base, 
                           taxon_name = x.name, 
                           ask_or_not = FALSE, 
                           tries = 5) 
  
  
  # get the relevant taxonomic information from the chosen database
  
  # if(): if the taxon name is a direct match then we proceed to:
  # 1. get the upwards classification
  # 2. obtain the taxonomic rank of the focal taxon name
  
  # else(): if the taxon name is not a direct match then:
  # 1. we assign a taxon rank to "taxa not found"
  if (attr(x.taxa.a, "match") == "found") {
    
    x.taxa <- x.taxa.a[[1]]
    
    # upwards classification
    # get the upwards classification
    x.class <- taxize::classification(sci_id = x.taxa, db = data.base)
    x.class <- x.class[[1]]
    
    # get the taxonomic rank
    x.rank <- x.class[x.class$name == x.name, "rank"]
    
  } else { x.rank <- "taxa not found" }
  
  
  # if() the focal taxon name's rank is within the chose set of ranks (i.e. df. tax$rank) then:
  # 1. proceed to get all taxa down from and including the focal taxa
  if ( (x.rank %in% df.tax$rank) ) {
    
    # for this rank, we exact the numeric rank
    x.rank.num <- df.tax[df.tax$rank == x.rank, "rank_number"]
    
    # we take one rank back to get the downstream of the focal individal
    x.one.back <- x.class[which(x.class$rank == x.rank)-1,]
    
    # get all species up using the classification
    x.up <- df.tax[near_equal(rank_number, (x.rank.num - rank.diff), mode = "ne.gt") & near_equal(rank_number, x.rank.num, mode = "ne.lt"), "rank"]
    x.up <- x.up[x.up != x.rank]
    
    # downwards taxa
    # get the rank to get children
    down.to.rank <- df.tax[ near_equal(rank_number, x.rank.num, mode = "ne.gt") & near_equal(rank_number, (x.rank.num + rank.diff), mode = "ne.lt"), ]
    down.to.rank <- down.to.rank$rank[nrow(down.to.rank)]
    
    # this makes sure we don't duplicate if we only focus on the focal taxa
    if (down.to.rank == x.rank) {
      
      x.down <- NULL
      
    } else {
      
      # get all species down (excluding focal)
      x.down <- taxize::downstream(sci_id = x.taxa, db = data.base, 
                           downto = down.to.rank, 
                           intermediate = TRUE)
      
      x.down <- x.down[[1]][[2]]
      
    }
    
    # get all species that share the rank with the focal species
    x.down.one.back <- taxize::downstream(sci_id = x.one.back$id, db = data.base,
                                  downto = x.rank, intermediate = FALSE)
    x.down.one.back <- x.down.one.back[[1]]
    
  }
  
  
  # process the extracted taxonomic information
  
  # if(): the taxon name is not a match or the rank is not included in df.tax$rank then
  # 1. assign x.df i.e. taxonomic information to NULL
  
  # else(): the taxon name is a match and the rank is included in df.tax$rank then
  # 1. upstream and downstream taxanomic information are processed into x.df dataframes
  # 2. processing differs based on the output from the database
  if (attr(x.taxa.a, "match") %in% c("not found", error_NA) | !(x.rank %in% df.tax$rank) ) {
    
    x.df <- NULL
    
  } else {
    
    if (data.base == "itis") {
      
      x.down <- dplyr::bind_rows(x.down)
      x.down <- x.down[, c(5, 4, 6)]
      names(x.down) <- c("name", "rank", "id")
      
      x.down.one.back <- x.down.one.back[, c(5, 4, 6)]
      names(x.down.one.back) <- c("name", "rank", "id")
      
      x.down.merge <- rbind(x.down, x.down.one.back)
      
    } else if (data.base == "gbif") {
      
      x.down.merge <- rbind( dplyr::bind_rows(x.down), x.down.one.back)
      x.down.merge$id <- paste(x.down.merge$key, x.down.merge$name_type, sep = ".")
      x.down.merge <- x.down.merge[, c("name", "rank", "id")]
      
    } else {
      
      x.down.merge <- rbind(dplyr::bind_rows(x.down), x.down.one.back)
      
    }
    
    # merge the upstream and downstream species
    x.merge <- rbind(x.class[x.class$rank %in% x.up, ], x.down.merge )
    
    # merge with the df.tax
    x.df <- dplyr::right_join(df.tax, x.merge, by = "rank")
    
    # add a focal taxa column
    x.df$focal_taxa <- ifelse(x.df$name == x.name, 1, 0) 
    
    # make sure that only the considered ranks are in
    x.df <- x.df[x.df$rank %in% df.tax$rank, ]
    
    # add the taxonomic database as a variable
    x.df$taxonomic_database <- data.base
    
  }
  
  
  # add relevant identification information and data are packaged into an output list
  
  # get potential synonymns
  equ.syn <- taxize::synonyms(x.taxa.a[[1]], db = data.base)
  
  x.list <- 
    list(database = data.base,
         equation_id = ifelse(equ.id == 0, NA, equ.id),
         synonymns = ifelse(is.na(equ.syn[[1]]$syn_name), NA, equ.syn[[1]]$syn_name),
         life_stage = life_stage,
         taxonomic_information = x.df)
  
  return(x.list)
  
}

### END
