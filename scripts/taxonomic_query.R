
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# load relevant libraries
library(taxize)
library(readr)
library(here)
library(dplyr)

rm(list = ls() )

# load important functions
source(here("scripts/functions/nearly_equal_function.R"))
source(here("scripts/functions/get_functions_error_handling.R"))

# load the equation and data input databases
equ.dat <- readxl::read_xlsx(here("raw_data/equation_data.xlsx"))
equ.dat <- equ.dat[!is.na(equ.dat$equation_id),]

# taxon list
tax.list <- equ.dat[, c("equation_id", "equation_target_taxon") ]
tax.list <- split(tax.list, tax.list$equation_id)

# function arguments
# we will use three databases
data.base <- "itis"

# arguments for the function
rank.diff <- 1

# set a name
# x.name <- "Mesostigmata"
# x.name <- "Collembola" # absent from the itis database
x.name <- "Lumbriculidae" # present in the itis database

ask_or_not <- FALSE

tries <- 5

# add the equation data as an argument
# add equation id and taxon name additionally


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
  x.class <- classification(sci_id = x.taxa, db = data.base)
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
      x.down <- downstream(sci_id = x.taxa, db = data.base, 
                           downto = down.to.rank, 
                           intermediate = TRUE)
      
      x.down <- x.down[[1]][[2]]
      
    }
    
    # get all species that share the rank with the focal species
    x.down.one.back <- downstream(sci_id = x.one.back$id, db = data.base,
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
    
    x.down <- bind_rows(x.down)
    x.down <- x.down[, c(5, 4, 6)]
    names(x.down) <- c("name", "rank", "id")
    
    x.down.one.back <- x.down.one.back[, c(5, 4, 6)]
    names(x.down.one.back) <- c("name", "rank", "id")
    
    x.down.merge <- rbind(x.down, x.down.one.back)
    
  } else if (data.base == "gbif") {
    
    x.down.merge <- rbind(bind_rows(x.down), x.down.one.back)
    x.down.merge$id <- paste(x.down.merge$key, x.down.merge$name_type, sep = ".")
    x.down.merge <- x.down.merge[, c("name", "rank", "id")]
    
  } else {
    
    x.down.merge <- rbind(bind_rows(x.down), x.down.one.back)
    
  }
  
  # merge the upstream and downstream species
  x.merge <- rbind(x.class[x.class$rank %in% x.up, ], x.down.merge )
  
  # merge with the df.tax
  x.df <- right_join(df.tax, x.merge, by = "rank")
  
  # add a focal taxa column
  x.df$focal_taxa <- ifelse(x.df$name == x.name, 1, 0) 
  
  # make sure that only the considered ranks are in
  x.df <- x.df[x.df$rank %in% df.tax$rank, ]
  
  # add the taxonomic database as a variable
  x.df$taxonomic_database <- data.base
  
}


# add relevant identification information and data are packaged into an output list

# add the equation label
equ.id <- equ.dat[equ.dat$equation_target_taxon == x.name, ]$equation_id

# get potential synonymns
equ.syn <- synonyms(x.taxa.a[[1]], db = data.base)

x.list <- 
  list(database = data.base,
       equation_id = ifelse(equ.id == 0, NA, equ.id),
       synonymns = ifelse(is.na(equ.syn[[1]]$syn_name), NA, equ.syn[[1]]$syn_name),
       taxonomic_information = x.df)

x.list

### END
