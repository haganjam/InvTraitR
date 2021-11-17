
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# next steps:

# run all of these steps for each species in our equation database for
# all three databases
# at each run, we must add the equation id
# thus, when there is a match, we can immediately get to the relevant dataset
# and calculate the taxonomic proximity

# load relevant libraries
library(taxize)
library(readr)
library(here)
library(dplyr)

rm(list = ls() )

# load important functions
source(here("scripts/functions/nearly_equal_function.R"))

# to do:

# before implementing this function, we can have a list of all names in the database
# this will then mean that before anyone puts a name in, we can test if it is in the database
# this means we don't have to search them all sequentially

# for implementing the functions with not just weights but also lengths
# ask for dataframes with species names and each measurement
# then the function can merge them so it fills in NA's when a measurement is not necessary

# load the equation and data input databases
equ.dat <- readxl::read_xlsx(here("raw_data/equation_data.xlsx"))
equ.dat <- equ.dat[!is.na(equ.dat$equation_id),]
in.dat <- readxl::read_xlsx(here("raw_data/variable_input_data.xlsx"))
in.dat <- in.dat[!is.na(in.dat$equation_id),]

# taxon list
tax.list <- equ.dat[, c("equation_id", "equation_target_taxon") ]
tax.list <- split(tax.list, tax.list$equation_id)

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


# we will use three databases
data.base <- "itis"

# arguments for the function
rank.diff <- 1

# set a name
# x.name <- "Mesostigmata"
# x.name <- "Collembola" # absent from the itis database
x.name <- "Lumbriculidae" # present in the itis database

ask_or_not <- FALSE

# get the taxon id from the database
if (data.base == "bold") {
  
  x.taxa.a <- get_boldid(sci = x.name, ask = ask_or_not )
  
} else if (data.base == "itis") {
  
  x.taxa.a <- get_tsn(sci_com = x.name, ask = ask_or_not )
  
} else if (data.base == "gbif") {
  
  x.taxa.a <- get_gbifid(sci = x.name, ask = ask_or_not )
  
}


# get taxonomic information
if (attr(x.taxa.a, "match") == "found") {
  
  x.taxa <- x.taxa.a[[1]]
  
  # upwards classification
  # get the upwards classification
  x.class <- classification(sci_id = x.taxa, db = data.base)
  x.class <- x.class[[1]]
  
  # get the taxonomic rank
  x.rank <- x.class[x.class$name == x.name, "rank"]
  
} else { x.rank <- "a" }

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


# now we fork the x.df to create the null database if the taxa is not found
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


# add the equation label
equ.id <- equ.dat[equ.dat$equation_target_taxon == x.name, ]$equation_id

# get potential synonymns
equ.syn <- synonyms(x.taxa.a[[1]], db = data.base)

x.list <- 
  list(database = data.base,
       equation_id = ifelse(equ.id == 0, NA, equ.id),
       synonymns = ifelse(is.na(equ.syn[[1]]$syn_name), NA, equ.syn[[1]]$syn_name),
       taxonomic_information = x.df)

### END
