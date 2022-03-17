
# Search taxon database

# Next steps:

# use the best equation to calculate the mass from some given length
# implement spread-sheet output

# build the default length data into this

# load relevant libraries
library(here)
library(dplyr)

# function to search database for best equations

# args
target.name <- "Toxomerini"
length.only <- FALSE
data.base <- "itis"
max_tax_dist <- 6
sim_level <- 1

# load the equation database
if (!exists("equ_id")) {
  equ_id <- readRDS(file = here("database/equation_vars_database.rds"))
}

# load the taxon information database
if (!exists("d.td")) {
  d.td <- readRDS(file = here("database/itis_taxon_identifiers.rds"))
}

# read in the taxonomic distance database
if (!exists("d.dist")) {
  d.dist <- readRDS(file = here("database/itis_order_taxon_information.rds"))
}

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

### EQUATION DATA

# if length.only = TRUE then subset equations with only length data
if (length.only) {
  
  d.td1 <- d.td[ sapply(d.td, function(x) x$equation_id %in% equ_id$id_only_equ_ID) ]
  
} else { d.td1 <- d.td }

# extract the orders present in the equation database
equ_orders <- sapply( d.td1, function(x) x$order)

# subset the distance matrices for the equations
equ_dist <- d.dist[ sapply(d.dist, function(x) x$order %in%  equ_orders) ]

# select the taxonomic distance matrix consistent with the taxonname
equ_dist <- equ_dist[ sapply(equ_dist, function(x) search.name %in% x$tax_names$taxonname) ]

# select the correct taxonomic distance matrix
equ_dist.m <- equ_dist[[1]]$tax_distance
print(paste("dimension: ", paste(dim(equ_dist.m), collapse = " x ") ))

# get the list of taxonomic names
equ_tax.names <- equ_dist[[1]]$tax_names$taxonname

# extract the order from the list
equ.order <- equ_dist[[1]]$order

# get the equations with the correct order
y <- sapply(d.td1, function(y) y$order == equ.order)
if (sum(y) == 0) {
  
  equation_id <- NA
  stop("No suitable equation in the data for the target taxon name")
  
}
d.td1 <- d.td1[y]


### DEFAULT LENGTH DATA


### WRITE INTO A FUNCTION to use on both the equations and default length data

# get the equation distance
equ_tax.dist <- 
  
  sapply(d.td1, function(x) { 
    
  if (x$name == target.name) {
    
    d <- 0
    
    } else if ( (search.name %in% equ_tax.names) | (search.name %in% x$synonymns) ) {
      
    z <- extract_genus( x$name)
    d1 <- equ_dist.m[which(row.names(equ_dist.m) == z), which(colnames(equ_dist.m) == search.name) ]
    d2 <- equ_dist.m[which(row.names(equ_dist.m) == search.name), which(colnames(equ_dist.m) == z ) ]
    d3 <- max(c(d1, d2))
    
    # add extra distances for the species level
    if (x$rank == "species") { d3 <- d3+0.25 } # add species level distance
    if (length(unlist( strsplit(x = target.name, split = " ", fixed = TRUE) )) > 1) { d3 <- d3+0.25 }
    
    } 
    
    else {d3 <- NA}
  
  return(d3)
  
  } )
rm(equ_dist.m)

# add equation range data i.e. min length individual and max length individual...

# choose output type i.e. single mass value or a spreadsheet with options?

equ.min <- which(near(min(equ_tax.dist), equ_tax.dist))
print(equ_tax.dist[equ.min])

# check if the minimum taxonomic is within the chosen threshold: max_tax_dist
if ( any(equ.min > max_tax_dist ) ) {
  
  equation_id <- NA
  stop(paste( "No suitable equation in the data for the target taxon name", "as the maximum taxonomic distance is greater than ", max_tax_dist, dep = "") )
  
}

# get the taxonomic rank table
tax.hier <- c("order", "suborder", "infraorder", "section", "subsection", "superfamily",
              "family", "subfamily", "tribe", "subtribe", "genus")

# if there are more than one suitable equation with the same minimum taxonomic distance
# choose the higher taxonomic level
if (length(equ.min) > 1 ) {
  
  equ.ranks <- sapply(td[equ.min], function(x) { x$rank } )
  x <- which(tax.hier %in% equ.ranks )
  y <- which(x == min(x))
  best.equ <- d.td1[equ.min[y]][[1]]$id 
  
  } else {
    
    best.equ <- d.td1[equ.min][[1]]$id
    
  }

# print the best.equation ID
print(best.equ)

### END
