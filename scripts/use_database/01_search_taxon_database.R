
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





## Equation data

# if length.only = TRUE then subset equations with only length data
if (length.only) {
  
  e_ti1 <- e_ti[ sapply(e_ti, function(x) x$id %in% equ_id$id_only_equ_ID) ]
  
} else { e_ti1 <- e_ti }

# extract the orders present in the equation database
equ_orders <- sapply( e_ti1, function(x) x$order)

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
y <- sapply(e_ti1, function(y) y$order == equ.order)
if (sum(y) == 0) {
  
  equation_id <- NA
  stop("No suitable equation in the data for the target taxon name")
  
}
e_ti1 <- e_ti1[y]


## Default length data


### WRITE INTO A FUNCTION to use on both the equations and default length data

# get the equation distance
equ_tax.dist <- 
  
  sapply(e_ti1, function(x) { 
    
  if (x$name == target.name) {
    
    d <- 0
    
    } else if ( (search.name %in% equ_tax.names) ) {
      
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
  
  equ.ranks <- sapply(e_ti1[equ.min], function(x) { x$rank } )
  x <- which(tax.hier %in% equ.ranks )
  y <- which(x == min(x))
  best.equ <- e_ti1[equ.min[y]][[1]]$id 
  
  } else {
    
    best.equ <- e_ti1[equ.min][[1]]$id
    
  }

# print the best.equation ID
print(best.equ)

### END
