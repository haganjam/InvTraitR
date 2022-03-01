
# Search taxon database

# next steps:

# need to save memory:

# 1. use a sparse matrix to store the taxonomic information
# 2. do not store copies of each order's taxonomic information i.e. only unique orders
# - link the order information when required

# load relevant libraries
library(here)
library(dplyr)

# read in the equation data
equ.dat <- readxl::read_xlsx(here("raw_data/equation_data.xlsx"))
equ.dat <- equ.dat[!is.na(equ.dat$equation_id),]
View(equ.dat)

# read in the variable inputs
in.dat <- readxl::read_xlsx(here("raw_data/variable_input_data.xlsx"))
in.dat <- in.dat[!is.na(in.dat$equation_id),]
View(in.dat)

# read in the taxonomic data
td <- readRDS(file = here("database/itis_taxon_database.rds"))

# list of equation IDs with only length data
x <- aggregate(in.dat$equation_id, by = list(in.dat$equation_id), length, simplify = TRUE)
y <- in.dat[in.dat$equation_id %in% x[x$x == 1, ]$Group.1, ]
id.length <- y[y$size_measurement == "body_length", ]$equation_id


# function to search database for best equations

# args
target.name <- "Borboropsinae"
length.only <- TRUE
data.base <- "itis"

# if the input name is a species then extract the genus
z <- unlist( strsplit(x = target.name, split = " ", fixed = TRUE) )
if (length(z) > 1) {
  search.name <- z[1]
} else {search.name <- target.name}

# if length.only = TRUE then subset equations with only length data
if (length.only == TRUE) {
  x <- lapply(td, function(x) x$equation_id %in% id.length)
  td <- td[unlist(x)]
} 

# get the order of the target
source(here("scripts/functions/08_get_taxonomic_information_function.R"))
search.order <- get_taxon_order(equ.name.input = search.name, equ.id = NA, data.base = data.base, life.stage = NA)
search.order <- search.order$order

# get the equations with the correct order
y <- lapply(td, function(y) y$order == search.order)
y <- unlist(y)
if (sum(y) == 0) {
  equation_id <- NA
  stop("No suitable equation in the data for the target taxon name")
}
td <- td[y]

td.dist <- 
  lapply(td, function(x) { 
  
  if (x$equation_name == target.name) {
    d <- 0
  } else if ( (search.name %in% x$tax_names$taxonname) | (search.name %in% x$synonymns) ) {
    dist.m <- x$tax_distance
    d1 <- dist.m[which(row.names(dist.m) == x$equation_name), which(colnames(dist.m) == search.name) ]
    d2 <- dist.m[which(row.names(dist.m) == search.name), which(colnames(dist.m) == x$equation_name) ]
    d <- max(c(d1, d2))
  } else {d <- NA}
  
  return(d)
  
  } )

# unlist the distances
td.dist <- unlist(td.dist)
td.min <- which(near(min(td.dist), td.dist))

# get the taxonomic rank table
tax.hier <- c("order", "suborder", "infraorder", "section", "subsection", "superfamily",
              "family", "subfamily", "tribe", "subtribe", "genus")
if (length(td.min) > 1 ) {
  equ.ranks <- lapply(td[td.min], function(x) {x$equation_rank} )
  x <- which(tax.hier %in% unlist(equ.ranks) )
  y <- which(x == min(x))
  best.equ <- td[td.min[y]][[1]]$equation_id 
  } else {
    best.equ <- td[td.min]$equation_id
  }



# rank of all the different names
