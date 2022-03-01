
# Search taxon database

# Next steps:

# implement some maximum taxonomic distance cut-off
# use the best equation to calculate the mass from some given length

# load relevant libraries
library(here)
library(dplyr)

# read in the equation data
equ.dat <- readxl::read_xlsx(here("raw_data/equation_data.xlsx"))
equ.dat <- equ.dat[!is.na(equ.dat$equation_id),]

# read in the variable inputs
in.dat <- readxl::read_xlsx(here("raw_data/variable_input_data.xlsx"))
in.dat <- in.dat[!is.na(in.dat$equation_id),]

# read in the taxonomic identifiers
td.dat <- readRDS(file = here("database/itis_taxon_identifiers.rds"))

# read in the taxonomic information for each order
td.dist <- readRDS(file = here("database/itis_order_taxon_information.rds"))

# list of equation IDs with only length data
x <- aggregate(in.dat$equation_id, by = list(in.dat$equation_id), length, simplify = TRUE)
y <- in.dat[in.dat$equation_id %in% x[x$x == 1, ]$Group.1, ]
id.length <- y[y$size_measurement == "body_length", ]$equation_id


# function to search database for best equations

# args
target.name <- "Ptilonyssus"
length.only <- TRUE
data.base <- "itis"

# load the relevant functions
source(here("scripts/functions/08_get_taxonomic_information_function.R"))

# if the input name is a species then extract the genus
search.name <- extract_genus(binomial = target.name)

# if length.only = TRUE then subset equations with only length data
if (length.only == TRUE) {
  x <- lapply(td.dat, function(x) x$equation_id %in% id.length)
  td <- td.dat[unlist(x)]
} else {td <- td.dat}

# get the order of the target
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

# select the correct taxonomic distance matrix
# select the correct distance matrix
order.dist <- td.dist[unlist(lapply(td.dist, function(x) x$order == search.order))][[1]]
dist.m <- order.dist$tax_distance
tax.names <- order.dist$tax_names
rm(order.dist)

equ.dist <- 
  lapply(td, function(x) { 
    
  if (x$equation_name == target.name) {
    d <- 0
  } else if ( (search.name %in% tax.names$taxonname) | (search.name %in% x$synonymns) ) {
    z <- extract_genus( x$equation_name)
    d1 <- dist.m[which(row.names(dist.m) == z), which(colnames(dist.m) == search.name) ]
    d2 <- dist.m[which(row.names(dist.m) == search.name), which(colnames(dist.m) == z ) ]
    d <- max(c(d1, d2))
    if (x$equation_rank == "species") {d <- d+1} # add species level distance
  } else {d <- NA}
  
  return(d)
  
  } )
rm(dist.m)

# unlist the distances
equ.dist <- unlist(equ.dist)
equ.min <- which(near(min(equ.dist), equ.dist))

# get the taxonomic rank table
tax.hier <- c("order", "suborder", "infraorder", "section", "subsection", "superfamily",
              "family", "subfamily", "tribe", "subtribe", "genus")
if (length(equ.min) > 1 ) {
  equ.ranks <- lapply(td[equ.min], function(x) {x$equation_rank} )
  x <- which(tax.hier %in% unlist(equ.ranks) )
  y <- which(x == min(x))
  best.equ <- td[equ.min[y]][[1]]$equation_id 
  } else {
    best.equ <- td[equ.min][[1]]$equation_id
  }
print(best.equ)

