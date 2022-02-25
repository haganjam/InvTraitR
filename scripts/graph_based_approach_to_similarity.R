
# Graph-based approach to similarity

# next steps:

# add multiple tries to the downstream functions (see taxon_query_functions)

# this approach gets me a distance matrix and list of taxa to genus level

# code algorithmic how to incorporate species-level data into this (i.e. one step down from genus)

# package this into a function and get a distance matrix and list of taxa to genus level
# for each taxa in the equation database

# finally, this just needs to be searched...

# check for the correct packages:
# add a check based on the current version
source("scripts/functions/version_package_warnings.R")

# load relevant libraries
library(igraph)
library(taxize)
library(here)
library(dplyr)

# load relevant functions
source(here("scripts/functions/get_gbifid2_function.R"))
source(here("scripts/functions/gbif_downstream_function.R"))
source(here("scripts/functions/itis_downstream_function.R"))
source(here("scripts/functions/get_taxon_id.R"))


# example code

# set the input name
equ.name.input <- "Sinantherina socialis"
z <- unlist( strsplit(x = equ.name.input, split = " ", fixed = TRUE) )

if (length(z) > 1) {
  equ.name <- z[1]
}

# get taxon id: database specific function i.e. get taxon ID from the database
#
y <- get_taxon_id(database_function = "itis", taxon_name = equ.name, ask_or_not = FALSE, tries = 5)
#

# get upwards classification: general function
#
y.c <- classification(y)[[1]]
y.c$ranknum <- 1:nrow(y.c)

# get target rank
equ.rank <- y.c[y.c$name == equ.name,][["ranknum"]]

# get target rank name
if (length(z) > 1) {
  equ.name.name <- "species"
} else { equ.rank.name <- y.c[y.c$name == equ.name,][["rank"]] }

# get order rank
ord.rank <- y.c[y.c$rank == "order", ][["ranknum"]]
ord.name <- y.c[y.c$rank == "order", ][["name"]]

if (equ.rank < ord.rank) {
  
  message("equation above the rank of order")
  equ.suitable <- FALSE
  
} else { equ.suitable <- TRUE }
#

### stop process if the equation is not suitable
#

# get taxon id: database specific function i.e. get the taxon ID of the order of the taxa from the database
ord.id <- get_taxon_id(database_function = "itis", 
                        taxon_name =  ord.name, 
                        ask_or_not = FALSE, tries = 5)
ord.id <- ord.id[[1]]

data_base = "itis"
# get a distance matrix and list of taxa downstream of the order
if (data_base == "gbif") {
  
  dmat.taxlist <- downstream_gbif(ord.id = ord.id, ord.name = ord.name)
  
} else if (data_base == "itis") {
  
  dmat.taxlist <- downstream_itis(ord.id = ord.id, ord.name = ord.name)
  
  # add synonymns here
  equ.syn <- taxize::synonyms(y, db = data_base)
  
  # combine into dmat.taxlist
  dmat.taxlist <- c(dmat.taxlist, synonymns = equ.syn[[1]]$syn_name)
  
} else {
  
  stop("error, specify an appropriate database: itis or gbif")
  
}

# add the equation taxon name and the rank
dmat.taxlist <- c(equation_name = equ.name.input, 
                  equation_rank = equ.rank.name,
                  dmat.taxlist)



# these functions give distance matrices to the genus level

# next step is to incorporate species level information
# to this...

# i.e. add one level of distance when it matches the genus

# this needs to be coupled to the equation database

# the search functions should be very simple

# but, will need some decisions e.g. if same distance but family then choose genus

# taxon name database
db.name <- "Pediciini"
db.name

# target name
tar.name <- "Limnophila costata"
tar.name

# this calculation across all the different distance matrices
d[which(row.names(d) == db.name), which(colnames(d) == tar.name) ]

# rank of all the different names



### END

