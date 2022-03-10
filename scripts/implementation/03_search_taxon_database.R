
# Search taxon database

# Next steps:

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
target.name <- "Gammarus roeseli"
length.only <- FALSE
data.base <- "itis"
max_tax_dist <- 6

# load the equation database
if (!exists("equ_id")) {
  equ_id <- readRDS(file = here("database/equation_vars_database.rds"))
}

# load the taxon information database
if (!exists("td")) {
  td <- readRDS(file = here("database/itis_taxon_identifiers.rds"))
}

# read in the taxonomic distance database
if (!exists("td.dist")) {
  td.dist <- readRDS(file = here("database/itis_order_taxon_information.rds"))
}

# load the relevant functions
source(here("scripts/functions/08_get_taxonomic_information_function.R"))

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

# if length.only = TRUE then subset equations with only length data
if (length.only == TRUE) {
  
  x <- lapply(td, function(x) x$equation_id %in% equ_id$id_only_equ_ID)
  td <- td[unlist(x)]
  
}

# subset the distance matrices that contain the taxon name
x <- lapply(td.dist, function(x) search.name %in% x$tax_names$taxonname)
y <- unlist(x)
td.dist.sel <- td.dist[y]

# select the correct taxonomic distance matrix
dist.m <- td.dist.sel[[1]]$tax_distance
print(paste("dimension: ", paste(dim(dist.m), collapse = " x ") ))

# get the list of taxonomic names
tax.names <- td.dist.sel[[1]]$tax_names$taxonname

# extract the order from the list
search.order <- unlist(lapply(td.dist.sel, function(x) x$order))[1]

# get the equations with the correct order
y <- lapply(td, function(y) y$order == search.order)
y <- unlist(y)

if (sum(y) == 0) {
  
  equation_id <- NA
  stop("No suitable equation in the data for the target taxon name")
  
}
td <- td[y]

# get the equation distance
equ.dist <- 
  lapply(td, function(x) { 
    
  if (x$equation_name == target.name) {
    
    d <- 0
    
    } else if ( (search.name %in% tax.names) | (search.name %in% x$synonymns) ) {
      
    z <- extract_genus( x$equation_name)
    d1 <- dist.m[which(row.names(dist.m) == z), which(colnames(dist.m) == search.name) ]
    d2 <- dist.m[which(row.names(dist.m) == search.name), which(colnames(dist.m) == z ) ]
    d3 <- max(c(d1, d2))
    
    # add extra distances for the species level
    if (x$equation_rank == "species") { d3 <- d3+0.25 } # add species level distance
    if (length(unlist( strsplit(x = target.name, split = " ", fixed = TRUE) )) > 1) { d3 <- d3+0.25 }
    
    } 
    
    else {d <- NA}
  
  return(d)
  
  } )
rm(dist.m)

# unlist the distances
equ.dist <- unlist(equ.dist)
equ.min <- which(near(min(equ.dist), equ.dist))
print(equ.dist[equ.min])

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
  
  equ.ranks <- lapply(td[equ.min], function(x) { x$equation_rank } )
  x <- which(tax.hier %in% unlist(equ.ranks) )
  y <- which(x == min(x))
  best.equ <- td[equ.min[y]][[1]]$equation_id 
  
  } else {
    
    best.equ <- td[equ.min][[1]]$equation_id
    
  }

# print the best.equation ID
print(best.equ)

### END
