
# Get taxonomic information for each equation in the database

# load relevant functions
library(here)
source(here("scripts/create_database/06_get_taxonomic_information_function.R"))

# set database to search
db <- "itis"


## Equation database

# implement the function for all taxa in the equation database
equ_id <- readRDS(file = here("database/equation_vars_database.rds"))
equ_dat <- equ_id$equation_data

# get the order for each equation in the database
order.equ <- vector("list", length = nrow(equ.dat))
for (i in 1:nrow(equ.dat)) {
  
  order.equ[[i]] <- 
    get_taxon_order(name.input = equ.dat$db_taxon[i], 
                    id = equ.dat$id[i], 
                    data.base = db) 
  
}

# get a list of suitable equations
order.equ <- order.equ[ sapply(order.equ, function(x) ifelse(x[["suitable"]], TRUE, FALSE))]


## Default length database

# implement the function for all taxa in the default length database
len_dat <- readRDS(file = here("database/default_length_database.rds"))

# run the function for a few taxa on the list
len.test <- len_dat[sample(1:nrow(len_dat), 15), ]
print(len.test)

# get the order for each equation in the database
order.len <- vector("list", length = nrow(len.test))
for (i in 1:nrow(len.test)) {
  
  order.len[[i]] <- 
    get_taxon_order(name.input = len.test$db_taxon[i], 
                    id = len.test$id[i], 
                    data.base = db) 
  
}

# get a list of suitable default lengths
order.len <- order.len[ sapply(order.len, function(x) ifelse(x[["suitable"]], TRUE, FALSE))]


## Get taxonomic distances for each order
uni.order <- unique(c(sapply(order.len, function(x) x[["order"]]), 
                      sapply(order.equ, function(x) x[["order"]])
                      ))

# for each unique order, get the taxonomic distance matrix
dmat.order <- vector("list", length = length(uni.order))
for (j in 1:length(uni.order)) {
  
  u <- get_taxon_distance(ord.name = uni.order[j], data.base = db)
  
  dmat.order[[j]] <- u
  
}

# subset out the orders where taxonomic information is available for the order
x <- lapply(dmat.order, function(x) x$tax_available)
dmat.order <- dmat.order[unlist(x)]

# get the orders from the unique order list where taxonomic information is available
uni.order.in <- uni.order[unlist(x)]


## Equation data

# subset out any equations where the order information is not available
order.equ <- order.equ[ sapply(order.equ, function(x) ifelse(x[["order"]] %in% uni.order.in, TRUE, FALSE)) ]

# save an equation database with taxonomic identifiers
saveRDS(order.equ, file = here("database/itis_taxon_identifiers_equation.rds") )


## Default length data

# subset out any equations where the order information is not available
order.len <- order.len[ sapply(order.len, function(x) ifelse(x[["order"]] %in% uni.order.in, TRUE, FALSE)) ]

# save an equation database with taxonomic identifiers
saveRDS(order.len, file = here("database/itis_taxon_identifiers_length.rds") )


## Taxon order matrices

# save an equation database with the taxonomic information for each order
saveRDS(dmat.order, file = here("database/itis_order_taxon_information.rds") )

### END
