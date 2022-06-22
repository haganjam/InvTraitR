
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
order.equ <- vector("list", length = nrow(equ_dat))
for (i in 1:nrow(equ_dat)) {
  
  order.equ[[i]] <- 
    get_taxon_order(name.input = equ_dat$db_taxon[i], 
                    id = equ_dat$id[i], 
                    life.stage = equ_dat$life_stage[i],
                    data.base = db) 
  
}

# get a list of suitable equations
order.equ <- order.equ[ sapply(order.equ, function(x) ifelse(x[["suitable"]], TRUE, FALSE))]


## Default length database

# implement the function for all taxa in the default length database
len_dat <- readRDS(file = here("database/default_length_database.rds"))
names(len_dat) <- c("id", "db_taxon", "Database", "life_stage", "length_mid_mm")

# get the order for each equation in the database
order.len <- vector("list", length = nrow(len_dat))
for (i in 1:nrow(len_dat)) {
  
  order.len[[i]] <- 
    get_taxon_order(name.input = len_dat$db_taxon[i], 
                    id = len_dat$id[i], 
                    life.stage = len_dat$life_stage[i],
                    data.base = db) 
  
}

# get a list of suitable default lengths
order.len <- order.len[ sapply(order.len, function(x) ifelse(x[["suitable"]], TRUE, FALSE))]


## Get taxonomic distances for each order
uni.order <- unique(c(sapply(order.len, function(x) x[["order"]]), 
                      sapply(order.equ, function(x) x[["order"]])
                      ))

# renew entire database?
renew <- FALSE
if (renew) {
  
  # for each unique order, get the taxonomic distance matrix
  dmat.order <- vector("list", length = length(uni.order))
  for (j in 1:length(uni.order)) {
    
    u <- get_taxon_distance(ord.name = uni.order[j], data.base = db)
    
    dmat.order[[j]] <- u
    
  }
  
}

# supplement the existing taxonomic distance matrix

# check which orders are already present in the taxonomic distance database
dmat.order <- readRDS(file = here("database/itis_order_taxon_information.rds"))
  
# remove any missing values
dmat.order <- dmat.order[sapply(dmat.order, function(x) !is.null(x$order) )]

# any of the unique orders are not in the d.dist matrix
if ( any( !(uni.order %in% sapply(dmat.order, function(x) x$order)) ) ) {
  
  mis.order <- uni.order[ which( !(uni.order %in% sapply(dmat.order, function(x) x$order)) ) ]
  
  # for each unique order, get the taxonomic distance matrix
  dmat.order.mis <- vector("list", length = length(mis.order))
  for (j in 1:length(mis.order)) {
    
    u <- get_taxon_distance(ord.name = mis.order[j], data.base = db)
    
    dmat.order.mis[[j]] <- u
    
  }
  
}

# join these databases
dmat.order <- c(dmat.order, dmat.order.mis)

# subset out the orders where taxonomic information is available for the order
x <- lapply(dmat.order, function(x) x$tax_available)
dmat.order <- dmat.order[unlist(x)]

# subset any orders where tax-distance information is missing
x <- lapply(dmat.order, function(x) identical(x$tax_distance, NA) )
dmat.order <- dmat.order[!unlist(x)]

# get the orders from the unique order list where taxonomic information is available
uni.order.in <- uni.order[uni.order %in% sapply(dmat.order, function(x) x$order)]
any( !(uni.order.in %in% sapply(dmat.order, function(x) x$order)) )


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
