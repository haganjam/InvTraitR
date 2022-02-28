
# Get taxonomic information for each equation in the database

# load relevant functions
library(here)
source(here("scripts/functions/08_get_taxonomic_information_function.R"))

# implement the function for all taxa in the equation database
equ.dat <- readxl::read_xlsx(here("raw_data/equation_data.xlsx"))
equ.dat <- equ.dat[!is.na(equ.dat$equation_id),]

# taxon list
tax.list <- equ.dat[, c("equation_id", "equation_target_taxon", "life_stage")]

# run the function for a few taxa on the list
x.samp <- sample(1:nrow(tax.list), 12)
df.test <- tax.list[x.samp,]
print(df.test)

# create the itis database
db <- "itis"

# get the order for each equation in the database
order.info <- vector("list", length = nrow(df.test))
for (i in 1:nrow(df.test)) {
  v <- 
    get_taxon_order(equ.name.input = df.test$equation_target_taxon[i], 
                       equ.id = df.test$equation_id[i], 
                       data.base = db) 
  order.info[[i]] <- v
}

# get a list of suitable equations
x <- sapply(order.info, function(x) if (x[["equation_suitable"]] == TRUE) { TRUE } else { FALSE }  )
print(sum(x))
order.info <- order.info[x]

# get list of unique orders  
uni.order <- unique(sapply(order.info, function(x) x[["order"]]))

# for each unique order, get the taxonomic distance matrix
dmat.order <- vector("list", length = length(uni.order))
for (j in 1:length(uni.order)) {
  
  u <- get_taxon_distance(ord.name = uni.order[j], data.base = db)
  
  dmat.order[[j]] <- u
  
}

# list of orders in dmat
tax.database <- lapply(order.info, function(x) {
  
  y <- uni.order == x[["order"]]
  z <- c(x, dmat.order[y][[1]])
  return(z)
  
} )

saveRDS(tax.database, file = here("database/itis_taxon_database.rds") )

### END
