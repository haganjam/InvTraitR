
# Get taxonomic information for each equation in the database

# load relevant functions
library(here)
library(dplyr)
source(here("scripts/create_database/06_get_taxonomic_information_function.R"))

# load the taxon data
tax.dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/taxon_database.xlsx")
head(tax.dat)

# remove the missing columns
tax.dat <- 
  tax.dat %>%
  select(group1, group2, db_taxon, db_taxon_gt_order)

# use the distinct() function to get unique records
dim(tax.dat)
tax.dat <- distinct(tax.dat)
# View(tax.dat)

# convert character NAs to real NAs
tax.dat[tax.dat == "NA"] <- NA

# remove 'special' names which don't need to be searched
tax.spe <- 
  tax.dat %>%
  filter(db_taxon_gt_order == "yes")

tax.dat <- 
  tax.dat %>%
  filter(db_taxon_gt_order != "yes")

# get highest taxonomic level equal to or below order
gbif_order <- 
  lapply(tax.dat$db_taxon, function(x) {
  
  y <- 
    get_taxon_order(name.input = x, 
                    data.base = "gbif")
  y <- bind_cols(y)
  y <- y[, c(1, 2, 4, 5, 6, 7)]
  names(y) <- c("db_taxon", "db_taxon_accepted", "suitable", "db_order_source", "db_taxon_higher_rank", "db_taxon_higher")
  
  return(y)
  
} )

# bind into a data.frame
gbif_order <- bind_rows(gbif_order)
View(gbif_order)

# remove the non-suitable equations
gbif_order <- 
  gbif_order %>%
  filter(suitable == TRUE)

# remove the equations where the higher tax.rank is not family or order
gbif_order <- 
  gbif_order %>%
  filter(db_taxon_higher_rank %in% c("family", "order"))


# get the taxon distance matrices

# create a database with the unique ranks and names
db <- distinct(gbif_order[, c("db_taxon_higher_rank", "db_taxon_higher")])

gbif_dm <- 
  mapply(function(x, y) {
  
  z <- get_taxon_distance(higher.tax.name = x,
                          higher.tax.rank = y,
                          data.base = "gbif")
  
  return(z) }, 
  
  db$db_taxon_higher, 
  db$db_taxon_higher_rank, SIMPLIFY = FALSE )

gbif_dm[[1]]










