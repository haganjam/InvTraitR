
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

# test the get.order function
a <- 
  lapply(tax.dat$db_taxon, function(x) {
  
  y <- 
    get_taxon_order(name.input = x, 
                    data.base = "itis")
  y <- bind_cols(y)
  y <- y[, c(1, 2, 4, 5, 6, 7)]
  names(y) <- c("db_taxon", "db_taxon_accepted", "suitable", "db_order_source", "db_taxon_higher_rank", "db_taxon_higher")
  
  return(y)
  
} )

View(bind_rows(a))

taxon_id <- get_taxon_id(database_function = "gbif", 
             taxon_name = "Chilina", ask_or_not = FALSE, tries = 5)

# get the upwards classification for the taxon_id
y.c <- classification(taxon_id)[[1]]
y.c$ranknum <- 1:nrow(y.c)
y.c

z <- get_taxon_order(name.input = "Chilina patagonica", 
                data.base = "gbif")

z

tax.order <- 
  full_join(tax.dat, 
            bind_rows(x),
            by = "db_taxon")
View(tax.order)  

y <- get_taxon_distance(ord.name = x$order, data.base = "gbif")
y$order
y$tax_available
y$tax_distance
unique(y$tax_names$rank)



