
# Describe the characteristics of the database

# load the relevant libraries
library(here)
library(dplyr)
library(readr)
library(igraph)
library(ggplot2)

# load the taxon databases from the different taxonomic backbones
tax.dat <- 
  
  lapply(c("col", "gbif", "itis"), function(x) {
  
  y <- readRDS(file = paste0(here::here("database"), "/", x, "_taxon", "_database.rds"))
  return(y)
  
} )
tax.dat <- bind_rows(tax.dat)
View(tax.dat)

# select the relevant columns for describing the taxonomic distribution in the database
names(tax.dat)
tax.dat <- 
  tax.dat %>%
  select(group1, group2, db_higher_rank_source, id, scientificName, db_taxon_higher_rank, db_taxon_higher)
View(tax.dat)


