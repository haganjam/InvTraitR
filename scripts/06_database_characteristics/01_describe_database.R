
# Describe the characteristics of the database

# load the relevant libraries
library(here)
library(dplyr)
library(readr)
library(igraph)
library(ggplot2)

# describe the taxonomic coverage

# need to decide on some higher classification scheme so that we can quantify the
# coverage... How detailed should it be?

# load the taxon databases from the different taxonomic backbones
tax.dat <- 
  
  lapply(c("col", "gbif", "itis"), function(x) {
  
  y <- readRDS(file = paste0(here::here("database"), "/", x, "_taxon", "_database.rds"))
  y[y$group1 == "Acari", ][["group1"]] <- "Collembola"
  y[y$group2 == "Cladocera", ][["group2"]] <- "Diplostraca"
  
  return(y)
  
} )
tax.dat <- bind_rows(tax.dat)
View(tax.dat)

# select the relevant columns for describing the taxonomic distribution in the database
names(tax.dat)
tax.dat <- 
  tax.dat %>%
  select(group1, group2, db_taxon, db_higher_rank_source, id, scientificName, family, db_taxon_higher_rank, db_taxon_higher)
View(tax.dat)

# summarise the taxon data
tax.dat %>%
  group_by(group1) %>%
  summarise(n_family = length(unique(family)),
            n_taxa = length(unique(scientificName)))

tax.dat %>%
  filter(db_taxon_higher_rank == "order") %>%
  group_by(group1) %>%
  summarise(n_order = length(unique(db_taxon_higher)))

tax.dat %>%
  filter(db_taxon_higher_rank == "order") %>%
  View()

tax.dat %>%
  filter(db_taxon_higher_rank == "order") %>%
  pull(db_taxon_higher) %>%
  unique()

tax.dat %>%
  filter(db_taxon_higher_rank == "order") %>%
  group_by(db_higher_rank_source, group1, group2, db_taxon_higher) %>%
  summarise(n_family = length(unique(family))) %>%
  View()

tax.dat %>%
  filter(group1 == "Acari") %>%
  View()


# describe the geographical coverage




