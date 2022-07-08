
# Search the database

# load relevant spatial libraries
library(sp)
library(sf)
library(raster)
library(dplyr)
library(here)
library(igraph)
library(bdc)
library(taxadb)

# load the taxon data
if (!exists("itis_td")) {
  itis_td <- readRDS(file = here("database/itis_taxon_database.rds"))
}

if (!exists("itis_htm")) {
  itis_htm <- readRDS(file = here("database/itis_higher_taxon_matrices.rds"))
}

if (!exists("gbif_td")) {
  gbif_td <- readRDS(file = here("database/gbif_taxon_database.rds"))
}

if (!exists("gbif_htm")) {
  gbif_htm <- readRDS(file = here("database/gbif_higher_taxon_matrices.rds"))
}

if (!exists("col_td")) {
  col_td <- readRDS(file = here("database/col_taxon_database.rds"))
}

if (!exists("col_htm")) {
  col_htm <- readRDS(file = here("database/col_higher_taxon_matrices.rds"))
}

# load the freshwater habitat map data
if (!exists("fw_er")) {
  fw_er <- readRDS(file = here("database/freshwater_ecoregion_data.rds"))
}

# load the freshwater habitat map
if (!exists("fw_map")) {
  fw_map <- readRDS(file = here("database/freshwater_ecoregion_map.rds"))
}

# load the equation data
if (!exists("equ.dat")) {
  equ.dat <- readRDS(file = here("database/equation_database.rds"))
}
head(equ.dat)

# test data.frame
t.dat <- data.frame(name = "Coloburiscus",
                    life_stage = "nymph",
                    length_mm = 10,
                    lat = c(-25.18),
                    lon = c(147.69))
t.dat

# choose the taxonomic database: "gbif", "itis", "col"
database <- "gbif"

# create the local database
td_create(
  provider = database,
  overwrite = FALSE)

# clean the names for typos etc.
x <- bdc_clean_names(sci_names = t.dat$name, save_outputs = FALSE)

# check if any names were changed
if ( !any(x$scientificName != x$names_clean) ) {
  message("No names were changed")
}

# replace the names in tax.dat with these cleaned names
t.dat$db_taxon <- x$names_clean

t.dat <- 
  t.dat %>%
  mutate(row_id = 1:n()) %>%
  dplyr::select(row_id, db_taxon, life_stage, length_mm, lat, lon)

# harmonise the names to the gbif database
harm.tax <- 
  bdc_query_names_taxadb(sci_name = t.dat$db_taxon,
                         db = database,
                         rank_name = "Animalia",
                         rank = "kingdom"
  )

# write some code to remove the output file
unlink("Output", recursive=TRUE)

harm.tax <- 
  harm.tax %>%
  mutate(db_taxon_higher_rank = ifelse(is.na(order) & is.na(family), NA, 
                                       ifelse(is.na(order) & !is.na(family), "family", "order") ) ) %>%
  mutate(db_taxon_higher = ifelse(is.na(order) & is.na(family), NA, 
                                  ifelse(is.na(order), family, order) ) ) %>%
  mutate(db_higher_rank_source = database) %>%
  mutate(row_id = 1:n()) %>%
  select(row_id, original_search, scientificName, acceptedNameUsageID, db_higher_rank_source, db_taxon_higher_rank, db_taxon_higher)

# remove the names that we were not able to resolve
harm.tax <- 
  harm.tax %>%
  filter(!(is.na(scientificName) |is.na(db_taxon_higher_rank) | is.na(db_taxon_higher) ) ) %>%
  rename(db_taxon = original_search)

# join these data to the tax.dat data
t.clean <- right_join(t.dat, harm.tax, by = c("row_id", "db_taxon") )

# check that the join worked correctly
if ( nrow(harm.tax) == nrow(t.clean) ) {
  message("Join worked correctly")
}

# remove the row_id column
t.clean <- 
  t.clean %>%
  select(-row_id)

gbif_htm

v.x <- V(d.g) 
v.x[which(attr(v.x, "names") == "Silvanidae")]

distances(d.g, 
          v.x[which(attr(v.x, "names") == "Silvanidae")],
          v.x[which(attr(v.x, "names") == "Anexantha")],
          mode = c("all"),
          algorithm = c("bellman-ford"))

