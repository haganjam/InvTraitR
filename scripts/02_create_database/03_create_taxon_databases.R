
#' @title itis, gbif and col databases
#' 
#' @description create itis, gbif and col taxon database and higher-level taxon matrices
#' 
#' @details This script uses the taxadb and bdc packages to clean the taxonomic
#' names in the taxon database, get the higher classification for each taxon in the
#' taxon database and generate higher-level taxon matrices for each taxon in the taxon
#' database. It supports three different taxonomic frameworks: itis, gbif and col.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 

# load relevant libraries
library(taxadb)
library(dplyr)
library(readr)
library(bdc)
library(tidyr)
library(Matrix)
library(here)

# check for the correct packages
source(here("scripts/02_create_database/01_version_package_warnings.R"))

# load the taxonomic distance matrix
source(here("scripts/02_create_database/02_taxon_matrix_function.R"))

# choose the taxonomic database: "gbif", "itis", "col"
database <- "gbif"

# create the local database
td_create(
  provider = database,
  overwrite = FALSE)

# load the taxon data
tax.dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/taxon_database.xlsx")
head(tax.dat)

# remove the empty columns
tax.dat <-
  tax.dat %>%
  select(-db_order_source, -db_taxon_higher, -db_taxon_higher_rank)

# add a row_id column
tax.dat <- 
  tax.dat %>%
  mutate(row_id = 1:n()) %>%
  select(row_id, group1, group2, db_taxon, db_taxon_gt_order)

# clean the names for typos etc.
x <- bdc_clean_names(sci_names = tax.dat$db_taxon, save_outputs = FALSE)

# check if any names were changed
if ( !any(x$scientificName != x$names_clean) ) {
  message("No names were changed")
}

# replace the names in tax.dat with these cleaned names
tax.dat$db_taxon <- x$names_clean

# harmonise the names to the gbif database
harm.tax <- 
  bdc_query_names_taxadb(sci_name = tax.dat$db_taxon,
                         db = database,
                         rank_name = "Animalia",
                         rank = "kingdom"
  )

# write some code to remove the output file
unlink("Output", recursive=TRUE)

# process the harmonised name taxa
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
tax.clean <- right_join(tax.dat, harm.tax, by = c("row_id", "db_taxon") )

# check that the join worked correctly
if ( nrow(harm.tax) == nrow(tax.clean) ) {
  message("Join worked correctly")
}

# remove the row_id column
tax.clean <- 
  tax.clean %>%
  select(-row_id)


# create the gbif taxon matrices

# get distinct higher taxa
d.ht <- 
  tax.clean %>%
  select(db_taxon_higher_rank, db_taxon_higher) %>%
  distinct()

d.dist <- vector("list", length = nrow(d.ht))
for (i in 6) {

  # get classification data for the higher taxon
  raw_class <- 
    filter_rank(name = d.ht[i, ]$db_taxon_higher, 
                rank =  d.ht[i, ]$db_taxon_higher_rank, 
                provider = database, 
                collect = FALSE
    ) %>%
    filter(!is.na(scientificName)) %>%
    filter(taxonomicStatus == "accepted") %>%
    select(order, family, genus) %>%
    distinct()
  
  # process data depending on whether the higher rank is order or family
  if (d.ht[i, ]$db_taxon_higher_rank == "order") {
    
    # some entries don't have proper classification data so we remove these
    raw_class <- raw_class[complete.cases(raw_class), ]
    
    proc_class <-
      
      bind_rows(
        
        raw_class %>%
          select(genus, family) %>%
          rename(name = genus, parentname = family) %>%
          mutate(rank = "genus") %>%
          mutate(parentrank = "family") %>%
          select(name, rank, parentname, parentrank), 
        
        raw_class %>%
          select(family, order) %>%
          rename(name = family, parentname = order) %>%
          mutate(rank = "family") %>%
          mutate(parentrank = "order") %>%
          select(name, rank, parentname, parentrank) 
        
      )
    
  } else if (d.ht[i, ]$db_taxon_higher_rank == "family") {
    
    raw_class <- 
      raw_class %>%
      select(-order)
    
    raw_class <- raw_class[complete.cases(raw_class), ]
    
    proc_class <- 
      raw_class %>%
      select(genus, family) %>%
      rename(name = genus, parentname = family) %>%
      mutate(rank = "genus") %>%
      mutate(parentrank = "family") %>%
      select(name, rank, parentname, parentrank)
    
  }
  
  # apply taxonomic weights
  weights <- mapply(function(x, y) { 
    
    tax.d[which(row.names(tax.d) == x), which(colnames(tax.d) == y) ] 
    
  } ,
  x = proc_class$rank,
  proc_class$parentrank)
  
  # add weights to the processed classification data
  proc_class$weights <- unlist(weights, use.names = FALSE)
  
  # create the distance matrix
  d.mat <- 
    proc_class %>%
    select(from = parentname, to = name, weights)
  
  # use igraph to create a graph from the matrix
  d.g <- graph_from_data_frame(d = d.mat, directed=FALSE)
  
  # write these sparse matrices into a list
  d.dist[[i]] <- d.g
  
}

# add the names of the higher taxa to the matrix
names(d.dist) <- d.ht$db_taxon_higher


# write these databases into the database folder

# write the taxon database
name1 <- paste(paste(database, "taxon", "database", sep = "_"), ".rds", sep = "")
saveRDS(tax.clean, file = paste(here("database"), "/", name1, sep = ""))

# write the higher taxon matrices
name2 <- paste(paste(database, "higher", "taxon", "matrices", sep = "_"), ".rds", sep = "")
saveRDS(tax.clean, file = paste(here("database"), "/", name2, sep = ""))

### END
