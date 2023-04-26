#' @title itis, gbif and col databases
#' @description create itis, gbif and col taxon database and higher-level taxon
#'  matrices
#' @details This script uses the taxadb and bdc packages to clean the taxonomic
#'  names in the taxon database, get the higher classification for each taxon
#'  in the taxon database and generate higher-level taxon matrices for each
#'  taxon in the taxon database. It supports three different taxonomic
#'  frameworks: itis, gbif and col.
#' @author James G. Hagan (james_hagan(at)outlook.com)

# load relevant libraries
library(taxadb)
library(dplyr)
library(readr)
library(bdc)
library(tidyr)
library(Matrix)

# load the special names function
source("R/special_names.R")

# load the taxonomic distance matrix
source("companion_scripts/02_create_database/helper-taxon-matrix-function.R")

# set-up a vector of taxonomic databases: "gbif", "itis", "col"
database <- c("gbif", "itis", "col")

# loop over the three databases
for (j in 1:length(database)) {
  # create the local database
  td_create(
    database[j]
    )

  # load the taxon data
  tax.dat <- readRDS(file = "database/taxon_database.rds")

  # remove the empty columns
  tax.dat <-
    tax.dat %>%
    dplyr::select(-db_higher_rank_source, -db_taxon_higher_rank, -db_taxon_higher)

  # remove the special names
  spec.names <- special_taxon_names()

  tax.dat <-
    tax.dat %>%
    filter(!(db_taxon %in% spec.names))

  # add a row_id column
  tax.dat <-
    tax.dat %>%
    mutate(row_id = 1:n()) %>%
    dplyr::select(row_id, database, id, group1, group2, db_taxon, db_taxon_gt_order)

  # clean the names for typos etc.
  x <- bdc_clean_names(sci_names = tax.dat$db_taxon, save_outputs = FALSE)

  # check if any names were changed
  if (!any(x$scientificName != x$names_clean)) {
    message("No names were changed")
  }

  # replace the names in tax.dat with these cleaned names
  tax.dat$db_taxon <- x$names_clean

  # harmonise the names to the database[j] backbone
  harm.tax <-
    bdc_query_names_taxadb(
      sci_name = tax.dat$db_taxon,
      db = database[j],
      rank_name = "Animalia",
      rank = "kingdom"
    )

  # process the harmonised name taxa
  harm.tax <-
    harm.tax %>%
    mutate(db_source = database[j]) %>%
    mutate(row_id = 1:n()) %>%
    dplyr::select(row_id, original_search, scientificName, acceptedNameUsageID, db_source, order, family)

  # remove the names that we were not able to resolve
  harm.tax <-
    harm.tax %>%
    filter(!(is.na(scientificName) | is.na(family) )) %>%
    rename(db_taxon = original_search)

  # join these data to the tax.dat data
  tax.clean <- right_join(tax.dat, harm.tax, by = c("row_id", "db_taxon"))

  # check that the join worked correctly
  if (nrow(harm.tax) == nrow(tax.clean)) {
    message("Join worked correctly")
  }

  # remove the row_id column
  tax.clean <-
    tax.clean %>%
    dplyr::select(-row_id)

  # create the taxon matrices

  # get distinct higher taxa
  d.ht <-
    tax.clean %>%
    dplyr::select(order, family) %>%
    distinct()
  
  higher_class <- vector("list", length = nrow(d.ht))
  for (i in 1:nrow(d.ht)) {
    # get classification data for the family
    raw_class <-
      taxadb::filter_rank(
        name = d.ht[i, ]$family,
        rank = "family",
        provider = database[j]
      )
    
    # clean the classification data from filter_rank
    raw_class <- 
      raw_class %>%
      filter(!is.na(scientificName)) %>%
      filter(scientificName != d.ht[i, ]$family) %>%
      filter(taxonomicStatus == "accepted") %>%
      dplyr::select(order, family, genus) %>%
      distinct()
    
    higher_class[[i]] <- raw_class 
    
    } 
  
  # bind into a data.frame
  higher_class <- bind_rows(higher_class)
  
  # split by order or family if order is missing
  higher_class <- split(higher_class, ifelse( is.na(higher_class$order), higher_class$family, higher_class$order ))
  
  # create the higher taxonomic distance matrices
  
  # create an output list
  d.dist <- vector("list", length = length(higher_class))
  
  for (k in 1:length(higher_class) ) {
    # initiate an input data.frame to process
    input_class <- higher_class[[k]]
    # process data depending on whether the higher rank is order or family
    if ( all(!is.na(input_class[["order"]]) ) ) {
      # some entries don't have proper classification data so we remove these
      input_class <- input_class[complete.cases(input_class), ]
      
      proc_class <-
        bind_rows(
          input_class %>%
            dplyr::select(genus, family) %>%
            rename(name = genus, parentname = family) %>%
            mutate(rank = "genus") %>%
            mutate(parentrank = "family") %>%
            dplyr::select(name, rank, parentname, parentrank),
          
          input_class %>%
            dplyr::select(family, order) %>%
            rename(name = family, parentname = order) %>%
            mutate(rank = "family") %>%
            mutate(parentrank = "order") %>%
            dplyr::select(name, rank, parentname, parentrank)
        )
      
    } else if ( all(is.na(input_class[["order"]])) ) {
      
      input_class <-
        input_class %>%
        dplyr::select(-order)
      
      input_class <- input_class[complete.cases(input_class), ]
      
      proc_class <-
        input_class %>%
        dplyr::select(genus, family) %>%
        rename(name = genus, parentname = family) %>%
        mutate(rank = "genus") %>%
        mutate(parentrank = "family") %>%
        dplyr::select(name, rank, parentname, parentrank)
    }
    
    # apply taxonomic weights
    weights <- mapply(
      function(x, y) {
        tax.d[which(row.names(tax.d) == x), which(colnames(tax.d) == y)]
      },
      x = proc_class$rank,
      proc_class$parentrank
    )
    
    # add weights to the processed classification data
    proc_class$weights <- unlist(weights, use.names = FALSE)
    
    # create the distance matrix
    d.mat <-
      proc_class %>%
      dplyr::select(from = parentname, to = name, weights)
    
    # use igraph to create a graph from the matrix
    d.g <- graph_from_data_frame(d = d.mat, directed = FALSE)
    
    # write these sparse matrices into a list
    d.dist[[k]] <- d.g
    
  }

  # add the names of the higher taxa to the matrix
  names(d.dist) <- names(higher_class)

  # write these databases into the database folder

  # write the taxon database
  name1 <- paste(paste(database[j], "taxon", "database", sep = "_"), ".rds", sep = "")
  saveRDS(tax.clean, file = paste("database", "/", name1, sep = ""))

  # write the higher taxon matrices
  name2 <- paste(paste(database[j], "higher", "taxon", "matrices", sep = "_"), ".rds", sep = "")
  saveRDS(d.dist, file = paste("database", "/", name2, sep = ""))

  }
