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
  taxadb::td_create(
    database[j]
    )

  # load the taxon data
  tax_dat <- readRDS(file = "database/taxon_database.rds")

  # remove the empty columns
  tax_dat <-
    tax_dat |>
    dplyr::select(-db_higher_rank_source, -db_taxon_higher_rank, -db_taxon_higher)

  # remove the special names
  spec_names <- special_taxon_names()

  tax_dat <-
    tax_dat |>
    dplyr::filter(!(db_taxon %in% spec_names))

  # add a row_id column
  tax_dat <-
    tax_dat |>
    dplyr::mutate(row_id = 1:n()) |>
    dplyr::select(row_id, database, id, group1, group2, db_taxon, db_taxon_gt_order)

  # clean the names for typos etc.
  x <- bdc::bdc_clean_names(sci_names = tax_dat$db_taxon, save_outputs = FALSE)

  # check if any names were changed
  if (!any(x$scientificName != x$names_clean)) {
    message("No names were changed")
  }

  # replace the names in tax_dat with these cleaned names
  tax_dat$db_taxon <- x$names_clean

  # harmonise the names to the database[j] backbone
  harm_tax <-
    bdc::bdc_query_names_taxadb(
      sci_name = tax_dat$db_taxon,
      db = database[j],
      rank_name = "Animalia",
      rank = "kingdom"
      )

  # process the harmonised name taxa
  harm_tax <-
    harm_tax |>
    dplyr::mutate(db_source = database[j]) |>
    dplyr::mutate(row_id = 1:n()) |>
    dplyr::select(row_id, original_search, scientificName, acceptedNameUsageID, db_source, order, family)

  # remove the names that we were not able to resolve
  harm_tax <-
    harm_tax |>
    dplyr::filter(!(is.na(scientificName) | is.na(family) )) |>
    dplyr::rename(db_taxon = original_search)

  # join these data to the tax_dat data
  tax_clean <- dplyr::right_join(tax_dat, harm_tax, by = c("row_id", "db_taxon"))

  # check that the join worked correctly
  if (nrow(harm_tax) == nrow(tax_clean)) {
    message("Join worked correctly")
  }

  # remove the row_id column
  tax_clean <-
    tax_clean |>
    dplyr::select(-row_id)

  # create the taxon matrices

  # get distinct higher taxa
  d_ht <-
    tax_clean |>
    dplyr::select(order, family) |>
    dplyr::distinct()
  
  higher_class <- vector("list", length = nrow(d_ht))
  for (i in 1:nrow(d_ht)) {
    # get classification data for the family
    raw_class <-
      taxadb::filter_rank(
        name = d_ht[i, ]$family,
        rank = "family",
        provider = database[j]
      )
    
    # clean the classification data from filter_rank
    raw_class <- 
      raw_class |>
      dplyr::filter(!is.na(scientificName)) |>
      dplyr::filter(scientificName != d_ht[i, ]$family) |>
      dplyr::filter(taxonomicStatus == "accepted") |>
      dplyr::select(order, family, genus) |>
      dplyr::distinct()
    
    higher_class[[i]] <- raw_class 
    
    } 
  
  # bind into a data.frame
  higher_class <- dplyr::bind_rows(higher_class)
  
  # split by order or family if order is missing
  higher_class <- split(higher_class, ifelse( is.na(higher_class$order), higher_class$family, higher_class$order ))
  
  # create the higher taxonomic distance matrices
  
  # create an output list
  d_dist <- vector("list", length = length(higher_class))
  
  for (k in 1:length(higher_class) ) {
    # initiate an input data.frame to process
    input_class <- higher_class[[k]]
    # process data depending on whether the higher rank is order or family
    if ( all(!is.na(input_class[["order"]]) ) ) {
      # some entries don't have proper classification data so we remove these
      input_class <- input_class[complete.cases(input_class), ]
      
      proc_class <-
        dplyr::bind_rows(
          input_class |>
            dplyr::select(genus, family) |>
            dplyr::rename(name = genus, parentname = family) |>
            dplyr::mutate(rank = "genus") |>
            dplyr::mutate(parentrank = "family") |>
            dplyr::select(name, rank, parentname, parentrank),
          
          input_class |>
            dplyr::select(family, order) |>
            dplyr::rename(name = family, parentname = order) |>
            dplyr::mutate(rank = "family") |>
            dplyr::mutate(parentrank = "order") |>
            dplyr::select(name, rank, parentname, parentrank)
        )
      
    } else if ( all(is.na(input_class[["order"]])) ) {
      
      input_class <-
        input_class %>%
        dplyr::select(-order)
      
      input_class <- input_class[complete.cases(input_class), ]
      
      proc_class <-
        input_class |>
        dplyr::select(genus, family) |>
        dplyr::rename(name = genus, parentname = family) |>
        dplyr::mutate(rank = "genus") |>
        dplyr::mutate(parentrank = "family") |>
        dplyr::select(name, rank, parentname, parentrank)
    }
    
    # apply taxonomic weights
    weights <- mapply(
      function(x, y) {
        tax_d[which(row.names(tax_d) == x), which(colnames(tax_d) == y)]
      },
      x = proc_class$rank,
      proc_class$parentrank
    )
    
    # add weights to the processed classification data
    proc_class$weights <- unlist(weights, use.names = FALSE)
    
    # create the distance matrix
    d_mat <-
      proc_class |>
      dplyr::select(from = parentname, to = name, weights)
    
    # use igraph to create a graph from the matrix
    d_g <- graph_from_data_frame(d = d_mat, directed = FALSE)
    
    # write these sparse matrices into a list
    d_dist[[k]] <- d_g
    
  }

  # add the names of the higher taxa to the matrix
  names(d_dist) <- names(higher_class)

  # write these databases into the database folder

  # write the taxon database
  name1 <- paste(paste(database[j], "taxon", "database", sep = "_"), ".rds", sep = "")
  saveRDS(tax_clean, file = paste("database", "/", name1, sep = ""))

  # write the higher taxon matrices
  name2 <- paste(paste(database[j], "higher", "taxon", "matrices", sep = "_"), ".rds", sep = "")
  saveRDS(d_dist, file = paste("database", "/", name2, sep = ""))

  }
