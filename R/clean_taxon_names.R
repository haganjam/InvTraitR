#' @title clean_taxon_names()
#' @description Get taxonomic distances of target names relative to the taxa
#'  databases
#' @details This function takes a data.frame with target_taxa and uses the bdc
#'  package to clean and harmonise the names.
#' @author James G. Hagan (james_hagan(at)outlook.com) and Ronald Bergmann
#' @param data - data.frame with a column containing target taxon names and
#'  life-stages
#' @param target_taxon - character string with the column name containing the
#'  taxon names
#' @param life_stage - character string with the column name containing the
#'  life_stages
#' @param database - taxonomic database to use: gbif (default), itis, col
#' @return data.frame with target taxon names and cleaned names for the chosen
#'  taxonomic backbone
#' @importFrom assertthat assert_that
#' @importFrom assertthat are_equal
#' @importFrom assertthat is.string
clean_taxon_names <- function(
    data,
    target_taxon,
    life_stage,
    database = "gbif") {
  # check that the database input is a supported taxonomic backbone
  assert_that(
    are_equal(database, "gbif") |
      are_equal(database, "itis") |
      are_equal(database, "col"),
    msg = paste(
      database,
      "is not a valid taxonomic backbone, pick: gbif, itis or col"
    )
  )
  
  # check that the data input is a data.frame or a tibble
  assert_that(
    is.data.frame(data) | dplyr::is.tbl(data),
    msg = paste(data, "is not a data.frame or tibble object")
  )
  
  # check that the target_taxon column is in the data object
  assert_that(
    is.string(target_taxon) & (target_taxon %in% names(data)),
    msg = paste(target_taxon, "is not a column in the supplied data object")
  )
  
  # check that the name column has a length of more than zero and that it is
  # a character variable
  data_target_taxon <- data[[target_taxon]]
  assert_that(
    is.character(data_target_taxon) & (length(data_target_taxon) > 0),
    msg = paste(
      data_target_taxon,
      "is not a character variable with length greater than zero"
    )
  )
  
  # check that the life-stage column is a character vector without NAs
  data_life_stage <- data[[life_stage]]
  assert_that(
    (is.character(data_life_stage) & all(data_life_stage %in% c(
      NA, "none", "larva", "pupa", "nymph",
      "adult", "nauplius", "copepodite"
    ))),
    msg = "one or more entries do not have appropriate life-stage classes: see documentation" # TODO: not yet in the docs
  )
  
  # add a targ_no column to the data input
  data[["row_id"]] <- 1:nrow(data)
  
  # clean the names for typos etc. using the bdc_clean_names function
  clean_names <- bdc::bdc_clean_names(
    sci_names = data[[target_taxon]],
    save_outputs = FALSE
  )
  
  # add the clean names to the data.frame
  clean_col <- paste("clean_", target_taxon, sep = "")
  data[[clean_col]] <- clean_names$names_clean
  
  # check that all the non-missing names were cleaned
  assert_that(
    length(data_target_taxon) == length(data[[clean_col]]),
    msg = "Length of clean names do not match length of original names"
  )
  
  # update the database if there is a valid internet connection
  if (curl::has_internet()) {
    taxadb::td_create(
      provider = database
    )
  }
  
  # add the database to the data
  data$db <- database
  
  # subset out taxa with special names
  spec_names <- special_taxon_names()
  data_spec <- 
    dplyr::filter(
      data,
      (eval(parse(text = clean_col)) %in% spec_names)
    )
  
  # change the database column to special
  if (nrow(data_spec) > 0) {
    data_spec[["db"]] <- "special"
  }
  
  # remove the special names from the data
  data <- dplyr::filter(data, !(row_id %in% data_spec[["row_id"]]))
  
  # if there are data points that are not special names
  # then we clean those names
  if (nrow(data) > 0) {
    # harmonise the names to the chosen data.base
    data_harm <- bdc::bdc_query_names_taxadb(
      sci_name = data[[clean_col]],
      db = database,
      rank_name = "Animalia",
      rank = "kingdom",
      export_accepted = FALSE
    )
    
    # add a row_id to this harm.tax object
    data_harm$row_id <- data$row_id
    
    data_harm <- 
      dplyr::rename(data_harm,
                    taxon_order = order,
                    taxon_family = family
      )
    
    # select the relevant columns
    data_harm <- data_harm[, c(
      "row_id",
      "scientificName",
      "taxonRank",
      "acceptedNameUsageID",
      "taxon_order",
      "taxon_family"
    )]
    
    # remove the names that we were not able to resolve
    data_harm <- dplyr::filter(
      data_harm,
      !(is.na(scientificName) | all( c(is.na(taxon_order), is.na(taxon_family)) ) )
    )
    
    # join these data to the tax.dat data
    data <- dplyr::left_join(data, data_harm, by = c("row_id"))
    
    # add the special names back
    data <- dplyr::bind_rows(data, data_spec)
    
  } else {
    # if we only have special names, then we only consider the special names and
    # add additional columns for consistency
    data <- data_spec
    data[["scientificName"]] <- NA
    data[["taxonRank"]] <- NA
    data[["acceptedNameUsageID"]] <- NA
    data[["taxon_order"]] <- NA
    data[["taxon_family"]] <- NA
  }
  
  # remove the row_id column
  data <- dplyr::select(data, -row_id)
  
  # convert to a tibble
  data <- dplyr::as_tibble(data)
  
  data
}
