#' @title clean_taxon_names()
#' @description Get taxonomic distances of target names relative to the taxa
#'  databases
#' @details This function takes a data.frame with target_taxa and uses the bdc
#'  package to clean and harmonise the names.
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' @param data - data.frame with a column containing target taxon names and
#'  life-stages
#' @param target_taxon - character string with the column name containing the
#'  taxon names
#' @param life_stage - character string with the column name containing the
#'  life_stages
#' @param database - taxonomic database to use: gbif (default), itis, col
#' @return data.frame with target taxon names and cleaned names for the chosen
#'  taxonomic backbone
#' @import bdc
#' @import dplyr
#' @import taxadb
#' @import curl
#' @import assertthat
clean_taxon_names <- function(
    data,
    target_taxon,
    life_stage,
    database = "gbif") {
    # check that the database input is a supported taxonomic backbone
    assert_that(
        assertthat::are_equal(database, "gbif") |
            assertthat::are_equal(database, "itis") |
            assertthat::are_equal(database, "col"),
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
        assertthat::is.string(target_taxon) & (target_taxon %in% names(data)),
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
    # TODO: this doesn't seem correct, there's NA in the list
    data_life_stage <- data[[life_stage]]
    assert_that(
        (is.character(data_life_stage) & all(data_life_stage %in% c(
            NA, "none", "larva", "pupa", "nymph",
            "adult", "nauplius", "copepodite", "tadpole"
        ))),
        msg = paste(
            data_life_stage,
            "one or more entries do not have appropriate life-stage classes: see documentation" # TODO: not yet in the docs
        )
    )

    # add a targ_no column to the data input
    data[["row_id"]] <- 1:nrow(data)

    # clean the names for typos etc. using the bdc_clean_names function
    clean.names <- bdc::bdc_clean_names(
        sci_names = data[[target_taxon]],
        save_outputs = FALSE
    )

    # add the clean names to the data.frame
    clean.col <- paste("clean_", target_taxon, sep = "")
    data[[clean.col]] <- clean.names$names_clean

    # check that all the non-missing names were cleaned
    assert_that(
        length(data_target_taxon) == length(data[[clean.col]]),
        msg = "Length of clean names do not match length of original names"
    )

    # update the database if there is a valid internet connection
    if (curl::has_internet()) {
        taxadb::td_create(
            provider = database,
            overwrite = FALSE
        )
    }

    # add the database to the data
    data$db <- ifelse(is.na(data[[clean.col]]), NA, database)

    # subset out taxa with special names
    spec.names <- special_taxon_names()
    data.spec <- dplyr::filter(
        data,
        (eval(parse(text = clean.col)) %in% spec.names)
    )

    # change the database column to special
    if (nrow(data.spec) > 0) {
        data.spec[["db"]] <- "special"
    }

    # remove the special names from the data
    data <- dplyr::filter(data, !(row_id %in% data.spec[["row_id"]]))

    # if the are data points that are not special names
    # then we clean those names
    if (nrow(data) > 0) {
        # harmonise the names to the chosen data.base
        data.harm <- bdc::bdc_query_names_taxadb(
            sci_name = data[[clean.col]],
            db = database,
            rank_name = "Animalia",
            rank = "kingdom",
            export_accepted = FALSE
        )

        # write some code to remove the output file
        unlink("Output", recursive = TRUE)

        # add a row_id to this harm.tax object
        data.harm$row_id <- data$row_id

        # higher taxon rank
        data.harm[["db_taxon_higher_rank"]] <-
            ifelse(is.na(data.harm[["order"]]) & is.na(data.harm[["family"]]), NA,
                ifelse(is.na(data.harm[["order"]]) & !is.na(data.harm[["family"]]), "family", "order")
            )
        
        # higher taxon names
        x.rows <- ifelse(is.na(data.harm[["order"]]) & is.na(data.harm[["family"]]), NA,
                         ifelse(is.na(data.harm[["order"]]), "family", "order"))

        # create a new column" (db_taxon_higher)
        data.harm[["db_taxon_higher"]] <- "A"

        # add either NA, family name or order name based on rows

        # add NA
        data.harm[which(is.na(x.rows)), ][["db_taxon_higher"]] <- NA

        # add family name
        x.rows.f <- which(x.rows == "family")
        if (length(x.rows.f) > 0) {
            data.harm[x.rows.f, ][["db_taxon_higher"]] <- data.harm[x.rows.f, ][["family"]]
        }

        # add order name
        x.rows.o <- which(x.rows == "order")
        if (length(x.rows.o) > 0) {
            data.harm[x.rows.o, ][["db_taxon_higher"]] <- data.harm[x.rows.o, ][["order"]]
        }

        # select the relevant columns
        data.harm <- data.harm[, c(
            "row_id",
            "scientificName",
            "acceptedNameUsageID",
            "db_taxon_higher_rank",
            "db_taxon_higher"
        )]

        # remove the names that we were not able to resolve
        data.harm <- dplyr::filter(
            data.harm,
            !(is.na(scientificName) | is.na(db_taxon_higher_rank) | is.na(db_taxon_higher))
        )

        # join these data to the tax.dat data
        data <- dplyr::left_join(data, data.harm, by = c("row_id"))

        # add the special names back
        data <- dplyr::bind_rows(data, data.spec)
    } else {
        # if we only have special names, then we only consider the special names and
        # add additional columns for consistency
        data <- data.spec
        data[["scientificName"]] <- NA
        data[["acceptedNameUsageID"]] <- NA
        data[["db_taxon_higher_rank"]] <- NA
        data[["db_taxon_higher"]] <- NA
    }

    # remove the row_id column
    data <- dplyr::select(data, -row_id)

    # convert to a tibble
    data <- dplyr::as_tibble(data)

    data
}
