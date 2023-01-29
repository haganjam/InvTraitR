#' @title extract_genus()
#' @description If a binomial taxa is supplied, extract the genus name only
#' @details This method works almost exclusively with genera and not species.
#'  This is because once the genus is known, the relationship among the
#'  species within that genus is known. It makes the method much more
#'  computationally efficient. This function is used to then extract the genus
#'  from a species name. The genus name is used to place the species in the
#'  correct taxonomic framework
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' @param binomial - binomial character string separated by a space (e.g.
#'  "Loxodonta loxodonta")
#' @return string with the genus name
#' @import assertthat
extract_genus <- function(binomial) {
  # input validation: string w/ only letters
  assert_that(
    is.string(binomial),
    msg = paste(binomial, "is not a string")
  )
  assert_that(
    !grepl(pattern = "[[:punct:]]", binomial),
    msg = paste(binomial, "contains special characters")
  )

  # split the binomial into separate parts
  binomial.1st <- unlist(strsplit(x = binomial, split = " ", fixed = TRUE))

  # calculate the length of the split object
  binomial.l <- length(binomial.1st)

  # if the resulting object has a length of greater than 1, then extract
  # first element
  if (binomial.l > 1) {
    binomial <- binomial.1st[1]
  }

  # add a word count attribute
  attr(binomial, "n") <- binomial.l

  # return the modified name
  binomial
}

#' @title select_traits_tax_dist()
#' @description Get taxonomic distances of target names relative to the taxa
#'  databases
#' @details This function searches the relevant trait or equation database for
#'  the best
#' matching trait or equation for a given target name based on three criteria:
#'  taxonomic distance life-stage match and habitat match
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' @param data - input data.frame exported from [get_habitat_data()] and
#'  clean_taxon_names() function
#' @param target_taxon - character string with the column name containing the
#'  taxon names
#' @param max_tax_dist - maximum taxonomic distance acceptable between the
#'  target and the taxa in the database (default = 3)
#' @param trait - trait to be searched for (default = "equation")
#' @param gen_sp_dist - taxonomic distance between a genus and a species
#'  (default = 0.5)
#' @return tibble of the input data with traits or equations within the maximum
#'  taxonomic distance
#' @import igraph
#' @import assertthat
select_traits_tax_dist <- function(data,
                                   target_taxon,
                                   max_tax_dist = 3,
                                   trait = "equation",
                                   gen_sp_dist = 0.5) {
    # make sure the max_tax_dist argument is a number > 0
    assert_that(is.number(max_tax_dist) & (max_tax_dist >= 0))
    assert_that(is.number(gen_sp_dist) & (gen_sp_dist >= 0))

    # make sure the trait chosen is supported
    test_3 <- function(x) {
        trait %in% c("equation", paste0("trait", 1:10))
    }
    assertthat::on_failure(test_3) <- function(call, env) {
        paste0(
            deparse(call$x),
            " is not a valid trait or equation, see documentation"
        )
    }
    assertthat::assert_that(test_3(trait))

    # load the trait data
    if (!exists(paste0(trait, "_db"))) {
        assign(
            paste0(trait, "_db"),
            readRDS(file = paste0(here::here("database"), "/", trait, "_database.rds"))
        )
    }

    # assign the object to trait_db
    trait_db <- get(paste0(trait, "_db"))

    # split the input data.frame into a list
    data.list <- split(data, 1:nrow(data))

    # for each entry in the input.list, select appropriate traits
    output <-
        lapply(data.list, function(input) {
            # get the higher-level taxon databases if the input database is one of
            # the taxonomic backbones
            if (input[["db"]] %in% c("gbif", "itis", "col")) {
                db_name <- paste0(input[["db"]], "_db")
                td_name <- paste0(input[["db"]], "_td")

                if (!exists(db_name)) {
                    assign(db_name, readRDS(file = paste0(here::here("database"), "/", input[["db"]], "_higher_taxon_matrices.rds")))
                }

                htm_db <- get(db_name)

                if (!exists(td_name)) {
                    assign(td_name, readRDS(file = paste0(here::here("database"), "/", input[["db"]], "_taxon_database.rds")))
                }

                td_db <- get(td_name)
            } else {
                htm_db <- NA
            }

            if (!is.na(input[["scientificName"]]) & !is.na(input[["db_taxon_higher"]]) & (input[["db_taxon_higher"]] %in% names(htm_db))) {
                # taxon matrix
                htm <- htm_db[(names(htm_db) == input[["db_taxon_higher"]])][[1]]

                # extract vertices
                v.x <- V(htm)

                # extract the equation entries from the taxon database
                td <- td_db[td_db$database == trait, ]

                # extract the entries from the equation taxon database matching the target higher taxon
                td <- td[td$db_taxon_higher == input[["db_taxon_higher"]], ]

                # extract the target.name
                target.name <- input[["scientificName"]]

                # taxonomic distance
                dist.df <-
                    mapply(function(db.name, id) {
                        if (db.name == target.name) {
                            tax.dist <- 0
                        } else {
                            # extract genus for species-level names
                            target.name2 <- extract_genus(target.name)
                            db.name2 <- extract_genus(db.name)

                            tax.dist <-
                                igraph::distances(htm,
                                    v.x[which(attr(v.x, "names") == target.name2)],
                                    v.x[which(attr(v.x, "names") == db.name2)],
                                    mode = c("all"),
                                    algorithm = c("bellman-ford")
                                )

                            # if length is zero then the distance is zero
                            if (length(tax.dist) == 0) {
                                tax.dist <- NA
                            } else {
                                tax.dist <- tax.dist[[1]]
                            }

                            # extra distance for species level: gen_sp_dist argument
                            sp.l <- sum(ifelse(c(attr(target.name2, "n"), attr(db.name2, "n")) > 1, gen_sp_dist, 0))

                            # add extra distance
                            tax.dist <- tax.dist + sp.l
                        }

                        dist.df <-
                            dplyr::tibble(
                                db.scientificName = db.name,
                                trait_out = trait,
                                id = id,
                                tax_distance = tax.dist
                            )

                        return(dist.df)
                    }, td[["scientificName"]], td[["id"]], SIMPLIFY = FALSE)

                # bind into a data.frame
                dist.df <- dplyr::bind_rows(dist.df)

                # remove the rows where tax_distance is too large
                dist.df <- dplyr::filter(dist.df, tax_distance <= max_tax_dist)
            } else if (is.na(input[["scientificName"]]) & (input[["db"]] == "special")) {
                # get row_id's from trait database matching the special names
                row_id <- which(trait_db[["db_taxon"]] == input[[paste0("clean_", target_taxon)]])

                # check if there are rows that outputted and if not return NA
                x <- if (length(row_id) == 0) {
                    NA
                } else {
                    trait_db[row_id, ][["db_taxon"]]
                }
                y <- if (length(row_id) == 0) {
                    NA
                } else {
                    trait_db[row_id, ][[paste0(trait, "_id")]]
                }

                # pull this into a data.frame
                dist.df <- dplyr::tibble(
                    db.scientificName = x,
                    trait_out = trait,
                    id = y,
                    tax_distance = NA
                )
            } else {
                dist.df <- dplyr::tibble(
                    db.scientificName = NA,
                    trait_out = trait,
                    id = NA,
                    tax_distance = NA
                )
            }

            # add metadata
            dist.df <- dplyr::bind_cols(input, dist.df)

            return(dist.df)
        })

    output
}
