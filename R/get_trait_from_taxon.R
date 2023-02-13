#' @title get_trait_from_taxon()
#' @description Use taxonomic distance, life-stages, habitat data etc. to
#'  select an appropriate trait/equation
#' @details This function uses the [clean_taxon_names()] function to clean and
#'  harmonise the names of a target dataset. The Get_Habitat_Data() function is
#'  then used to get relevant habitat information. The
#'  [select_traits_tax_dist()] function calculates taxonomic distance for each
#'  focal taxon and the relevant entries in the different higher-level taxon
#'  graphs generated from the GBIF, ITIS and COL taxonomic backbones. This
#'  information is combined to choose an appropriate trait or equation for each
#'  target taxon name.
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' @param data - data.frame with at least five columns: target taxon, life
#'  stage, latitude (dd), longitude (dd) and body size (mm)
#'  if trait == "equation"
#' @param target_taxon - character string with the column name containing the
#'  taxon names
#' @param life_stage - character string with the column name containing the
#'  life-stage information
#' @param latitude_dd - character string with the column name containing the
#'  latitude in decimal degrees
#' @param longitude_dd - character string with the column name containing the
#'  longitude in decimal degrees
#' @param body_size - character string with the column name containing the
#'  body size data if trait = "equation"
#' @param workflow - options are "workflow1" or "workflow2"
#'  (default = "workflow2)
#' @param max_tax_dist - maximum taxonomic distance acceptable between the
#'  target and the taxa in the database (default = 3)
#' @param trait - trait to be searched for (default = "equation")
#' @param gen_sp_dist - taxonomic distance between a genus and a species
#'  (default = 0.5)
#' @return tibble with chosen traits or equations based on the input parameters
#' @export
#' @import assertthat
get_trait_from_taxon <- function(data,
                                 target_taxon,
                                 life_stage,
                                 latitude_dd,
                                 longitude_dd,
                                 body_size,
                                 workflow = "workflow2",
                                 max_tax_dist = 3,
                                 trait = "equation",
                                 gen_sp_dist = 0.5) {
  # make sure the max_tax_dist argument is a number
  assert_that(
    is.character(workflow) & (workflow %in% c("workflow1", "workflow2")),
    msg = paste(
      workflow,
      "must be a character string corresponding to: workflow1 or workflow2"
    )
  )

  # if the trait is an equation then we calculate the maximum and minimum
  # body size
  # if not, then we simply get the unique values for the different parameters

  if (trait == "equation") {
    assert_that(
      !is.na(body_size) & is.character(body_size),
      msg = "if trait = equation then body size data must be provided"
    )

    # create a data.unique data.frame
    data.unique <- data

    # get a unique identifier
    data.unique[["group_id"]] <- apply(
      data.unique[, c(target_taxon, life_stage, latitude_dd, longitude_dd)],
      1,
      function(x) paste0(x, collapse = "_")
    )

    # get min and max of body size for each unique combination of name,
    # life-stage and latitude-longitude

    # get a formula to describe the grouping used for the calculation
    form.agg <- stats::reformulate(termlabels = c("group_id"), body_size)

    # calculate the minimum body size per group
    x.min <- stats::aggregate(form.agg,
      data = data.unique, FUN = function(x) {
        if (all(is.na(x))) {
          return(NA)
        } else {
          return(min(x, na.rm = TRUE))
        }
      },
      na.action = stats::na.pass
    )

    # rename the body_size variable to min_body_size_mm
    names(x.min)[names(x.min) == body_size] <- "min_body_size_mm"

    # calculate the maximum body size per group
    x.max <- stats::aggregate(form.agg,
      data = data.unique, FUN = function(x) {
        if (all(is.na(x))) {
          return(NA)
        } else {
          return(max(x, na.rm = TRUE))
        }
      },
      na.action = stats::na.pass
    )

    # rename the body_size variable to max_body_size_mm
    names(x.max)[names(x.max) == body_size] <- "max_body_size_mm"

    # bind x.min and x.max
    x.min.max <- dplyr::full_join(x.min, x.max, by = "group_id")
    data.unique <- dplyr::full_join(data.unique, x.min.max, by = "group_id")

    # join the max and minimum body_size data together
    data.unique <- dplyr::distinct(data.unique[, c(
      target_taxon,
      life_stage,
      latitude_dd,
      longitude_dd,
      "min_body_size_mm",
      "max_body_size_mm"
    )])
  } else {
    data.unique <- dplyr::distinct(data[, c(
      target_taxon,
      life_stage,
      latitude_dd,
      longitude_dd
    )])
  }

  # add a row_id variable
  data.unique <- dplyr::bind_cols(
    dplyr::tibble(taxon_id = 1:nrow(data.unique)),
    data.unique
  )

  # run the get_habitat_data() function
  x1 <- get_habitat_data(data = data.unique, latitude_dd = latitude_dd, longitude_dd = longitude_dd)

  # set-up a vector of taxonomic databases
  db_vec <- c("gbif", "itis", "col")

  # clean the taxon names for each of the three taxonomic databases
  y1 <-
    lapply(db_vec, function(database) {
      cl.tax <- clean_taxon_names(
        data = x1,
        target_taxon = target_taxon, life_stage = life_stage,
        database = database
      )

      return(cl.tax)
    })

  # bind these rows into a single data.frame
  y1 <- dplyr::bind_rows(y1)

  # arrange by taxon_name
  y1 <- dplyr::arrange(y1, taxon_id)

  # remove any duplicates that can arise from the special name procedure
  y1 <- dplyr::distinct(y1)

  # run the Select_Traits_Tax_Dist() function
  z1 <- select_traits_tax_dist(data = y1, target_taxon = target_taxon)

  # bind the rows together
  z1 <- dplyr::bind_rows(z1)

  # get equation match data

  # load the trait data
  if (!exists(paste0(trait, "_db"))) {
    assign(
      paste0(trait, "_db"),
      readRDS(file = paste0(here::here("database"), "/", trait, "_database.rds"))
    )
  }

  # assign the object to trait_db
  trait_db <- get(paste0(trait, "_db"))

  # life_stage match
  life_stage_match <-
    mapply(function(x, y) {
      if (!is.na(x)) {
        return((trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["db_life_stage"]] == y))
      } else {
        return(NA)
      }
    }, z1[["id"]], z1[[life_stage]])

  # add life-stage match column
  z1[["life_stage_match"]] <- life_stage_match

  # additional matches that are only relevant for the equation trait
  if (trait == "equation") {
    # life_stage match
    r2_match <-
      sapply(z1[["id"]], function(x) {
        if (!is.na(x)) {
          return(trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["r2"]])
        } else {
          return(NA)
        }
      })

    z1[["r2_match"]] <- r2_match

    # body size range match
    body_size_range_match <-
      mapply(function(x, y, z) {
        if (!is.na(x)) {
          trait_db_sel <- trait_db[trait_db[[paste0(trait, "_id")]] == x, ]

          t1 <- (trait_db_sel[["body_size_min"]] <= y) & (trait_db_sel[["body_size_max"]] >= z)

          return(t1)
        } else {
          return(NA)
        }
      }, z1[["id"]], z1[["min_body_size_mm"]], z1[["max_body_size_mm"]])

    z1[["body_size_range_match"]] <- body_size_range_match
  }

  # get habitat match data

  # load the habitat database
  if (!exists("hab_db")) {
    hab_db <- readRDS(file = here::here("database/freshwater_ecoregion_data.rds"))
  }

  # select the correct trait from the habitat database
  hab_db_sel <- hab_db[hab_db[["database"]] == trait, ]

  # realm match
  realm_match <-
    mapply(function(x, y) {
      if (!is.na(x)) {
        return((hab_db_sel[hab_db_sel[["id"]] == x, ][["realm"]] == y))
      } else {
        return(NA)
      }
    }, z1[["id"]], z1[["realm"]])

  z1[["realm_match"]] <- realm_match

  # major habitat type match
  mht_match <-
    mapply(function(x, y) {
      if (!is.na(x)) {
        # subset the hab_db_sel data to only include the correct id
        hab_db_sel.sub <- hab_db_sel[hab_db_sel[["id"]] == x, ]
      } else {
        return(NA)
      }

      # test whether the equation has an ID and if the scale is correct
      if (!is.na(x) & (hab_db_sel.sub[["accuracy"]] %in% c("approximate", "exact"))) {
        return((hab_db_sel.sub[["major_habitat_type"]] == y))
      } else {
        return(NA)
      }
    }, z1[["id"]], z1[["major_habitat_type"]])

  z1[["major_habitat_type_match"]] <- mht_match

  # ecoregion match
  ecoregion_match <-
    mapply(function(x, y) {
      if (!is.na(x)) {
        # subset the hab_db_sel data to only include the correct id
        hab_db_sel.sub <- hab_db_sel[hab_db_sel[["id"]] == x, ]
      } else {
        return(NA)
      }

      # test whether the equation has an ID and if the scale is correct
      if (!is.na(x) & (hab_db_sel.sub[["accuracy"]] %in% c("exact"))) {
        return((hab_db_sel.sub[["ecoregion"]] == y))
      } else {
        return(NA)
      }
    }, z1[["id"]], z1[["ecoregion"]])

  z1[["ecoregion_match"]] <- ecoregion_match

  # split into a list
  z1.list <- split(z1, z1[["taxon_id"]])

  z1.select <-
    lapply(z1.list, function(input) {
      # if none of the id values are present, then return any row or else
      # remove the NAs
      if (all(is.na(input[["id"]]))) {
        input <- input[sample(1:nrow(input), 1), ]

        return(input)
      } else {
        input <- input[!is.na(input[["id"]]), ]
      }

      # if workflow one is chosen then return this input value
      if (workflow == "workflow1") {
        return(input)
      }

      # get the equations with matching life-stages
      if (sum(input[["life_stage_match"]] == TRUE, na.rm = TRUE) == 0) {
        input <- input[sample(1:nrow(input), 1), ]
        input[1, c("id", "tax_distance", names(input)[grepl(pattern = "_match", x = names(input))])] <- NA

        return(input)
      } else {
        input <- input[input[["life_stage_match"]] == TRUE & !is.na(input[["life_stage_match"]]), ]
      }

      # get the minimum taxonomic distance as long as the
      # difference is greater than 0.5
      if (all(!is.na(input[["tax_distance"]])) & (sum(!is.na(input[["tax_distance"]])) >= 1)) {
        input <- input[input[["tax_distance"]] <= (min(input[["tax_distance"]], na.rm = FALSE) + 0.5), ]
      }

      # if there is still no clear decision, then use the additional matches
      if (nrow(input) > 1) {
        # get the correct columns to match for the chosen trait
        if (trait == "equation") {
          match_cols <- c(
            "r2_match",
            "body_size_range_match",
            "realm_match",
            "major_habitat_type_match",
            "ecoregion_match"
          )
          weights <- c(1, 2, 0.75, 0.75, 0.75)
        } else {
          match_cols <- input[, c(
            "realm_match",
            "major_habitat_type_match",
            "ecoregion_match"
          )]
          weights <- c(0.75, 0.75, 0.75)
        }

        # calculate the match score based on the different match categories
        match_score <-
          apply(input[, match_cols], 1, function(x) {
            sum((x * weights), na.rm = TRUE)
          })

        input <- input[which(match_score == max(match_score)), ]
      }

      # if there are still multiple equations then we pick them randomly
      input <- input[sample(1:nrow(input), 1), ]

      return(input)
    })

  # bind these choices into a data.frame
  z1.select <- dplyr::bind_rows(z1.select)

  # add the traits or equation to the selected id numbers
  z1.select[[trait]] <-
    sapply(z1.select[["id"]], function(x) {
      if (!is.na(x)) {
        trait_db[trait_db[[paste0(trait, "_id")]] == x, ][[trait]]
      } else {
        NA
      }
    })

  # combine these chosen traits or equations with the original data
  z1.select <- dplyr::full_join(data, z1.select, by = c(
    target_taxon,
    life_stage,
    latitude_dd,
    longitude_dd
  ))

  # remove the taxon_id column
  z1.select <- dplyr::select(z1.select, -taxon_id)

  # parse the equation to calculate dry_biomass_mg
  if (trait == "equation") {
    z1.select[["dry_biomass_mg"]] <-
      mapply(function(x, y) {
        if (any(is.na(c(x, y)))) {
          return(NA)
        } else {
          var1 <- x
          eval(parse(text = y))
        }
      }, z1.select[[body_size]], z1.select[[trait]])
  }

  assert_that(
    all(unique(data[[target_taxon]]) == unique(z1.select[[target_taxon]])),
    msg = "number of unique taxa in input and output do not match"
  )

  z1.select
}
