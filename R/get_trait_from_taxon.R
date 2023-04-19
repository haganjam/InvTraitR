#' @title get_trait_from_taxon()
#' @description Use taxonomic distance, life-stages, habitat data etc. to
#'  select an appropriate trait/equation
#' @details This function uses the [clean_taxon_names()] function to clean and
#'  harmonise the names of a target dataset. The [get_habitat_data()] function is
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
#' @importFrom assertthat assert_that
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
  
  # add a row_id variable
  data <- dplyr::bind_cols(
    dplyr::tibble(taxon_id = 1:nrow(data)),
    data
  )
  
  # run the get_habitat_data() function: x1
  hab_dat <- get_habitat_data(data = data, latitude_dd = latitude_dd, longitude_dd = longitude_dd)
  
  # set-up a vector of taxonomic databases
  db_vec <- c("gbif", "itis", "col")
  
  # clean the taxon names for each of the three taxonomic databases: y1
  clean_taxa <-
    lapply(db_vec, function(database) {
      cl_tax <- clean_taxon_names(
        data = hab_dat,
        target_taxon = target_taxon, life_stage = life_stage,
        database = database
      )
      
      return(cl_tax)
    })
  
  # bind these rows into a single data.frame
  clean_taxa <- dplyr::bind_rows(clean_taxa)
  
  # arrange by taxon_name
  clean_taxa <- dplyr::arrange(clean_taxa, taxon_id)
  
  # remove any duplicates that can arise from the special name procedure
  clean_taxa <- dplyr::distinct(clean_taxa)
  
  # run the select_traits_tax_dist() function: z1
  trait_sel <- select_traits_tax_dist(data = clean_taxa, 
                                      target_taxon = target_taxon,
                                      max_tax_dist = max_tax_dist,
                                      trait = trait,
                                      gen_sp_dist = gen_sp_dist
  )
  
  # bind the rows together
  trait_sel <- dplyr::bind_rows(trait_sel)
  
  # get equation match data
  
  # load the trait data
  if (!exists(paste0(trait, "_db"))) {
    assign(
      paste0(trait, "_db"),
      readRDS(file = get_db_file_path(paste0(trait, "_database.rds")))
    )
  }
  
  # assign the object to trait_db
  trait_db <- get(paste0(trait, "_db"))
  
  # life_stage match
  life_stage_match <- mapply(function(x, y) {
    if (!is.na(x)) {
      return((trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["db_life_stage"]] == y))
    } else {
      return(NA)
    }
  }, trait_sel[["id"]], trait_sel[[life_stage]])
  
  # add life-stage match column
  trait_sel[["life_stage_match"]] <- life_stage_match
  
  # additional matches that are only relevant for the equation trait
  if (trait == "equation") {
    # r2 value match
    r2_match <- sapply(trait_sel[["id"]], function(x) {
      if (!is.na(x)) {
        return(trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["r2"]])
      } else {
        return(NA)
      }
    })
    
    trait_sel[["r2_match"]] <- r2_match
    
    # sample size match
    N <- sapply(trait_sel[["id"]], function(x) {
      if (!is.na(x)) {
        return(trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["n"]])
      } else {
        return(NA)
      }
    })
    
    trait_sel[["N"]] <- N 
    
    }
  
  # add equation min body size
  min_bs <- sapply(trait_sel[["id"]], function(x) {
    if (!is.na(x)) {
      return(trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["body_size_min"]])
    } else {
      return(NA)
    }
  })
  
  trait_sel[["db_min_body_size_mm"]] <- min_bs
  
  # add equation max body size
  max_bs <- sapply(trait_sel[["id"]], function(x) {
    if (!is.na(x)) {
      return(trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["body_size_max"]])
    } else {
      return(NA)
    }
  })
  
  trait_sel[["db_max_body_size_mm"]] <- max_bs
  
  
  # get habitat match data
  
  # load the habitat database
  if (!exists("hab_db")) {
    hab_db <- readRDS(file = get_db_file_path("freshwater_ecoregion_data.rds"))
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
    }, trait_sel[["id"]], trait_sel[["realm"]])
  
  trait_sel[["realm_match"]] <- realm_match
  
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
    }, trait_sel[["id"]], trait_sel[["major_habitat_type"]])
  
  trait_sel[["major_habitat_type_match"]] <- mht_match
  
  # ecoregion match
  ecoregion_match <-
    mapply(function(x, y) {
      if (!is.na(x)) {
        # subset the hab_db_sel data to only include the correct id
        hab_db_sel_sub <- hab_db_sel[hab_db_sel[["id"]] == x, ]
      } else {
        return(NA)
      }
      
      # test whether the equation has an ID and if the scale is correct
      if (!is.na(x) & (hab_db_sel_sub[["accuracy"]] %in% c("exact"))) {
        return((hab_db_sel_sub[["ecoregion"]] == y))
      } else {
        return(NA)
      }
    }, trait_sel[["id"]], trait_sel[["ecoregion"]])
  
  trait_sel[["ecoregion_match"]] <- ecoregion_match
  
  # split into a list
  trait_sel_list <- split(trait_sel, trait_sel[["taxon_id"]])
  
  trait_sel_select <-
    lapply(trait_sel_list, function(input) {
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
      if (all(!is.na(input[["tax_distance"]])) && (sum(!is.na(input[["tax_distance"]])) >= 1)) {
        input <- input[input[["tax_distance"]] <= (min(input[["tax_distance"]], na.rm = FALSE) + 0.5), ]
      }
      
      # if there is still no clear decision, then use the additional matches
      if (nrow(input) > 1) {
        # get the correct columns to match for the chosen trait
        if (trait == "equation") {
          match_cols <- c(
            "r2_match",
            "realm_match",
            "major_habitat_type_match",
            "ecoregion_match"
          )
          weights <- c(1, (1/3), (1/3), (1/3))
        } else {
          match_cols <- input[, c(
            "realm_match",
            "major_habitat_type_match",
            "ecoregion_match"
          )]
          weights <- c((1/3), (1/3), (1/3))
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
  trait_sel_select <- dplyr::bind_rows(trait_sel_select)
  
  # if the trait is an equation then we add all the relevant equation information
  # if the trait is not an equation, then we simply add the trait value
  if(trait == "equation") {
    
    # extract the relevant columns
    cor_factors <- trait_db[, c("equation_id", 
                                "preservation",
                                "equation_form", "log_base", "a", "b",  
                                "lm_correction", "lm_correction_type",
                                "dry_biomass_scale")]
    names(cor_factors)[1] <- "id"
    
    # join these columns to the trait_sel_select
    trait_sel_select <- dplyr::left_join(trait_sel_select, cor_factors, by = "id")
    
  } else {
    
    # add the traits or equation to the selected id numbers
    trait_sel_select[[trait]] <-
      sapply(trait_sel_select[["id"]], function(x) {
        if (!is.na(x)) {
          trait_db[trait_db[[paste0(trait, "_id")]] == x, ][[trait]]
        } else {
          NA
        }
      })
  }
  
  # remove the taxon_id column
  trait_sel_select <- dplyr::select(trait_sel_select, -taxon_id)
  
  # parse the equation to calculate dry_biomass_mg
  if (trait == "equation") {
    
    # set-up a vector to capture the dry biomass values
    dry_biomass_mg <- vector(length = nrow(trait_sel_select))
    
    # loop over all the rows
    for(i in 1:nrow(trait_sel_select)) {
      
      # get the ith row of data
      L <- unlist(trait_sel_select[i, body_size], use.names = FALSE)
      model <- trait_sel_select[i,]$equation_form
      log_base <- trait_sel_select[i,]$log_base
      a <- trait_sel_select[i,]$a
      b <- trait_sel_select[i,]$b
      CF <- trait_sel_select[i,]$lm_correction
      scale <- trait_sel_select[i,]$dry_biomass_scale
      
      # evalulate the equation
      if( any(is.na(c(a, b))) )  {
        
        dry_biomass_mg[i] <- NA
        
      } else if (model == "model1") {
        
        # calculate the raw prediction on the log-scale
        x <- a + (b*logb(x = L, base = log_base))
        
        # convert to the natural scale
        x <- (log_base^x)
        
        # apply the correction factor
        dry_biomass_mg[i] <- ifelse(!is.na(CF), x*CF, x)*scale
        
      } else if (model == "model2") {
        
        # calculate the raw prediction
        dry_biomass_mg[i] <- a*(L^b)*scale
        
      }
      
    }
    
  }
  
  # add this dry biomass estimate to the data.frame
  trait_sel_select[["dry_biomass_mg"]] <- dry_biomass_mg
  
  assert_that(
    all(unique(data[[target_taxon]]) == unique(trait_sel_select[[target_taxon]])),
    msg = "number of unique taxa in input and output do not match"
  )
  
  trait_sel_select
  
}
