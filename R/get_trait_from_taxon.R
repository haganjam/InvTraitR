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
  
  # add a row variable
  data <- dplyr::bind_cols(
    dplyr::tibble(row = 1:nrow(data)),
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
  clean_taxa <- dplyr::arrange(clean_taxa, row)
  
  # remove any duplicates that can arise from the special name procedure
  clean_taxa <- dplyr::distinct(clean_taxa)
  
  # we use the function to generate a decision df which we always output
  # this allows the user to check how we made the decisions
  decision_df <- select_traits_tax_dist(data = clean_taxa, 
                                        target_taxon = "scientificName",
                                        life_stage = life_stage,
                                        body_size = body_size,
                                        max_tax_dist = max_tax_dist,
                                        trait = trait,
                                        gen_sp_dist = gen_sp_dist
  )
  
  # subset the data according to the workflow
  if(workflow == "workflow1") {
    
    trait_sel <- dplyr::filter(decision_df, tax_distance <= max_tax_dist)
    
  } else if(workflow == "workflow2") {
    
    trait_sel <-
      decision_df |>
      dplyr::filter(workflow2_choice == TRUE) |>
      dplyr::group_by(row) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup()
    
  } else {
    
    stop("choose appropriate workflow")
    
  }
  
  # remove redundant columns
  trait_sel <- dplyr::select(trait_sel, -recommend, -explanation, -workflow2_choice)
  
  # load the trait database
  if (!exists(paste0(trait, "_db"))) {
    assign(
      paste0(trait, "_db"),
      readRDS(file = get_db_file_path(paste0(trait, "_database.rds")))
    )
  }
  
  # assign the object to trait_db
  trait_db <- get(paste0(trait, "_db"))
  
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
    trait_sel <- dplyr::left_join(trait_sel, cor_factors, by = "id")
    
  } else {
    
    # add the traits or equation to the selected id numbers
    trait_sel[[trait]] <-
      sapply(trait_sel[["id"]], function(x) {
        if (!is.na(x)) {
          trait_db[trait_db[[paste0(trait, "_id")]] == x, ][[trait]]
        } else {
          NA
        }
      })
  }
  
  if(nrow(trait_sel) > 0) {
    
    # parse the equation to calculate dry_biomass_mg
    if ( (trait == "equation") && (workflow == "workflow2") ) {
      
      # set-up a vector to capture the dry biomass values
      dry_biomass_mg <- vector(length = nrow(trait_sel))
      
      # loop over all the rows
      for(i in 1:nrow(trait_sel)) {
        
        # get the ith row of data
        L <- unlist(trait_sel[i, body_size], use.names = FALSE)
        model <- trait_sel[i,]$equation_form
        log_base <- trait_sel[i,]$log_base
        a <- trait_sel[i,]$a
        b <- trait_sel[i,]$b
        CF <- trait_sel[i,]$lm_correction
        scale <- trait_sel[i,]$dry_biomass_scale
        
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
      
      # add this dry biomass estimate to the data.frame
      trait_sel[["dry_biomass_mg"]] <- dry_biomass_mg
      
    }
    
  }
  
  # bind to the original data
  trait_sel <- dplyr::full_join(data, trait_sel, by = names(data))
  
  # make an output list
  output_list <- list(data = trait_sel, decision_data = decision_df)
  
  output_list
  
}
