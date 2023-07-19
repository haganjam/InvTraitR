#' @title extract_genus()
#' @description If a binomial taxa is supplied, extract the genus name only
#' @details This method works exclusively with genera and not species.
#'  This is because once the genus is known, the relationship among the
#'  species within that genus is known. It makes the method much more
#'  computationally efficient. This function is used to then extract the genus
#'  from a species name. The genus name is used to place the species in the
#'  correct taxonomic framework
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' @param binomial - binomial character string separated by a space (e.g.
#'  "Loxodonta loxodonta")
#' @return string with the genus name
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.string
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

#' @title extract_body_size_range_match()
#' @description calculate whether target length matches equation range
#' @details This function generates a data.frame with three variables:
#' body_size_range_match: whether the target length is within the equation 
#' length range (TRUE or FALSE), body_size_min_dist: percentage away 
#' from the minimum length used to create the equation, body_size_max_dist: 
#' percentage away from the maximum length.
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' @param equation_id - id number of the equation
#' @param target_body_size - body size of the target taxon
#' @param prop - proportion above and below equation body size range
#' considered acceptable for using an equation
#' @param equation_db - equation database object
#' @return string with the genus name
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.string
extract_body_size_range_match <- function(equation_id, 
                                          target_body_size,
                                          prop = 0.30,
                                          equation_db) {
  
  if (!is.na(equation_id)) {
    
    # extract relevant equation from database
    equ_meta <- equation_db[equation_db[["equation_id"]] == equation_id, ]
    
    # calculate the minimum and maximum acceptable body size ranges
    equ_min <- equ_meta[["body_size_min"]]
    equ_min <- equ_min - (equ_min*prop)
    
    equ_max <- equ_meta[["body_size_max"]]
    equ_max <- equ_max + (equ_max*prop)
    
    # check whether target_body_length is within the equation range
    body_size_range_match <- (target_body_size >= equ_min) & (target_body_size <= equ_max)
    
    body_size_range_match
    
  } else {
    
    NA
    
  }
  
}

#' @title select_traits_tax_dist()
#' @description Get taxonomic distances of target names relative to the taxa
#'  databases
#' @details This function searches the relevant trait or equation database for
#'  the best matching trait or equation for a given target name based on 
#'  taxonomic distance
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' @param data - input data.frame exported from [get_habitat_data()] and
#'  [clean_taxon_names()] function
#' @param target_taxon - character string with the column name containing the
#'  taxon names
#' @param life_stage - character string with the column name containing the
#'  life-stage information
#' @param body_size - column name containing the body length data
#' @param max_tax_dist - maximum taxonomic distance acceptable between the
#'  target and the taxa in the database (default = 3)
#' @param trait - trait to be searched for (default = "equation")
#' @param gen_sp_dist - taxonomic distance between a genus and a species
#'  (default = 0.5)
#' @return tibble of the input data with traits or equations within the maximum
#'  taxonomic distance
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.number
select_traits_tax_dist <- function(data,
                                   target_taxon,
                                   body_size,
                                   life_stage,
                                   max_tax_dist = 3,
                                   trait = "equation",
                                   gen_sp_dist = 0.5
                                   ) {
  # make sure the max_tax_dist argument is a number > 0
  assert_that(is.number(max_tax_dist) & (max_tax_dist >= 0))
  assert_that(is.number(gen_sp_dist) & (gen_sp_dist >= 0))
  
  # make sure the trait chosen is supported
  assert_that(
    trait %in% c("equation", paste0("trait", 1:10)),
    msg = paste(trait, "is not a valid trait or equation, see documentation") # TODO: not yet in docs
  )
  
  # add a column to the dataset specifying the trait being searched
  data[["trait_out"]] <- trait
  
  # add columns to fill
  new_cols <- c("db_scientificName", "id", "tax_distance", "body_size_range_match",
                "life_stage_match", "r2_match", "n", "db_min_body_size_mm", "db_max_body_size_mm",
                "realm_match", "major_habitat_type_match", "ecoregion_match",
                "recommend", "explanation", "workflow2_choice")  
  for(i in 1:length(new_cols)) {
    data[[new_cols[i]]] <- NA
  }
  
  # load the trait data
  if (!exists(paste0(trait, "_db"))) {
    assign(
      paste0(trait, "_db"),
      readRDS(file = get_db_file_path(paste0(trait, "_database.rds")))
    )
  }
  
  # assign the object to trait_db
  trait_db <- get(paste0(trait, "_db"))
  
  # load the higher taxon matrices and taxon databases
  db_vec <- c("gbif", "itis", "col")
  for(i in 1:length(db_vec)) {
    
    db_name <- paste0(db_vec[i], "_db")
    if (!exists(db_name)) {
      assign(db_name, readRDS(file = get_db_file_path(
        paste0(db_vec[i], "_higher_taxon_matrices.rds")
      )))}
    
    td_name <- paste0(db_vec[i], "_td")
    if (!exists(td_name)) {
      assign(td_name, readRDS(file = get_db_file_path(
        paste0(db_vec[i], "_taxon_database.rds")
      )))
    }
    
  }
  
  # get the rows with special names
  data_spec <- dplyr::filter(data, db == "special")
  
  # if there are special names
  if(nrow(data_spec) > 0) {
    
    # split into a list
    data_spec_list <- split(data_spec, 1:nrow(data_spec))
    
    # check that the data split worked
    assert_that(
      length(data_spec_list) == nrow(data_spec),
      msg = paste(trait, "list conversion did not work")
    )
    
    # loop over the special names
    output_spec <- 
      
      lapply(data_spec_list, function(input) {
        
        # get row_id's from trait database matching the special names
        row_id <- which(trait_db[["db_taxon"]] == input[[paste0("clean_", target_taxon)]])
        
        # if there are rows present, then:
        # add the database taxon name and the trait or equation
        if (length(row_id) != 0) {
          input <- 
            lapply(row_id, function(x) {
              
              input[["db_scientificName"]] <- trait_db[x, ][["db_taxon"]]
              input[["id"]] <- trait_db[x, ][[paste0(trait, "_id")]]
              
              input
              
            })
          input <- dplyr::bind_rows(input)
        } else {
          input[["explanation"]] <- "no appropriate special names in database"
        }
        
        input
        
      })
    
  } else {
    output_spec <- NULL
  }
  
  # get a data.frame without special names
  data <- dplyr::filter(data, db != "special")
  
  # check that the data split worked
  assert_that(
    all(data[["db"]] %in% c("gbif", "itis", "col")),
    msg = paste(trait, "data contains entries in the db column that are not supported")
  )
  
  # if there are non-special names
  if(nrow(data) > 0) {
    
    # split the input data.frame into a list
    data_list <- split(data, 1:nrow(data))
    
    # check that the data split worked
    assert_that(
      length(data_list) == nrow(data),
      msg = paste(trait, "list conversion did not work")
    )
    
    output <- lapply(data_list, function(input) {
      
      # get the relevant taxonomic backbones
      htm_db <- get(paste0(input[["db"]], "_db"))
      td_db <- get(paste0(input[["db"]], "_td") )
      
      # extract the target_name
      target_name <- input[["scientificName"]]
      target_name_cond <- extract_genus(target_name)
      
      # test if the target name is present in any of the taxon matrices
      target_present <- sapply(htm_db, function(htm) {
        target_name_cond %in% names(igraph::V(htm))
      })
      
      # if the target taxon is not present in any backbone then return an explanation
      if(is.na(target_name_cond)) {
        
        input[["explanation"]] <- "target name not present in taxonomic backbone"
        
      } else if(all(target_present == FALSE)) {
        
        input[["explanation"]] <- "target name not found in taxonomic graph"
        
      } else if (any(target_present == TRUE)) {
        
        # get the relevant taxon matrix
        higher_taxon <- names(htm_db[target_present])
        htm <- htm_db[target_present][[1]]
        
        # extract vertices
        v_x <- igraph::V(htm)
        
        # extract the equation entries from the taxon database
        td <- td_db[td_db$database == trait, ]
        
        # extract the entries from the equation taxon database matching the target higher taxon
        td <- td[td$order == higher_taxon | td$family == higher_taxon , ]
        
        # remove the NA values
        td <- td[ !(is.na(td$order) & is.na(td$family)), ]
        
        # taxonomic distance
        td_vec <-
          mapply(function(db_taxon, id) {
            if (db_taxon == target_name) {
              tax_dist <- 0
            } else {
              # extract genus for species-level names
              db_taxon_cond <- extract_genus(db_taxon)
              
              tax_dist <-
                igraph::distances(htm,
                                  v_x[which(attr(v_x, "names") == target_name_cond)],
                                  v_x[which(attr(v_x, "names") == db_taxon_cond)],
                                  mode = c("all"),
                                  algorithm = c("bellman-ford")
                )
              
              # if length is zero then the distance is NA
              if (length(tax_dist) == 0) {
                tax_dist <- NA
              } else {
                tax_dist <- tax_dist[[1]]
              }
              
              # extra distance for species level: gen_sp_dist argument
              sp_l <- sum(ifelse(c(attr(target_name_cond, "n"), attr(db_taxon_cond, "n")) > 1, gen_sp_dist, 0))
              
              # add extra distance
              tax_dist <- tax_dist + sp_l
              
            }
            
          }, td[["scientificName"]], td[["id"]], SIMPLIFY = FALSE)
        
        # add the equation taxon name, id and taxonomic distance to the input data
        input_list <- vector("list", length = nrow(td))
        for(i in 1:nrow(td)) {
          # copy the input data row
          input_mod <- input
          
          # add relevant information
          input_mod[["db_scientificName"]] <- td[["scientificName"]][i]
          input_mod[["id"]] <- td[["id"]][i]
          input_mod[["tax_distance"]] <- td_vec[[i]]
          
          # add the input row to the input_list object
          input_list[[i]] <- input_mod
          
        }
        
        # bind the list into a data.frame
        input <- dplyr::bind_rows(input_list)
        
      }
      input
      
    })
    
  } else {
    output <- NULL
  }
  
  # bind the list of regular names and the list of special names
  if( (length(output_spec) > 0) && (length(output) > 0) ) {
    output <- c(output, output_spec)
  } else if( (length(output_spec) > 0) && (length(output) == 0) ) {
    output <- output_spec
  }
  
  # bind the output list into a data.frame
  output <- dplyr::bind_rows(output)
  
  # get the life-stage matches
  output[["life_stage_match"]] <- 
    mapply(function(x, y) {
      if (!is.na(x)) {
        trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["db_life_stage"]] == y
      } else {
        NA
      }
    }, output[["id"]], output[[life_stage]])
  
  # get the body size range matches if the trait is an equation
  output[["body_size_range_match"]] <- 
    mapply(function(x, y) {
      # extract the body size range match information
      extract_body_size_range_match(equation_id = x,
                                    target_body_size = y,
                                    equation_db = trait_db)
    }, output[["id"]], output[[body_size]])
  
  # additional matches
  
  # set-up a vector of the relevant columns and relevant names
  rel_cols <- c("r2", "n", "body_size_min", "body_size_max")
  rel_names <- c("r2_match", "n", "db_min_body_size_mm", "db_max_body_size_mm")
  
  # loop over these variables
  for (i in 1:length(rel_cols)) {
    
    output[[rel_names[i]]] <- 
      sapply(output[["id"]], function(x) {
        if (!is.na(x)) {
          return(trait_db[trait_db[[paste0(trait, "_id")]] == x, ][[rel_cols[i]]])
        } else {
          return(NA)
        }
      })
  }
  
  # get habitat match data
  
  # load the habitat database
  if (!exists("hab_db")) {
    hab_db <- readRDS(file = get_db_file_path("freshwater_ecoregion_data.rds"))
  }
  
  hab_cols <- c("realm", "major_habitat_type", "ecoregion")
  hab_names <- paste0(hab_cols, "_match")
  
  # select the correct trait from the habitat database
  hab_db_sel <- hab_db[hab_db[["database"]] == trait, ]
  
  # loop over these variables
  for (i in 1:length(hab_cols)) {
    
    output[[hab_names[i]]] <- 
      mapply(function(x, y) {
        
        # if there is a valid equation present
        if (!is.na(x)) {
          hab_sub <- hab_db_sel[hab_db_sel[["id"]] == x, ]
          
          # get the habitat match dat
          hab_match <- hab_sub[[hab_cols[i]]] == y
          
          # if the accuracy of the coordinates are approximate, then do not return
          # an ecoregion variable
          if( (hab_cols[i] == "ecoregion") & (hab_sub[["accuracy"]] == "approximate") ) {
            hab_match <- NA
          }
          
        } else {
          hab_match <- NA
        }
        
        return(hab_match)
        
      }, output[["id"]], output[[hab_cols[i]]])
    
  }
  
  # set-up explanations
  explan_vec <- c("equation not within max taxonomic distance", 
                  "body size is not within body size range", 
                  "life stage does not match equation")
  
  # add an explanation for why different equations were not chosen
  for(i in 1:nrow(output)) {
    
    if(!is.na(output[i,][["id"]])) {
      
      x <- (output[i,][["tax_distance"]] <= max_tax_dist)
      y <- (output[i,][["body_size_range_match"]])
      z <- (output[i,][["life_stage_match"]])
      
      if(output[i, ][["db"]] == "special") {
        
        output[i,][["explanation"]] <- if(any((!c(y, z)) == TRUE)) {
          paste(explan_vec[!c(y, z)], collapse = " & ")
        } else { 
          NA 
        } 
        
        output[i,][["recommend"]] <- all(c(y, z))
        
      } else{
        
        output[i,][["explanation"]] <- if(any((!c(x, y, z)) == TRUE)) {
          paste(explan_vec[!c(x, y, z)], collapse = " & ")
        } else { 
          NA 
        } 
        
        output[i,][["recommend"]] <- all(c(x, y, z))
        
      }
      
    } else {
      
      output[i,][["recommend"]] <- FALSE
      
    }
    
  }
  
  # convert output data.frame to a list
  output_list <- split(output, output[["row"]])
  
  # which names are special?
  w_spec <- sapply(output_list, function(x) all(x[["db"]] == "special"))
  
  # split the list into special and non-special names
  output_list_reg <- output_list[!w_spec]
  output_list_spec <- output_list[w_spec]
  
  # annotate each potential equation
  output_list_reg <- 
    lapply(output_list_reg, function(input) {
      
      # get the recommended equations with the lowest taxonomic distance
      min_td <- 
        if( all(is.na(input[["tax_distance"]])) ) {
          FALSE
        } else {
          with(input, (recommend == TRUE) & 
                 (tax_distance <= (min(tax_distance, na.rm = TRUE)+0.25) ) )
        }
      
      if(sum(min_td, na.rm = TRUE) > 1) {
        
        # get additional matches
        match_cols <- c(
          "r2_match",
          "realm_match",
          "major_habitat_type_match",
          "ecoregion_match"
        )
        weights <- c(1, (1/3), (1/3), (1/3))
        
        # calculate the match score based on the different match categories
        add_score <-
          apply(input[, match_cols], 1, function(x) {
            sum((x * weights), na.rm = TRUE)
          })
        
        # choose the best one equation
        input[["workflow2_choice"]] <- 
          dplyr::near(add_score, 
                      max(add_score[min_td], na.rm = TRUE)) & min_td
        
      } else {
        
        input[["workflow2_choice"]] <- min_td
        
      }
      
      input
      
    } )
  
  # annotate each potential equation
  output_list_spec <- 
    lapply(output_list_spec, function(input) {
      
      # get the recommended equations with the lowest taxonomic distance
      add_score <- 
        if( all(input[["recommend"]] == FALSE) ) {
          FALSE
        } else {
          with(input, body_size_range_match + life_stage_match)
        }
      
      # if both of those are true then it is a workflow choice
      input[["workflow2_choice"]] <- (add_score == 2)
      
      input
      
    })
  
  # bind the list into a data.frame
  output_df <- dplyr::bind_rows( c(output_list_reg, output_list_spec) )
  
  # return the data.frame
  output_df
  
}
