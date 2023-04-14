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
                                   max_tax_dist = 3,
                                   trait = "equation",
                                   gen_sp_dist = 0.5) {
  # make sure the max_tax_dist argument is a number > 0
  assert_that(is.number(max_tax_dist) & (max_tax_dist >= 0))
  assert_that(is.number(gen_sp_dist) & (gen_sp_dist >= 0))
  
  # make sure the trait chosen is supported
  assert_that(
    trait %in% c("equation", paste0("trait", 1:10)),
    msg = paste(trait, "is not a valid trait or equation, see documentation") # TODO: not yet in docs
  )
  
  # load the trait data
  if (!exists(paste0(trait, "_db"))) {
    assign(
      paste0(trait, "_db"),
      readRDS(file = get_db_file_path(paste0(trait, "_database.rds")))
    )
  }
  
  # assign the object to trait_db
  trait_db <- get(paste0(trait, "_db"))
  
  # split the input data.frame into a list
  data_list <- split(data, 1:nrow(data))
  
  # for each entry in the input.list, select appropriate traits
  output <- lapply(data_list, function(input) {
    
    # get the higher-level taxon databases if the input database is
    # one of the taxonomic backbones
    if (input[["db"]] %in% c("gbif", "itis", "col")) {
      db_name <- paste0(input[["db"]], "_db")
      td_name <- paste0(input[["db"]], "_td")
      
      if (!exists(db_name)) {
        assign(db_name, readRDS(file = get_db_file_path(
          paste0(input[["db"]], "_higher_taxon_matrices.rds")
        )))
        
      }
      
      htm_db <- get(db_name)
      
      if (!exists(td_name)) {
        assign(td_name, readRDS(file = get_db_file_path(
          paste0(input[["db"]], "_taxon_database.rds")
        )))
      }
      
      td_db <- get(td_name)
      
    } 
    
    # extract the target_name
    target_name <- input[["scientificName"]]
    
    target_name_cond <- extract_genus(target_name)
    
    # test if the target name is present in any of the taxon matrices
    target_present <- 
      
      sapply(htm_db, function(htm) {  
        
        target_name_cond %in% names(V(htm))
        
      })
    
    # target_present[1] <- TRUE
    if ( any( target_present == TRUE ) ) {
      
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
      dist_df <-
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
          
          dist_df <-
            dplyr::tibble(
              db_scientificName = db_taxon,
              trait_out = trait,
              id = id,
              tax_distance = tax_dist
            )
          
          return(dist_df)
          
        }, td[["scientificName"]], td[["id"]], SIMPLIFY = FALSE)
      
      # bind into a data.frame
      dist_df <- dplyr::bind_rows(dist_df)
      
      # filter out rows where the taxonomic distance is too large
      dist_df <- dplyr::filter(dist_df, tax_distance <= max_tax_dist)
      
    }
    
    # check if there are special names
    if (is.na(input[["scientificName"]]) & (input[["db"]] == "special")) {
      
      # get row_id's from trait database matching the special names
      row_id <- which(trait_db[["db_taxon"]] == input[[paste0("clean_", target_taxon)]])
      
      # check if there are rows that are outputted and if not return NA
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
      dist_df <- dplyr::tibble(
        db_scientificName = x,
        trait_out = trait,
        id = y,
        tax_distance = NA
      )
    } 
    
    if( all( target_present == FALSE ) ) {
      
      dist_df <- dplyr::tibble(
        db_scientificName = NA,
        trait_out = trait,
        id = NA,
        tax_distance = NA
      )
    }
    
    # add metadata
    dist_df <- dplyr::bind_cols(input, dist_df)
    
    dist_df
    
  })
  
  # test if all entries in the original data match the output list
  assert_that(
    all(unique(data[[target_taxon]]) == unique(unlist(sapply(output, function(x) {x[[target_taxon]]}))) ),
    msg = "number of unique taxa in input and output do not match"
  )
  
  output
  
}
