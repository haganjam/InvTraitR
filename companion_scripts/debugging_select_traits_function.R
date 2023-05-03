
data <- dat[sample(1:nrow(dat), 20),]
head(data)

target_taxon = "taxon"
life_stage = "life_stage"
latitude_dd = "lat"
longitude_dd = "lon"
body_size = "length_mm"
workflow = "workflow2"
max_tax_dist = 3
trait = "equation"
gen_sp_dist = 0.5

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


data <- clean_taxa

# add a special name row
x <- data[5,]
x$clean_taxon <- "Nematoda"
x$db <- "special"
x$scientificName <- NA
x$life_stage <- "none"
x$row <- 21

# add this to the data object
data <- dplyr::bind_rows(data, x)


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
      input[["db_scientificName"]] <- trait_db[row_id, ][["db_taxon"]]
      input[["id"]] <- trait_db[row_id, ][[paste0(trait, "_id")]]
    } else {
      input[["explanation"]] <- "no appropriate special names in database"
    }
    
    input
    
  })
  
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
  
  output <- 
    
    lapply(data_list, function(input) {
    
    # get the relevant taxonomic backbones
    htm_db <- get(paste0(input[["db"]], "_db"))
    td_db <- get(paste0(input[["db"]], "_td") )
    
    # extract the target_name
    target_name <- input[["scientificName"]]
    target_name_cond <- extract_genus(target_name)
    
    # test if the target name is present in any of the taxon matrices
    target_present <- 
      sapply(htm_db, function(htm) {  
        
        target_name_cond %in% names(V(htm))
        
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
    
  stop("data.frame must have some taxon names to search")
  
  }

# bind the list of regular names and the list of special names
if( exists(x = "output_spec") ) {
  output <- c(output, output_spec)
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
    
    output[i,][["explanation"]] <- if(any((!c(x, y, z)) == TRUE)) {
      paste(explan_vec[!c(x, y, z)], collapse = " & ")
    } else { 
      NA 
    } 
    
    output[i,][["recommend"]] <- all(c(x, y, z))
    
  } else {
    
    output[i,][["recommend"]] <- FALSE
    
  }
  
  
  }

# split by row id
output_list <- split(output, output[["row"]])

# deal with NAs in the id column
lapply(output_list, function(input) {
  
  # get the recommended equations with the lowest taxonomic distance
  min_td <- with(input, (recommend == TRUE) & 
                   (tax_distance <= (min(tax_distance, na.rm = TRUE)+0.25) ) )
  
  if(sum(min_td) > 1) {
    
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
    input[["workflow2_choice"]] <- (dplyr::near(add_score, 
                                                max(add_score[min_td], na.rm = TRUE)) & min_td)
    
  } else {
    
    input[["workflow2_choice"]] <- min_td
    
  }
  
  input
  
  } )

data

#






d# for each entry in the input.list, select appropriate traits
output <- lapply(data_list, function(input) {
  
    input <- data_list[[1]]
  
    htm_db <- get(paste0(input[["db"]], "_db"))
    td_db <- get(paste0(input[["db"]], "_td") )
    
    # extract the target_name
    target_name <- input[["scientificName"]]
    
    target_name_cond <- extract_genus(target_name)
    
    # test if the target name is present in any of the taxon matrices
    target_present <- 
      
      sapply(htm_db, function(htm) {  
        
        target_name_cond %in% names(V(htm))
        
      }) 
  
  # check 
  if (all(target_present == FALSE) & (input[["db"]] != "special")) {
    
    dist_df <- 
      dplyr::tibble(
        db_scientificName = NA,
        trait_out = trait,
        id = NA,
        tax_distance = NA,
        body_size_range_match = NA,
        life_stage_match = NA,
        explanation = "target name not found in any taxonomic backbone"
        
      )
    
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
              tax_distance = tax_dist,
              explanation = NA
            )
          
          return(dist_df)
          
        }, td[["scientificName"]], td[["id"]], SIMPLIFY = FALSE)
      
      # bind into a data.frame
      dist_df <- dplyr::bind_rows(dist_df)
      
      # get the body size range matches if the trait is an equation
      if (trait == "equation") {
        
        dist_df[["body_size_range_match"]] <- 
          
          sapply(dist_df[["id"]], function(x) {
            
            # extract the body size range match information
            bs_df <- extract_body_size_range_match(equation_id = x,
                                                   target_body_size = input[[body_size]],
                                                   equation_db = trait_db)
            
            return(bs_df)
            
          })
        
      }
      
      # get the life-stage matches
      dist_df[["life_stage_match"]] <- mapply(function(x, y) {
        if (!is.na(x)) {
          return((trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["db_life_stage"]] == y))
        } else {
          return(NA)
        }
      }, dist_df[["id"]], dist_df[[life_stage]])
      
      
      
      
  }
    
  }
    
    
    
    
    
    
    # check if the life-stages match
    
    
    
    if ( any(target_present == TRUE) & (body_size_filter == TRUE) ) {
      explan <- 
        ifelse(sum(dist_df[["body_size_range_match"]]) == 0,
               "no equation within max taxonomic distance has appropriate body size range",
               NA)
      dist_df <- dplyr::filter(dist_df, body_size_range_match == TRUE)
      
    }
    
    if( any(target_present == TRUE) & ( nrow(dist_df) > 0 ) ) {
      
      # remove the rows where the taxonomic distance is too great
      explan <- 
        ifelse(sum(dist_df[["tax_distance"]] <= max_tax_dist) == 0,
               "no equation within max taxonomic distance", 
               NA)
      dist_df <- dplyr::filter(dist_df, tax_distance <= max_tax_dist)
      
    }
    
    
    
    
  } 
  
  
  
   else if ( any(target_present == TRUE) ) {
    
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
            tax_distance = tax_dist,
            explanation = NA
          )
        
        return(dist_df)
        
      }, td[["scientificName"]], td[["id"]], SIMPLIFY = FALSE)
    
    # bind into a data.frame
    dist_df <- dplyr::bind_rows(dist_df)
    
  } else if ( all( target_present == FALSE ) ) {
    
    dist_df <- 
      dplyr::tibble(
        db_scientificName = NA,
        trait_out = trait,
        id = NA,
        tax_distance = NA,
        explanation = "target name not found in any taxonomic backbone"
        
      )
    
  }
  
  # get the body size range matches if the trait is an equation
  if (trait == "equation") {
    
    body_size_range_match <- 
      
      lapply(dist_df[["id"]], function(x) {
        
        # extract the body size range match information
        bs_df <- extract_body_size_range_match(equation_id = x,
                                               target_body_size = input[[body_size]],
                                               equation_db = trait_db)
        
        return(bs_df)
        
      })
    
    dist_df[["body_size_range_match"]] <- unlist(lapply(body_size_range_match, function(x) x[["body_size_range_match"]]))
    
  }
  
  if ( any(target_present == TRUE) & (body_size_filter == TRUE) ) {
    explan <- 
      ifelse(sum(dist_df[["body_size_range_match"]]) == 0,
             "no equation within max taxonomic distance has appropriate body size range",
             NA)
    dist_df <- dplyr::filter(dist_df, body_size_range_match == TRUE)
    
  }
  
  if( any(target_present == TRUE) & ( nrow(dist_df) > 0 ) ) {
    
    # remove the rows where the taxonomic distance is too great
    explan <- 
      ifelse(sum(dist_df[["tax_distance"]] <= max_tax_dist) == 0,
             "no equation within max taxonomic distance", 
             NA)
    dist_df <- dplyr::filter(dist_df, tax_distance <= max_tax_dist)
    
  }
  
  # if all equation are removed, we provide a data.frame and explanation
  if ( nrow(dist_df) == 0 ) {
    
    dist_df <- 
      dplyr::tibble(
        db_scientificName = NA,
        trait_out = trait,
        id = NA,
        tax_distance = NA,
        explanation = explan
      )
    
    if (trait == "equation") {
      dist_df[["body_size_range_match"]] <- NA
    }
    
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









