
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

# split the input data.frame into a list
data_list <- split(data, 1:nrow(data))

# for each entry in the input.list, select appropriate traits
output <- lapply(data_list, function(input) {
  
  input <- data_list[[1]]
  
  # if the database is gbif, itis or col, determine if target name is present
  # in any of the higher taxonomic graphs
  if (input[["db"]] %in% c("gbif", "itis", "col")) {
    
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
    
  }
  
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
    
  } else {
    
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
      dist_df <- 
        dplyr::tibble(
          db_scientificName = x,
          trait_out = trait,
          id = y,
          tax_distance = NA,
          body_size_range_match = NA,
          life_stage_match = NA,
          explanation = NA
        ) 
      
    } else if any(target_present == TRUE) {
      
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
      
    }
    
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
    life_stage_match <- mapply(function(x, y) {
      if (!is.na(x)) {
        return((trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["db_life_stage"]] == y))
      } else {
        return(NA)
      }
    }, dist_df[["id"]], dist_df[[life_stage]])
    
    # add life-stage match column
    dist_df[["life_stage_match"]] <- life_stage_match
    
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









