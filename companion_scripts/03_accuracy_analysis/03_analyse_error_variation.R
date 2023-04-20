
# generate a test database that uses many different types of equations

# load the relevant libraries
library(dplyr)

# load the use-scripts
source("companion_scripts/helper-plot-theme.R")

# load the required functions
source("R/clean_taxon_names.R")
source("R/get_habitat_data.R")
source("R/select_traits_tax_dist.R")
source("R/special_names.R")
source("R/helpers.R")

# load libraries required for those function
library(igraph)
library(assertthat)

# load the test data
dat <- readRDS("database/test_a_data_compilation.rds")

# set-up the relevant arguments
data <- dat
head(data)
target_taxon = "taxon"
life_stage = "life_stage"
latitude_dd = "lat"
longitude_dd = "lon"
body_size = "length_mm"
workflow = "workflow2"
max_tax_dist = 5
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
data_list <- split(clean_taxa, 1:nrow(clean_taxa))

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
    
    if (trait == "equation") {
      
      dist_df[["body_size_range_match"]] <- 
        
        sapply(dist_df[["id"]], function(x) {
          
          extract_body_size_range_match(equation_id = x, 
                                        target_body_size = input[[body_size]],
                                        equation_db = trait_db)
          
        })
      
    }
    
    # remove the rows where the taxonomic distance is too great
    dist_df <- dplyr::filter(dist_df, tax_distance <= max_tax_dist)
    
  } else {
      
      dist_df <- dplyr::tibble(
        db_scientificName = NA,
        trait_out = trait,
        id = NA,
        tax_distance = NA
      )
    
  }
  
  # add metadata
  dist_df <- dplyr::bind_cols(input, dist_df)
  
  return(dist_df)
  
  })

# if the algorithm selected the same equation with a different database
# we collapse these
output_df <- 
  bind_rows(output) %>%
  select(taxon_id, reference, order, taxon, lat, lon,
         life_stage, habitat_id, realm, major_habitat_type, ecoregion, tax_distance, body_size_range_match, id,
         db_scientificName, id, length_mm, obs_dry_biomass_mg) %>% 
  distinct()

  nrow(output_df)
  
  
  # run the select_traits_tax_dist() function: z1
  trait_sel <- select_traits_tax_dist(data = clean_taxa, 
                                      target_taxon = target_taxon,
                                      body_size = body_size,
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
  
  

  
  