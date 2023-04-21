
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
max_tax_dist = 4
trait = "equation"
gen_sp_dist = 0.5

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

# add a row_id variable
clean_taxa <- dplyr::bind_cols(
  dplyr::tibble(row_id = 1:nrow(clean_taxa)),
  clean_taxa
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
data_list <- split(clean_taxa, clean_taxa$row_id)

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

# remove rows without equation ids
output_df <- 
  output_df %>%
  filter(!is.na(id))

# life_stage match
output_df[["life_stage_match"]] <- 
  
  mapply(function(x, y) {
    if (!is.na(x)) {
      return((trait_db[trait_db[[paste0(trait, "_id")]] == x, ][["db_life_stage"]] == y))
    } else {
      return(NA)
    }
  }, output_df[["id"]], output_df[[life_stage]])
  
# only keep rows with matching life-stage information
output_df <- 
  output_df %>%
  filter(life_stage_match == TRUE)
  
# additional matches that are only relevant for the equation trait

# set-up a vector of the relevant columns and relevant names
rel_cols <- c("r2", "n", "body_size_min", "body_size_max")
rel_names <- c("r2_match", "n", "db_min_body_size_mm", "db_max_body_size_mm")

# loop over these variables
for (i in 1:length(rel_cols)) {
  
  output_df[[rel_names[i]]] <- 
    sapply(output_df[["id"]], function(x) {
      if (!is.na(x)) {
        return(trait_db[trait_db[[paste0(trait, "_id")]] == x, ][[rel_cols[i]]])
      } else {
        return(NA)
      }
    })
}

# add habitat match data
# set-up a vector of the relevant columns and relevant names
hab_cols <- c("realm", "major_habitat_type", "ecoregion")
hab_names <- paste0(hab_cols, "_match")

# load the habitat database
hab_db <- readRDS("database/freshwater_ecoregion_data.rds")

# loop over these habitat variables
for (i in 1:length(hab_cols)) {
  
  output_df[[hab_names[i]]] <-
    mapply(function(x, y) {
      if (!is.na(x)) {
        return((hab_db[hab_db[["id"]] == x, ][[hab_cols[i]]] == y))
      } else {
        return(NA)
      }
    }, output_df[["id"]], output_df[[hab_cols[i]]])
  
}

# add the correction factor data
# extract the relevant columns
cor_factors <- trait_db[, c("equation_id", 
                            "preservation",
                            "equation_form", "log_base", "a", "b",  
                            "lm_correction", "lm_correction_type",
                            "dry_biomass_scale")]
names(cor_factors)[1] <- "id"

# join these columns to the trait_sel_select
output_df <- dplyr::left_join(output_df, cor_factors, by = "id")

# set-up a vector to capture the dry biomass values
dry_biomass_mg <- vector(length = nrow(output_df))

# loop over all the rows
for(i in 1:nrow(output_df)) {
  
  # get the ith row of data
  L <- unlist(output_df[i, body_size], use.names = FALSE)
  model <- output_df[i,]$equation_form
  log_base <- output_df[i,]$log_base
  a <- output_df[i,]$a
  b <- output_df[i,]$b
  CF <- output_df[i,]$lm_correction
  scale <- output_df[i,]$dry_biomass_scale
  
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

# add the dry biomass mg column to the data
output_df[["dry_biomass_mg"]] <- dry_biomass_mg

# how many taxon id datapoints do we have compared to all datapoints
length(unique(output_df$taxon_id))
nrow(output_df)

# create a data.frame for modelling the error
output_lm <- 
  output_df %>%
  select(taxon_id, order, taxon, db_scientificName,
         tax_distance, body_size_range_match, r2_match, n,
         realm_match, major_habitat_type_match, ecoregion_match,
         length_mm, obs_dry_biomass_mg, dry_biomass_mg
         )

# get the percentage prediction error
output_lm <- 
  output_lm %>%
  mutate(error_perc = ((obs_dry_biomass_mg - dry_biomass_mg)/obs_dry_biomass_mg)*100,
         abs_error_perc = (abs(obs_dry_biomass_mg - dry_biomass_mg)/obs_dry_biomass_mg)*100) 


# get a habitat match variable
output_lm[["habitat_match"]] <- 
  apply(output_lm[, paste0(c("realm", "major_habitat_type", "ecoregion"), "_match") ], 1,
        function(x) sum(x) )

# absolute error

# check the distribution of these errors
hist(output_lm$abs_error_perc)
summary(output_lm$abs_error_perc)

# what about the log-transformed abs error perc
hist(log(output_lm$abs_error_perc))

# error

# check the distribution
hist(output_lm$error_perc)
summary(output_lm$error_perc)

# how many data points do we have per taxon id
output_lm %>%
  group_by(taxon_id) %>%
  summarise( n = n()) %>%
  pull(n)

# fit a massive interaction model

mod_dat <- 
  output_lm %>%
  select(taxon_id, 
         taxon,
         tax_distance, body_size_range_match, habitat_match,
         r2_match, n, 
         length_mm, dry_biomass_mg,
         abs_error_perc) %>%
  mutate(taxon_id = as.character(taxon_id))

lm1 <- lm(log(abs_error_perc) ~ 
            taxon_id + 
            taxon_id:tax_distance + 
            taxon_id:body_size_range_match +
            taxon_id:r2_match,
          data = mod_dat)
summary(lm1)
anova(lm1)

x <- lm1$coefficients[grepl("body_size_range_match", names(lm1$coefficients) )]
hist(x)
summary(x)

lm1$coefficients[lm1$coefficients > 40]




  
  