
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

output <-
  get_trait_from_taxon(
    data = dat,
    target_taxon = "taxon",
    life_stage = "life_stage",
    body_size = "length_mm",
    latitude_dd = "lat",
    longitude_dd = "lon",
    workflow = "workflow1",
    trait = "equation",
    max_tax_dist = 4,
    gen_sp_dist = 0.5
  )

# set-up a vector to capture the dry biomass values
dry_biomass_mg <- vector(length = nrow(output))

# loop over all the rows
for(i in 1:nrow(output )) {
  
  # get the ith row of data
  L <- unlist(output[i, "length_mm"], use.names = FALSE)
  model <- output[i,]$equation_form
  log_base <- output[i,]$log_base
  a <- output[i,]$a
  b <- output[i,]$b
  CF <- output[i,]$lm_correction
  scale <- output[i,]$dry_biomass_scale
  
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
output[["dry_biomass_mg"]] <- dry_biomass_mg

# remove rows where the dry biomass estimate is less than 0.0001
output <- 
  output %>%
  filter(dry_biomass_mg > 0.0001)

# create a data.frame for modelling the error
output_df <- 
  output %>%
  select(row, order, taxon, db_scientificName,
         tax_distance, body_size_range_match, equation_form, r2_match, n,
         realm_match, major_habitat_type_match, ecoregion_match,
         length_mm, obs_dry_biomass_mg, dry_biomass_mg
         )

# get the percentage prediction error
output_df <- 
  output_df %>%
  mutate(error_perc = ((obs_dry_biomass_mg - dry_biomass_mg)/obs_dry_biomass_mg)*100,
         abs_error_perc = (abs(obs_dry_biomass_mg - dry_biomass_mg)/obs_dry_biomass_mg)*100) 


# get a habitat match variable
output_df[["habitat_match"]] <- 
  apply(output_df[, paste0(c("realm", "major_habitat_type", "ecoregion"), "_match") ], 1,
        function(x) sum(x) )

# absolute error

# check the distribution of these errors
hist(output_df$abs_error_perc)
summary(output_df$abs_error_perc)

# what about the log-transformed abs error perc
hist(log(output_df$abs_error_perc))

# error

# check the distribution
hist(output_df$error_perc)
summary(output_df$error_perc)

# how many data points do we have per taxon id
output_df %>%
  group_by(row) %>%
  summarise( n = n()) %>%
  pull(n)

# check the most extreme data points
output_df %>%
  filter(error_perc < -300) %>%
  View()

# does the model type make a difference
output_df %>%
  filter(order == "Bivalvia") %>%
  group_by(equation_form) %>%
  summarise(mean_error_perc = mean(abs_error_perc))

# fit a massive interaction model
mod_dat <- 
  output_df %>%
  select(row, 
         taxon,
         tax_distance, body_size_range_match, habitat_match,
         r2_match, n, 
         length_mm, dry_biomass_mg,
         error_perc) %>%
  mutate(row = as.character(row),
         log_error_perc = log(1700 + error_perc) )

lm1 <- lm(error_perc ~ 
            row + 
            row:tax_distance + 
            row:body_size_range_match +
            row:r2_match,
          data = mod_dat)
summary(lm1)
anova(lm1)

# plot the observed versus the predicted values
plot(mod_dat$error_perc, predict(lm1))
abline(0, 1)





  
  