
# Test the method for calculating species biomasses

library(here)
library(dplyr)
library(readr)
library(ggplot2)
source(here("scripts/03_use_database/01_search_database.R"))

data.test <- data.frame(name = c("Aedes",
                                 "Daphnia pulex", NA, "Daphnia"),
                        life_stage = c("larva", "adult", NA, "adult"),
                        lat = c(57.5, 57.5, 60, NA),
                        lon = c(119.3, 112.8, -60, NA),
                        length = c(4, 8, 10, 20)
                        )

# test the data by getting biomass conversion data for the Inselberg and Korranneberg data
dat.loc <- "C:/Users/james/Documents/github/predicting_trait_responses/data/biomass_conversions/"

# Korranneberg data
kor <- read_csv(file = paste(dat.loc, "kor_bio.csv", sep = "") )
View(kor)

# remove the missing data
kor <- 
  kor %>%
  filter(!is.na(life_stage))
View(kor)

# add biomass data
kor$body_mm <- NA

# add latitude longitude data
kor$lat <- -28.87
kor$lon <- 27.21

kor.output <- 
  Get_Trait_From_Taxon(input_data = kor, 
                       target_taxon = "taxon", 
                       life_stage = "life_stage", 
                       body_size = "body_mm",
                       latitude_dd = "lat", 
                       longitude_dd = "lon",
                       trait = "equation", 
                       max_tax_dist = 3,
                       gen_sp_dist = 0.5
                       )

plot(kor.output$biomass_mg, kor.output$weight_mg)

### END
