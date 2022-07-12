
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

data.test2 <- data.test
data.test2$lat <- NA
data.test2$lon <- NA

x <- Get_Habitat_Data(data = data.test2, latitude_dd = "lat", longitude_dd = "lon")
View(x)

test.output <- 
  Get_Trait_From_Taxon(input_data = data.test2, 
                       target_taxon = "name", 
                       life_stage = "life_stage", 
                       body_size = "length",
                       latitude_dd = "lat", 
                       longitude_dd = "lon",
                       trait = "equation", 
                       max_tax_dist = 3,
                       gen_sp_dist = 0.5
  )

View(test.output)




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


# load the equation data
test.dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/test_data.xlsx")
View(test.dat)
test.dat$lat <- -28.87
test.dat$lat <- as.numeric(test.dat$lat)
test.dat$lon <- 27.21
test.dat$lon <- as.numeric(test.dat$lon)
library(dplyr)
test.dat <- 
  test.dat %>%
  filter(Reference == "Dumont_1975")

test.dat <- dplyr::filter(test.dat, !is.na(Life_stage), Life_stage != "NA")
unique(test.dat$Life_stage)

test.dat <- test.dat[sample(1:nrow(test.dat), 20), ]
View(test.dat)

test.output <- 
  Get_Trait_From_Taxon(input_data = test.dat, 
                       target_taxon = "Taxa", 
                       life_stage = "Life_stage", 
                       body_size = "Length_mm",
                       latitude_dd = "lat", 
                       longitude_dd = "lon",
                       trait = "equation", 
                       max_tax_dist = 4,
                       gen_sp_dist = 0.5
  )

View(test.output)

### END
