
# Test the method for calculating species biomasses

library(here)
library(dplyr)
library(readr)
library(ggplot2)
source(here("scripts/use_database/01_search_equ_len_database_functions.R"))

# test the data by getting biomass conversion data for the Inselberg and Korranneberg data
dat.loc <- "C:/Users/james/Documents/github/predicting_trait_responses/data/biomass_conversions/"

# Korranneberg data
kor <- read_csv(file = paste(dat.loc, "kor_bio.csv", sep = "") )

# use the method to get biomass data using default length data
kor$length_dat <- NA
sort(unique(kor$taxon))

View(kor)

# add life-stage data
c("")

kor.x <- 
  get_taxa_mass(data.base = "itis",
                max_tax_dist = 6,
                data = kor,
                target.name.col = "taxon",
                life.stage.col = "life_stage",
                length.col = "length_dat" )

View(kor.x)

# write this into a .csv file that can be supplemented
write_csv(x = kor.x, file = paste(dat.loc, "kor_bio_input.csv", sep = ""))

# Inselberg data
ins <- read_csv(file = paste(dat.loc, "aus_ins_bio.csv", sep = "") )

# use the method to get biomass data using default length data
ins$length_dat <- NA
sort(unique(ins$taxon))

ins.x <- 
  get_taxa_mass(data.base = "itis",
                max_tax_dist = 8,
                data = ins,
                target.name.col = "taxon",
                life.stage.col = "life_stage",
                length.col = "length_dat" )

View(ins.x$mass_data)

# write this into a .csv file that can be supplemented
write_csv(x = ins.x$mass_data, file = paste(dat.loc, "ins_aus_input.csv", sep = ""))

### END
