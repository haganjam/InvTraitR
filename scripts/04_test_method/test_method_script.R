
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

test.output <- 
  Get_Trait_From_Taxon(input_data = data.test, 
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

# load the equation data
test.dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/test_data.xlsx")
test.dat$lat <- 50.45
test.dat$lon <- 5.266

# remove cases where life-stages are NA
test.dat <- dplyr::filter(test.dat, !is.na(Life_stage), Life_stage != "NA")
unique(test.dat$Life_stage)

# test the method
test.output <- 
  Get_Trait_From_Taxon(input_data = test.dat, 
                       target_taxon = "Taxa", 
                       life_stage = "Life_stage", 
                       body_size = "Length_mm",
                       latitude_dd = "lat", 
                       longitude_dd = "lon",
                       trait = "equation", 
                       max_tax_dist = 3,
                       gen_sp_dist = 0.5
  )

View(test.output)

# remove rows where the weight is not there
test.output <- 
  test.output %>%
  filter(!is.na(weight_mg) )

length(unique(test.output$Taxa))

library(ggplot2)
ggplot(data = test.output %>% 
         mutate(beyond_length_range = ifelse(is.na(flags), FALSE, TRUE)),
       mapping = aes(x = log(Dry_weight_mg), y = log(weight_mg)) ) +
  geom_point() +
  ylab("log estimated weight (mg)") +
  xlab("log actual weight (mg)") +
  geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "bottom")

# calculate error
test.output %>%
  mutate(error = (abs(Dry_weight_mg - weight_mg)/Dry_weight_mg)*100 ) %>%
  ggplot(data = .,
         mapping = aes(x = tax_distance, error)) +
  ylab("absolute error (%)") +
  xlab("taxonomic distance") +
  geom_jitter(width = 0.1) +
  theme_bw()

### END
