
# Test the method for calculating species biomasses

# load the relevant libraries
library(here)
library(dplyr)
library(readr)
library(ggplot2)

# load the use-scripts
source(here("scripts/03_use_database/01_search_database.R"))


# test 1: actual body size data and actual length data from the literature

# load the literature test data
test1 <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/test_data.xlsx")
test1$lat <- 50.45
test1$lon <- 5.266

# remove cases where life-stages are NA
test1 <- dplyr::filter(test1, !is.na(Life_stage), Life_stage != "NA")
unique(test1$Life_stage)

# test the method
test1.output <- 
  Get_Trait_From_Taxon(input_data = test1, 
                       target_taxon = "Taxa", 
                       life_stage = "Life_stage", 
                       body_size = "Length_mm",
                       latitude_dd = "lat", 
                       longitude_dd = "lon",
                       trait = "equation", 
                       max_tax_dist = 3,
                       gen_sp_dist = 0.5
  )

View(test1.output)

# remove rows where the weight is not there
test1.output <- 
  test1.output %>%
  filter(!is.na(weight_mg) )

length(unique(test1.output$Taxa))

library(ggplot2)
ggplot(data = test1.output %>% 
         mutate(beyond_length_range = ifelse(is.na(flags), FALSE, TRUE)),
       mapping = aes(x = log(Dry_weight_mg), y = log(weight_mg)) ) +
  geom_point() +
  ylab("log estimated weight (mg)") +
  xlab("log actual weight (mg)") +
  geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "bottom")

# calculate error
test1.output %>%
  mutate(error = (abs(Dry_weight_mg - weight_mg)/Dry_weight_mg)*100 ) %>%
  ggplot(data = .,
         mapping = aes(x = tax_distance, error)) +
  ylab("absolute error (%)") +
  xlab("taxonomic distance") +
  geom_jitter(width = 0.1) +
  theme_bw()


# test 2: equations selected by expert

# load data compiled by Vincent
test2 <- read_csv("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/test_data_vincent.csv")
test2$lat <- NA
test2$lon <- NA

# test the method
test2.output <- 
  Get_Trait_From_Taxon(input_data = test2, 
                       target_taxon = "Focal_taxon", 
                       life_stage = "Life_stage", 
                       body_size = "length_mm",
                       latitude_dd = "lat", 
                       longitude_dd = "lon",
                       trait = "equation", 
                       max_tax_dist = 3.5,
                       gen_sp_dist = 0.5
  )

# remove rows where the weight is not there
test2.output <- 
  test2.output %>%
  filter(!is.na(weight_mg) )

length(unique(test2.output$Focal_taxon))

library(ggplot2)
ggplot(data = test2.output %>% 
         mutate(beyond_length_range = ifelse(is.na(flags), FALSE, TRUE)),
       mapping = aes(x = log(Biomass_mg), y = log(weight_mg)) ) +
  geom_point() +
  ylab("log estimated weight (mg)") +
  xlab("log Vincent's weight (mg)") +
  geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "bottom")

### END
