
# Test the method for calculating species biomasses

# load the relevant libraries
library(here)
library(dplyr)
library(readr)
library(ggplot2)

# load the use-scripts
source(here("scripts/03_use_database/01_search_database.R"))
source(here("scripts/02_function_plotting_theme.R"))


# test 1: actual body size data and actual length data from the literature

# load the literature test data
test1 <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/test_data.xlsx")
str(test1)

# convert the lat-lon variables to numeric variables
test1$lat <- as.numeric(test1$lat)
test1$lon <- as.numeric(test1$lon)

# remove cases where life-stages are NA
test1 <- dplyr::filter(test1, !is.na(Life_stage), Life_stage != "NA")

# test the method
test1.output <- 
  Get_Trait_From_Taxon(input_data = test1, 
                       target_taxon = "Taxa", 
                       life_stage = "Life_stage", 
                       body_size = "Length_mm",
                       latitude_dd = "lat", 
                       longitude_dd = "lon",
                       trait = "equation", 
                       max_tax_dist = 3.5,
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
  mutate(flags = as.character(if_else(is.na(flags), 0, 1)) ) %>%
  ggplot(data = .,
         mapping = aes(x = tax_distance, y = error, colour = flags)) +
  ylab("absolute error (%)") +
  xlab("taxonomic distance") +
  geom_jitter(width = 0.1) +
  theme_bw()

# view the outliers
test1.output %>%
  mutate(error = (abs(Dry_weight_mg - weight_mg)/Dry_weight_mg)*100 ) %>%
  filter(error > 1000) %>%
  View()

# replot without the major outliers
test1.output %>%
  mutate(error = ((Dry_weight_mg - weight_mg)/Dry_weight_mg)*100 ) %>%
  mutate(flags = as.character(if_else(is.na(flags), 0, 1)) ) %>%
  filter(error < 500 & error > -500) %>%
  ggplot(data = .,
         mapping = aes(x = tax_distance, y = error, colour = flags)) +
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
  xlab("log expert weight (mg)") +
  geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") +
  theme_meta() +
  theme(legend.position = "bottom")

### END
