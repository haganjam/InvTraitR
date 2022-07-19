
# Clean data from Hebert et al. 2016, Ecology
# https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/15-1275.1#support-information-section

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(here)

# load the relevant data
z.dat <- 
  read_csv2(file = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/input_data/papers_to_get_test_data_from/Hebert_2016a/zooplankton_traits.csv")
str(z.dat)
View(z.dat)

# subset the freshwater species
z.dat <- 
  z.dat %>%
  filter(Habitat == "Freshwater") %>%
  select(Genus, Species, Body.length, Dry.mass)

# get only the complete.cases()
z.dat <- z.dat[complete.cases(z.dat),]

# remove the sp. values in the species column
z.dat <- 
  z.dat %>%
  mutate(Species = ifelse(Species == "sp.", "", Species))
View(z.dat)

# make a taxa name column
z.dat <- 
  z.dat %>%
  mutate(Taxa = paste(Genus, Species, sep = " ")) %>%
  mutate(Taxa = trimws(Taxa, "right"))

# reorder the columns
z.dat <- 
  z.dat %>%
  select(Taxa, Body.length, Dry.mass) %>%
  rename(Length_mm = Body.length, Dry_weight_mg = Dry.mass)

# write this as a .csv file
write_csv(x = z.dat,
          file = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/papers_to_get_test_data_from/Hebert_2016a/zooplankton_traits_test.csv")

### END
