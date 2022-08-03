
# Clean Vincent's biomass conversions

# load relevant libraries
library(readr)
library(dplyr)
library(ggplot2)

# load the biomass data
vd.bio <- read_delim("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/input_data/papers_to_get_test_data_from/Vincent_2022/test_data_Vincent_biomass.csv",
                     delim = ";")
head(vd.bio)

# rename the taxon column
vd.bio <- 
  vd.bio %>%
  rename(Focal_taxon = ...4) %>%
  filter(!is.na(Focal_taxon)) %>%
  select(-starts_with("..."), -`Formulae of`) %>%
  rename(Taxon = `Length of `)
View(vd.bio)

# load the length data
vd.len <- read_delim("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/input_data/papers_to_get_test_data_from/Vincent_2022/test_data_Vincent_lengths.csv",
                     delim = ";")
head(vd.len)
View(vd.len)

# remove any length data that are not BL
vd.len <- 
  vd.len %>%
  filter(Kenmerk == "BL")

# rename the notes column
vd.len <- 
  vd.len %>%
  rename(notes = ...5)

# remove the additional, empty columns
vd.len <- 
  vd.len %>%
  select(-starts_with("..."))
View(vd.len)

# left join the length and biomass data
vd <- inner_join(vd.bio, vd.len, by = "Taxon")

# move the columns around
names(vd)
vd <- 
  vd %>%
  select(Focal_taxon, Taxon, Locatie, Kenmerk, `Lengte (mm)`, Biomass) %>%
  rename(Length_location = Locatie, Length_taxon = Taxon, Measurement = Kenmerk, length_mm = `Lengte (mm)`,
         Biomass_mg = Biomass)

# convert the length_mm data to a numeric variable
vd$length_mm <- as.numeric(vd$length_mm)
View(vd)
str(vd)

# two entries with erroneously large lengths
# Paradiaptomus
# Diaptomidae
vd[vd$Focal_taxon %in% c("Diaptomidae", "Paradiaptomus"), "length_mm"] <- c(0.336, 0.37)

write_csv(x = vd,
          "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/test_data_vincent.csv")

### END
