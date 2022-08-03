
# Clean data from Gorman et al. 2017

# load relevant libraries
library(readr)
library(dplyr)
library(ggplot2)

# set the seed
set.seed(549854)

# load the biomass data
go.bio <- read_delim("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/input_data/papers_to_get_test_data_from/Gorman_2017/test_data_Hengill_2008_raw.csv",
                     delim = ",")
head(go.bio)

# remove data that uses hw
go.bio <- go.bio[!grepl(pattern = "HW", x = go.bio$formula), ]

# sample 10 from each stream from each taxa as having too much replication just
# leads to unnecesary extra information
go.bio <- 
  go.bio %>%
  group_by(species) %>%
  sample_n(size = 15, replace = TRUE) %>%
  distinct()

# output as a .csv file
write_csv(go.bio, "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/input_data/papers_to_get_test_data_from/Gorman_2017/test_data_Gorman_2017.csv")

# life-stage data was then manually added to this dataset

### END
