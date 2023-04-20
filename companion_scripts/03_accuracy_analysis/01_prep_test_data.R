
# clean the test data

# load relevant libaries
library(dplyr)
library(readr)

# load datasets compiled from the literature
test_names <- list.files("database/")
test_names <- test_names[grepl(pattern = "test_a_", x = test_names)]

# load the datasets and name them dat1 ... datN
dat_list <- vector("list", length = length(test_names))
for(i in 1:length(test_names)) {
  
  x <- read_csv(paste0("database/", test_names[i]))
  dat_list[[i]] <- x
  
}

# rename the dat_list
names(dat_list) <- paste0("dat_", 1:length(test_names))

# check the datasets
lapply(dat_list, head)

# extract the relevant columns from each dataset
dat <- 
  lapply(dat_list, function(x) {
    
    x %>%
      dplyr::select(reference, order, taxon, lat, lon, sex, gravid, life_stage, length_mm, dry_weight_mg, dry_weight_type) %>%
      rename(obs_dry_biomass_mg = dry_weight_mg)  
    
  })

# bind into a single large data.frame
dat <- bind_rows(dat)

# how many species do we have here
length(unique(dat$taxon))

# remove cases where life-stages are NA
dat <- dplyr::filter(dat, !is.na(life_stage), life_stage != "NA")

# check the summary statistics
summary(dat)

# filter these negative values
dat <- 
  dat %>%
  filter( !(obs_dry_biomass_mg < 0) )

# remove the rows without lat and lon data
dat <- 
  dat %>%
  filter(!is.na(lat))

# remove gravid individuals
dat <- 
  dat %>%
  filter(gravid == "yes" | is.na(gravid))

# remove the Mahrlein datapoints
dat <- 
  dat %>%
  filter(reference != "Maehrlein_2016")

# sample from these data to make sure we don't pseudoreplicate too much
head(dat)
dat <- 
  dat %>%
  group_by(reference, taxon) %>%
  sample_n(size = ifelse(min(n()) < 5, min(n()), 5), replace = FALSE) %>%
  ungroup()

# how many samples are left?
nrow(dat)

# check the summary statistics
summary(dat)

# write out into an .rds file
saveRDS(dat, file = paste("database", "/", "test_a_data_compilation.rds", sep = ""))
