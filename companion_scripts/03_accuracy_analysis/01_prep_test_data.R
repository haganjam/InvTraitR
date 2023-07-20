
# clean the test data

# load relevant libaries
library(dplyr)
library(readr)

# set the seed because we do some random drawing
set.seed(3597655)

# test 1: prep

# load datasets compiled from the literature
test_names <- list.files("database/")
test_names <- test_names[grepl(pattern = "test_a_", x = test_names)]
test_names <- test_names[grepl(pattern = ".csv", x = test_names)]

# load the datasets and name them dat1 ... datN
dat_list <- vector("list", length = length(test_names))
for(i in 1:length(test_names)) {
  
  x <- readr::read_csv(paste0("database/", test_names[i]))
  dat_list[[i]] <- x
  
}

# rename the dat_list
names(dat_list) <- paste0("dat_", 1:length(test_names))

# check the datasets
lapply(dat_list, head)

# extract the relevant columns from each dataset
dat <- 
  lapply(dat_list, function(x) {
    
    x |>
      dplyr::select(reference, order, taxon, lat, lon, sex, gravid, life_stage, length_mm, dry_weight_mg, dry_weight_type) |>
      dplyr::rename(obs_dry_biomass_mg = dry_weight_mg)  
    
  })

# bind into a single large data.frame
dat <- dplyr::bind_rows(dat)

# how many species do we have here
length(unique(dat$taxon))

# remove cases where life-stages are NA
dat <- dplyr::filter(dat, !is.na(life_stage), life_stage != "NA")

# check the summary statistics
summary(dat)

# filter these negative values
dat <- 
  dat |>
  dplyr::filter( !(obs_dry_biomass_mg < 0) )

# remove the rows without lat and lon data
dat <- 
  dat |>
  dplyr::filter(!is.na(lat))

# remove gravid individuals
dat <- 
  dat |>
  dplyr::filter(gravid == "yes" | is.na(gravid))

# remove the Mahrlein datapoints
dat <- 
  dat |>
  dplyr::filter(reference != "Maehrlein_2016")

# sample from these data to make sure we don't pseudoreplicate too much
head(dat)
dat <- 
  dat |>
  dplyr::group_by(reference, taxon) |>
  dplyr::sample_n(size = ifelse(min(n()) < 5, min(n()), 5), replace = FALSE) |>
  dplyr::ungroup()

# how many samples are left?
nrow(dat)

# check the summary statistics
summary(dat)

# write out into an .rds file
saveRDS(dat, file = paste("database", "/", "test_a_data_compilation.rds", sep = ""))


# test 2: prep

# load dolmans data
dat_x <- read_csv("database/test_b_dolmans_2022.csv")
names(dat_x)

# select relevant columns
dat_x <- 
  dat_x |>
  dplyr::mutate(reference = "Dolmans 2022",
                lat_dd = NA,
                lon_dd = NA) |>
  dplyr::select(reference, Focal_taxon, Life_stage, lat_dd, lon_dd, length_mm, Biomass_mg)

# rename the columns
names(dat_x) <- c("reference", "taxon", "life_stage", "lat_dd", "lon_dd", "length_mm", "obs_dry_biomass_mg")

# load the o'gorman data
dat_y <- readr::read_csv("database/test_b_gorman_2017.csv")
names(dat_y)

# select relevant columns
dat_y <- 
  dat_y |>
  dplyr::mutate(reference = "OGorman 2017") |>
  dplyr::select(reference, species, life_stage, lat_dd, lon_dd, length, mass)

# rename the columns
names(dat_y) <- c("reference", "taxon", "life_stage", "lat_dd", "lon_dd", "length_mm", "obs_dry_biomass_mg")

# combine these two datasets
dat_z <- dplyr::bind_rows(dat_x, dat_y)

# write out into an .rds file
saveRDS(dat_z, file = paste("database", "/", "test_b_data_compilation.rds", sep = ""))

