# clean the equation data

# load relevant libraries
library(bdc)
library(stringdist)
library(dplyr)

# load special names function
source("R/special_names.R")

# load the equation data
equ_dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_FreshInvTraitR/data/allometry_database_ver4/equation_database.xlsx")
head(equ_dat)

# clean the names for typos etc.
x <- bdc::bdc_clean_names(sci_names = equ_dat$db_taxon, save_outputs = FALSE)

# check if any names were changed
if (!any(x$scientificName != x$names_clean)) {
  message("No names were changed")
}

# replace the names in tax.dat with these cleaned names
equ_dat$db_taxon <- x$names_clean

# fix the special names
spec_names <- special_taxon_names()

# replace incorrectly spelled special names
for (i in 1:length(spec_names)) {
  x <-
    sapply(equ_dat$db_taxon, function(y) {
      ain(x = spec_names[i], table = y, method = "lv", maxDist = 2)
    })

  equ_dat[x, "db_taxon"] <- spec_names[i]
}

# convert relevant columns to numeric variables

# maximum and minimum body size
equ_dat[["body_size_min"]] <- round(as.numeric(equ_dat[["body_size_min"]]), 4)
equ_dat[["body_size_max"]] <- round(as.numeric(equ_dat[["body_size_max"]]), 4)

# number of data points
equ_dat[["n"]] <- round(as.numeric(equ_dat[["n"]]), 0)

# r2 of the log-linear equation
equ_dat[["r2"]] <- round(as.numeric(equ_dat[["r2"]]), 2)

# maximum and minimum dry biomass
equ_dat[["dry_biomass_min"]] <- round(as.numeric(equ_dat[["dry_biomass_min"]]), 4)
equ_dat[["dry_biomass_max"]] <- round(as.numeric(equ_dat[["dry_biomass_max"]]), 4)

# residual mean squared error
equ_dat[["RMS"]] <- round(as.numeric(equ_dat[["RMS"]]), 4)

# back-transformation correction factor
equ_dat[["lm_correction"]] <- round(as.numeric(equ_dat[["lm_correction"]]), 4)

# preservation correction factor
equ_dat[["correction_percentage"]] <- round(as.numeric(equ_dat[["correction_percentage"]]), 4)

# convert the log-base to a numeric factor
equ_dat[["log_base"]] <- round(as.numeric(equ_dat[["log_base"]]), 5)

# convert the equation parameters
equ_dat[["a"]] <- round(as.numeric(equ_dat[["a"]]), 5)
equ_dat[["b"]] <- round(as.numeric(equ_dat[["b"]]), 5)

# calculate the correction factors

# function to calculate the BC3 correction factor sensu Strimbu et al. (2018)
BC_correction <- function(r2, a, ymin, ymax) {
  
  x <- exp( (0.5* (1 - r2)) * ( (1/(log(a))) * (( log(ymax, b=a)-log(ymin, b=a) )/6))^2)
  return(x)
}

# 1. calculate the BC-corrections
equ_dat <- 
  equ_dat |>
  dplyr::mutate(lm_correction = ifelse(lm_correction_type == "BC_correction",
                                       BC_correction(r2 = r2, 
                                                     a = log_base, 
                                                     ymin = dry_biomass_min, 
                                                     ymax = dry_biomass_max),
                                       lm_correction))

# 2. calculate the RMS_corrections
equ_dat <- 
  equ_dat |>
  dplyr::mutate(lm_correction = ifelse(lm_correction_type == "RMS_correction",
                                       log_base^(RMS/2),
                                       lm_correction))

# replace the character NAs with true NAs as interpreted by R
equ_dat[equ_dat == "NA"] <- NA

# write this into a .rds file
saveRDS(equ_dat, file = paste("database", "/", "equation_database.rds", sep = ""))
