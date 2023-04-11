# clean the equation data

# load relevant libraries
library(bdc)
library(stringdist)

# load special names function
source("R/special_names.R")

# load the equation data
equ.dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_FreshInvTraitR/data/allometry_database_ver3/equation_database.xlsx")
head(equ.dat)

# clean the names for typos etc.
x <- bdc_clean_names(sci_names = equ.dat$db_taxon, save_outputs = FALSE)

# check if any names were changed
if (!any(x$scientificName != x$names_clean)) {
  message("No names were changed")
}

# replace the names in tax.dat with these cleaned names
equ.dat$db_taxon <- x$names_clean

# fix the special names
spec.names <- special_taxon_names()

# replace incorrectly spelled special names
for (i in 1:length(spec.names)) {
  x <-
    sapply(equ.dat$db_taxon, function(y) {
      ain(x = spec.names[i], table = y, method = "lv", maxDist = 2)
    })

  equ.dat[x, "db_taxon"] <- spec.names[i]
}

# convert relevant columns to numeric variables

# maximum and minimum body size
equ.dat[["body_size_min"]] <- round(as.numeric(equ.dat[["body_size_min"]]), 4)
equ.dat[["body_size_max"]] <- round(as.numeric(equ.dat[["body_size_max"]]), 4)

# number of data points
equ.dat[["n"]] <- round(as.numeric(equ.dat[["n"]]), 0)

# r2 of the log-linear equation
equ.dat[["r2"]] <- round(as.numeric(equ.dat[["r2"]]), 2)

# maximum and minimum dry biomass
equ.dat[["dry_biomass_min"]] <- round(as.numeric(equ.dat[["dry_biomass_min"]]), 4)
equ.dat[["dry_biomass_max"]] <- round(as.numeric(equ.dat[["dry_biomass_max"]]), 4)

# residual mean squared error
equ.dat[["RMS"]] <- round(as.numeric(equ.dat[["RMS"]]), 4)

# back-transformation correction factor
equ.dat[["lm_correction"]] <- round(as.numeric(equ.dat[["lm_correction"]]), 4)

# check the variables structure
str(equ.dat)

# write this into a .rds file
saveRDS(equ.dat, file = paste("database", "/", "equation_database.rds", sep = ""))
