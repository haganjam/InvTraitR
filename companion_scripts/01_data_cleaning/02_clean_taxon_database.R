# Clean taxon database

# load relevant libraries
library(bdc)
library(stringdist)

# load special names function
source("R/special_names.R")

# load the equation data
t_dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_FreshInvTraitR/data/allometry_database_ver4/taxon_database.xlsx")
head(t_dat)

# clean the names for typos etc.
x <- bdc::bdc_clean_names(sci_names = t_dat$db_taxon, save_outputs = FALSE)

# check if any names were changed
if (!any(x$scientificName != x$names_clean)) {
  message("No names were changed")
}

# replace the names in tax.dat with these cleaned names
t_dat$db_taxon <- x$names_clean

# fix the special names
spec_names <- special_taxon_names()

# replace incorrectly spelled special names
for (i in 1:length(spec_names)) {
  x <-
    sapply(t_dat$db_taxon, function(y) {
      ain(x = spec_names[i], table = y, method = "lv", maxDist = 2)
    })

  t_dat[x, "db_taxon"] <- spec_names[i]
}

# write this into a .rds file
saveRDS(t_dat, file = paste("database", "/", "taxon_database.rds", sep = ""))
