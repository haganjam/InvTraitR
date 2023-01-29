# Clean taxon database

# load relevant libraries
library(bdc)
library(here)
library(stringdist)

# load special names function
source(here("R/01_special_names_func.R"))

# load the equation data
t.dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/taxon_database.xlsx")
head(t.dat)

# clean the names for typos etc.
x <- bdc_clean_names(sci_names = t.dat$db_taxon, save_outputs = FALSE)

# check if any names were changed
if (!any(x$scientificName != x$names_clean)) {
  message("No names were changed")
}

# replace the names in tax.dat with these cleaned names
t.dat$db_taxon <- x$names_clean

# write some code to remove the output file
unlink("Output", recursive = TRUE)

# fix the special names
spec.names <- special_taxon_names()

# replace incorrectly spelled special names
for (i in 1:length(spec.names)) {
  x <-
    sapply(t.dat$db_taxon, function(y) {
      ain(x = spec.names[i], table = y, method = "lv", maxDist = 2)
    })

  t.dat[x, "db_taxon"] <- spec.names[i]
}

# write this into a .rds file
saveRDS(t.dat, file = paste(here("database"), "/", "taxon_database.rds", sep = ""))
