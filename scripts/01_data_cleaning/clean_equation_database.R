
# clean the equation data

# load relevant libraries
library(bdc)
library(here)

# load the equation data
equ.dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/equation_database.xlsx")
head(equ.dat)
View(equ.dat)

# clean the names for typos etc.
x <- bdc_clean_names(sci_names = equ.dat$db_taxon, save_outputs = FALSE)

# check if any names were changed
if ( !any(x$scientificName != x$names_clean) ) {
  message("No names were changed")
}

# replace the names in tax.dat with these cleaned names
equ.dat$db_taxon <- x$names_clean

# write some code to remove the output file
unlink("Output", recursive=TRUE)

# write this into a .rds file
saveRDS(equ.dat, file = paste(here("database"), "/", "equation_database.rds", sep = ""))

### END
