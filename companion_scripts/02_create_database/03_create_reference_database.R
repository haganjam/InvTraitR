#' @title reference list
#' @description create a data.frame of the equation references
#' @details This script compiles the list of publications from which the
#'   equations in the equation database were compiled
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 

# load relevant libraries
library(readxl)

# load the reference list
ref_dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_FreshInvTraitR/data/allometry_database_ver4/reference_database.xlsx")
head(ref_dat)

# check the database
View(ref_dat)

# replace the NA characters with true NAs
ref_dat[ref_dat == "NA"] <- NA

# export the database as a .rds file
saveRDS(ref_dat, "database/reference_database.rds")
