#' @title preservation correction factors
#' @description create a data.frame with the preservation correction factor data
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 

# load relevant libraries
library(readxl)

# load the reference list
pre_dat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_FreshInvTraitR/data/allometry_database_ver4/dry_biomass_correction_data.xlsx")
head(pre_dat)

# check the database
View(pre_dat)

# replace the NA characters with true NAs
pre_dat[pre_dat == "NA"] <- NA

# export the database as a .rds file
saveRDS(pre_dat, "database/preservation_correction_database.rds")
