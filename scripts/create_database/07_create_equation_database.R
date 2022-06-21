
# Create a merged data.list from the equation data and variable input databases

# load relevant libraries
library(dplyr)
library(here)

# read in the equation data
equ.dat <- readxl::read_xlsx("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/equation_data.xlsx")
equ.dat <- equ.dat[!is.na(equ.dat$id),]
equ.dat[equ.dat == "NA"] <- NA

# read in the variable inputs
in.dat <- readxl::read_xlsx("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/variable_input_data.xlsx")
in.dat <- in.dat[!is.na(in.dat$id),]
in.dat[in.dat == "NA"] <- NA

# list of equation IDs with only length data
x <- aggregate(in.dat$id, by = list(in.dat$id), length, simplify = TRUE)
y <- in.dat[in.dat$id %in% x[x$x == 1, ]$Group.1, ]
id.length <- y[y$size_measurement == "body_length", ]$id

# pull these data into a list
equ_vars <- 
  list(equation_data = as_tibble(equ.dat),
       variable_input_data = as_tibble(in.dat),
       id_only_equ_ID = id.length)

# merge into a list
saveRDS(equ_vars, file = here("database/equation_vars_database.rds") )


# Create a database of default lengths
dl.dat <- readr::read_csv("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/default_length_data.csv")
head(dl.dat)

# load dplyr
library(dplyr)

# add a target taxon column
dl.dat <- 
  dl.dat %>%
  rename(db_taxon = Taxon)

# add a length id column
dl.dat$length_id <- 1:nrow(dl.dat)

# reorder the columns
dl.dat <- 
  dl.dat %>%
  select(length_id, db_taxon, Database, LifeStage, length_mid_mm)

# rename the columns
dl.dat <- 
  dl.dat %>%
  rename(id = length_id, life_stage = LifeStage)

# replace NA characters with true NAs
dl.dat[dl.dat == "NA"] <- NA

# save as a .RDS file
saveRDS(dl.dat, file = here("database/default_length_database.rds") )

### END
