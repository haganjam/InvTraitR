
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

# rename length_id to id
dl.dat <-
  dl.dat %>%
  rename(id = length_id)

# save as a .RDS file
saveRDS(dl.dat, file = here("database/default_length_database.rds") )


# Supplementary equation and default length data

# read the equation data
equ_supp <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/equation_data_supp.xlsx")

# read the variable input data
equ_var_supp <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/variable_input_data_supp.xlsx")

# join these data into a list and save this as a .rds file
equ_vars_supp <- 
  list(equation_data = as_tibble(equ_supp),
       variable_input_data = as_tibble(equ_var_supp))

# merge into a list
saveRDS(equ_vars_supp, file = here("database/equation_vars_supp_database.rds") )


# read the default length data supplementary
len_supp <- readr::read_csv("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/default_length_data_supp.csv")

# save as a .RDS file
saveRDS(len_supp, file = here("database/default_length_supp_database.rds") )

### END
