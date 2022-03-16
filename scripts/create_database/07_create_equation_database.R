
# Create a merged data.list from the equation data and variable input databases

# read in the equation data
equ.dat <- readxl::read_xlsx("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/equation_data.xlsx")
equ.dat <- equ.dat[!is.na(equ.dat$equation_id),]

# read in the variable inputs
in.dat <- readxl::read_xlsx("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/variable_input_data.xlsx")
in.dat <- in.dat[!is.na(in.dat$equation_id),]

# list of equation IDs with only length data
x <- aggregate(in.dat$equation_id, by = list(in.dat$equation_id), length, simplify = TRUE)
y <- in.dat[in.dat$equation_id %in% x[x$x == 1, ]$Group.1, ]
id.length <- y[y$size_measurement == "body_length", ]$equation_id

# pull these data into a list
equ_vars <- 
  list(equation_data = equ.dat,
       variable_input_data = in.dat,
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
  mutate(target_taxon = if_else(is.na(Species), Genus, Species) )

# save as a .RDS file
saveRDS(dl.dat, file = here("database/default_length_database.rds") )

### END
