
# Search taxon database

# next steps:

# run all of these steps for each species in our equation database for
# all three databases
# at each run, we must add the equation id
# thus, when there is a match, we can immediately get to the relevant dataset
# and calculate the taxonomic proximity

# to do:

# before implementing this function, we can have a list of all names in the database
# this will then mean that before anyone puts a name in, we can test if it is in the database
# this means we don't have to search them all sequentially

# for implementing the functions with not just weights but also lengths
# ask for dataframes with species names and each measurement
# then the function can merge them so it fills in NA's when a measurement is not necessary

in.dat <- readxl::read_xlsx(here("raw_data/variable_input_data.xlsx"))
in.dat <- in.dat[!is.na(in.dat$equation_id),]

## some possible search function options

# taxon name database
db.name <- "Pediciini"
db.name

# target name
tar.name <- "Limnophila costata"
tar.name

# this calculation across all the different distance matrices
d[which(row.names(d) == db.name), which(colnames(d) == tar.name) ]

# rank of all the different names
