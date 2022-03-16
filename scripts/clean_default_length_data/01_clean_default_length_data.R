
# Clean trait data

# load relevant libraries
library(here)
library(dplyr)
library(tidyr)
library(readr)


## Axell 2021:

# load the source list
axell.s <- read_delim(file = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/papers_to_get_test_data_from/Axell_2021/Database_references.txt", 
                      delim = "\t")

# check a random reference
axell.s[axell.s$`Source code` == 37, ]$Reference

# load the data
axell.d <- read_delim(file = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/papers_to_get_test_data_from/Axell_2021/Database_invertebrates.txt", 
                      delim = "\t")
head(axell.d)

# check parsing problems
problems(axell.d)

# fix the data loading problems
axell.d[1046, ]$`NFE (value)` <- 21.9
axell.d[1095, ]$`NFE (value)` <- 0.9

# check the added logical column
names(axell.d)

# the logical column is only FALSE values
any(!is.na(axell.d$...41))

# remove the ...41 variable which was a parsing failure
axell.d <- axell.d[, -41]

# subset the relevant rows of data
axell.lw <- 
  axell.d %>%
  filter( (!is.na(Genus)) & (!is.na(`TL (value)`)) ) %>%
  select(Genus, Species, LifeStage,
         `TL (value)`, `TL (reference)`,
         `DW (value)`, `DW (reference)`)

## TL variable
# split the TL (value) variable
x <- strsplit(axell.lw$`TL (value)`, split = "-") 

# add a min and max column
axell.lw$TL_min <- 
  unlist( 
    lapply(x, function(y) { 
      z <- gsub(pattern = "\\+/", replacement = "", y[1])
      gsub( pattern = "\\ .*", replacement = "", z) } )
  ) %>%
  as.numeric()

axell.lw$TL_max <- 
  unlist( 
    lapply(x, function(y) gsub(pattern = "]", replacement = "", y[2]) )
    ) %>%
  as.numeric()

# add a median column
axell.lw$TL_middle <- 
  ifelse(is.na(axell.lw$TL_max), axell.lw$TL_min, NA) %>%
  as.numeric()

# remove the duplicate minima 
axell.lw$TL_min <- 
  ifelse(!is.na(axell.lw$TL_middle), NA, axell.lw$TL_min) %>%
  as.numeric()

## DW variable 
# split the DW variable
x <- strsplit(axell.lw$`DW (value)`, split = "-") 

# create a middle variable
axell.lw$DW_middle <- 
  unlist( 
  lapply(x, function(y) ifelse(length(y) == 1, y[1], NA )) 
  ) %>%
  as.numeric()

# create a min variable
axell.lw$DW_min <- 
  unlist( 
    lapply(x, function(y) y[1])
    ) %>%
  as.numeric()

# create a max variable
axell.lw$DW_max <- 
  unlist( 
    lapply(x, function(y) y[2]) 
    ) %>%
  as.numeric()

# remove the duplicate minima 
axell.lw$DW_min <- 
  ifelse(!is.na(axell.lw$DW_middle), NA, axell.lw$DW_min) %>%
  as.numeric()

# add a database column
axell.lw$Database <- "Axell_2021"

## Reorder
# reorganise the columns
axell.lw <- 
  axell.lw %>%
  select(Database, Genus, Species, LifeStage, TL_reference = `TL (reference)`,
         TL_middle, TL_min, TL_max, DW_reference = `DW (reference)`,
         DW_middle, DW_min, DW_max)


## Tachet

# load the Tachet database
tach <- readxl::read_xlsx("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/freshwater_trait_databases/Tachet_2010_database/Tachet_invertebrate_database.xlsx",
                          skip = 3)
head(tach)

# subset the relevant columns
tach <- tach[, 3:14]

# subset more relevant columns
tach <- 
  tach %>%
  select(-...4, -Modalities, -`> 8 cm`)

# rename the columns
names(tach) <- c("Species", "Genus", "Taxon", "less_0.25_cm", "bet_0.25_0.5_cm",
                 "bet_0.5_1_cm", "bet_1_2_cm", "bet_2_4_cm", "bet_4_8_cm")

# add a Database column
tach$Database <- "Tachet_2010"

# add a LifeStage column
tach$LifeStage <- NA

# reorder the columns
tach <- tach[, c(10, 3, 2, 1, 11, 4:9)]

# remove rows with missing data
any(is.na(tach$Taxon))
tach <- 
  tach %>%
  filter(!is.na(Taxon))

# remove the crayfish species
tach <- 
  tach %>%
  filter(Genus != "Pacifastacus")

# calculate a midpoint of total length
r.seq <- list(c(0, 0.25), 
              c(0.25, 0.5),
              c(0.5, 1), 
              c(1, 2), 
              c(2, 4),
              c(4, 8) )

# get a sequence of length increments within each size class
r.seq <- lapply(r.seq, function(x) seq(x[1], x[2], 0.01))

# calculate the density of those increments assuming a uniform distribution within the size class
r.prob <- lapply(r.seq, function(x) dunif(x = x, min = min(x), max = max(x)) )

# assume that the density of each increments determines the frequency in the distribution
tach$TL_middle <- 
  
  apply(tach[, -c(1:5)], 1, function(z) { 
  
  a <- mapply(function(x, y) x*y, x1, z) 
  b <- rep(unlist(r.seq), round(unlist(a), 0)) 
  mean ( b[b>0] ) * 10 # convert to mm
  
  } )


## Default length data

# subset the tachet database
tach <- tach[, c(1, 3, 4, 5, 12)]
str(tach)

# subset the axell 2021 database
axell.def <- axell.lw[, c(1, 2, 3, 4, 6, 7, 8)]
str(axell.def)

# add TL_middle from the midpoint of the min_max values
axell.def <- 
  axell.def %>%
  mutate(TL_middle = if_else( is.na(TL_middle), (TL_min + TL_max)/2, TL_middle ) ) %>%
  select(Database, Genus, Species, LifeStage, TL_middle)
head(axell.def)

# bind these databases into a single default length data
def.length <- bind_rows(tach, axell.def)

# rename the TL_middle variable
def.length <- 
  def.length %>%
  rename(length_mid_mm = TL_middle)

# write this into a database
write_csv(x = def.length, 
          file = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/default_length_data.csv")

### END
