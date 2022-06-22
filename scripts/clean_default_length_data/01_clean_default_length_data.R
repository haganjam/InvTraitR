
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

# add life stage information to the Axell 2021 database
ls_df <- 
  data.frame(Species = c("Erpobdella octoculata",
                       "Glossiphonia complanata",
                       "Helobdella stagnalis",
                       "Aulodrilus pigueti",
                       "Limnodrilus",
                       "Streptocephalus sealii",
                       "Acantholeberis curvirostris",
                       "Bosmina hagmanni",
                       "Bosmina longirostris",
                       "Acroperus harpae",
                       "Alona affinis",
                       "Alona costata",
                       "Alona rectangula",
                       "Alonella exigua",
                       "Alonella nana",
                       "Chydorus sphaericus",
                       "Dunhevedia crassa",
                       "Eurycercus lamellatus",
                       "Eurycercus lamellatus",
                       "Graptoleberis testudinaria",
                       "Leydigia leydigi",
                       "Pleuroxus aduncus",
                       "Caenestheriella setosa",
                       "Ceriodaphnia quadrangula",
                       "Ceriodaphnia reticulata",
                       "Ceriodaphnia silvestrii",
                       "Daphnia ambigua",
                       "Daphnia galeata",
                       "Daphnia gessneri",
                       "Daphnia magna",
                       "Daphnia magna",
                       "Daphnia pulex",
                       "Daphnia pulex",
                       "Diaphanosoma spinulosum",
                       "Scapholeberis mucronata",
                       "Simocephalus vetulus",
                       "Simocephalus vetulus",
                       "Simocephalus vetulus",
                       "Moina macrocopa",
                       "Moina micrura",
                       "Moina mongolica",
                       "Polyphemus pediculus",
                       "Diaphanosoma brachyurum",
                       "Hyalella",
                       "Aegla",
                       "Samastacus spinifrons",
                       "Acanthocyclops robustus",
                       "Mesocyclops ogunnus",
                       "Argyrodiaptomus azevedoi",
                       "Argyrodiaptomus azevedoi",
                       "Diaptomus arcticus",
                       "Diaptomus arcticus",
                       "Diaptomus leptopus",
                       "Diaptomus leptopus",
                       "Diaptomus siciloides",
                       "Diaptomus siciloides",
                       "Notodiaptomus cearensis",
                       "Notodiaptomus cearensis",
                       "Notodiaptomus evaldus",
                       "Notodiaptomus evaldus",
                       "Notodiaptomus iheringi",
                       "Notodiaptomus iheringi",
                       "Mesocyclops edax",
                       "Hedriodiscus truquii",
                       "Baetis",
                       "Ephemerella",
                       "Ecdyonurus",
                       "Epeorus",
                       "Rhithrogena",
                       "Ephoron album",
                       "Microvelia macgregori",
                       "Microvelia macgregori",
                       "Leuctra",
                       "Isoperla",
                       "Perlodes",
                       "Chilina gibbosa"),
           life.stage = c(NA,
                          NA,
                          NA,
                          NA,
                          NA,
                          "adult", 
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "adult",
                          "larva",
                          "nymph",
                          "nymph",
                          "nymph",
                          "nymph",
                          "nymph",
                          "nymph",
                          "adult",
                          "adult",
                          "nymph",
                          "nymph",
                          "nymph",
                          "adult"))

# rows to remove because cannot be confirmed
rows.remove <- c(64)
axell.lw <- axell.lw[-rows.remove, ]

# join these data
all(axell.lw$Species == ls_df$Species)
axell.lw$LifeStage <- ls_df$life.stage

# subset the axell 2021 database
axell.def <- axell.lw[, c(1, 2, 3, 4, 6, 7, 8)]
str(axell.def)

# add TL_middle from the midpoint of the min_max values
def.length <- 
  axell.def %>%
  mutate(TL_middle = if_else( is.na(TL_middle), (TL_min + TL_max)/2, TL_middle ) ) %>%
  select(Database, Genus, Species, LifeStage, TL_middle)
head(def.length)

# rename the TL_middle variable
def.length <- 
  def.length %>%
  rename(length_mid_mm = TL_middle)

# make a taxon variable instead of a genus-species classification
def.length <- 
  def.length %>%
  mutate(Taxon = ifelse(is.na(Species), Genus, Species)) %>%
  select(Database, Taxon, LifeStage, length_mid_mm) %>%
  rename(db_taxon = Taxon, life_stage = LifeStage)


## Joren 2022 data

# read in the data
jdat <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/freshwater_trait_databases/Joren_2022_length_data/Joren_2020_length_database.xlsx")
head(jdat)

# make a database column
jdat$Database <- "Joren_2022"

# subset the columns to match with the axell data
jdat <- 
  jdat %>%
  select(Database, Taxon, Life_stage, Length_mm)

# replace fake NAs with real NAs
jdat[jdat == "NA"] <- NA

# match the names
names(jdat) <- names(def.length)
View(jdat)

# which taxa to remove and generate separate length database
jdat_supp <- 
  jdat %>%
  filter(db_taxon %in% all_of(c("Nematoda", "Oligochaeta", "Turbellaria")) | life_stage == "tadpole" )

# remove these rows from the jdat data
jdat <- 
  jdat %>%
  filter(!(db_taxon %in% all_of(c("Nematoda", "Oligochaeta", "Turbellaria")) | life_stage == "tadpole"))


# bind these databases by row
def.length <- rbind(def.length, jdat)

# add a length id column
def.length$length_id <- 1:nrow(def.length)

# reorder the columns
def.length <- 
  def.length %>%
  select(length_id, db_taxon, Database, life_stage, length_mid_mm)

# write this into a database
write_csv(x = def.length, 
          file = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/default_length_data.csv")

# write the supplementary length data into a database
# add a length id column
jdat_supp$length_id <- 1:nrow(jdat_supp)

# reorder the columns
jdat_supp<- 
  jdat_supp %>%
  select(length_id, db_taxon, Database, life_stage, length_mid_mm)

write_csv(x = jdat_supp, 
          file = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/default_length_data_supp.csv")

### END
