
# set a seed to make it reproducible
set.seed(347854)

# load the search functions for testing
source(here::here("scripts/03_use_database/01_search_database_ver2.R"))

# Get_Habitat_Data tests

# set an error message
error_string <- "Get_Habitat_Data() did not correctly extract habitat information from test data"

# set-up the test data
df.test1 <- data.frame(latitude = c(49.76, # regular lat-lon point (Europe)
                                    47.19, # regular lat-lon point (North America)
                                    -33.56, # regular lat-lon point (Southern Africa)
                                    42.71, # lat-lon point in the ocean
                                    -21.55, # regular lat-lon point on an island (Madagascar)
                                    NA, # missing latitude data
                                    21.44, # missing longitude data
                                    NA ), # missing latitude and longitude data
                       longitude = c(3.12, 
                                     -105.79, 
                                     20.859, 
                                     -153.59, 
                                     44.78, 
                                     48.13, 
                                     NA, NA) )

# set-up the correct result
df.test1.out <- data.frame(habitat_id = c(404, 142, 578, NA, 579, NA, NA, NA),
                           major_habitat_type = c("Temperate floodplain rivers and wetlands", 
                                                  "Temperate upland rivers",
                                                  "Temperate coastal rivers",
                                                  NA,
                                                  "Xeric freshwaters and endorheic (closed) basins", 
                                                  NA, 
                                                  NA,
                                                  NA),
                           ecoregion = c("Central & Western Europe", 
                                         "Upper Missouri",
                                         "Cape Fold",
                                         NA,
                                         "Western Madagascar",
                                         NA, 
                                         NA, 
                                         NA)
)

# count how many NAs there should be in each row after running Get_Habitat_Data() on the df.test1 data.frame
test1.na.output <- c(0, 0, 0, 4, 0, 5, 5, 6)

# test1: Does the Get_Habitat_Data() function obtain the correct information?

# run the function on the test data
x <- Get_Habitat_Data(data = df.test1, latitude_dd = "latitude", longitude_dd = "longitude")

# test whether all derived entries are correct
y <- unlist(x[, names(df.test1.out)], use.names = FALSE) == unlist(df.test1.out, use.names = FALSE)

# all should either be true or NA
assertthat::assert_that(assertthat::see_if(all(y == TRUE | is.na(y))), 
                        msg = error_string)

assertthat::assert_that(assertthat::see_if( all( apply(x, 1, function(x) sum(is.na(x)) ) == test1.na.output ) ), 
                        msg = error_string)


# test2: Are the output columns correct?

# add an identifier column to the df.test1 data
df.test2 <- dplyr::mutate(df.test1, site = 1:nrow(df.test1))

# run the function
x <- Get_Habitat_Data(data = df.test2, latitude_dd = "latitude", longitude_dd = "longitude")

# test if the columns are there and whether they are correct
assertthat::assert_that(assertthat::see_if( all( names(x) == c("latitude", "longitude", "site", "habitat_id", "realm", "major_habitat_type", "ecoregion")) ), 
                        msg = error_string)

# test if the identifier column is correctly attached
assertthat::assert_that(assertthat::see_if(all( x[["site"]] == df.test2[["site"]] ) ),
                        msg = error_string)


# Clean_Taxon_Names() tests

# set an error message
error_string <- "Clean_Taxon_Names() does not properly clean the given taxon names"

# set-up the test data
df.test1 <- data.frame(taxon_name = c("Gammarus_", 
                                      "Daphnia", 
                                      "Triops granitica",
                                      "Triops",
                                      "Simocephalus vetulus",
                                      NA,
                                      "Turbellaria",
                                      "Nematoda", 
                                      "Bae.tidae"),
                       Life_stage = c("adult", 
                                      "adult", 
                                      "adult", 
                                      "adult",
                                      "adult",
                                      NA,
                                      "none",
                                      "none", 
                                      NA))

# set-up the correct result
clean_taxon_name <- c("Gammarus", "Daphnia", "Triops granitica", "Triops",
                      "Simocephalus vetulus", NA, NA, "Turbellaria", "Nematoda")
db <- c("gbif", "gbif", "gbif", "gbif", "gbif", NA, NA , "special", "special")
acceptedNameUsageID <- c("GBIF:2218440", "GBIF:2234785", NA, "GBIF:2235057",
                         "GBIF:2234807", NA, NA, NA, NA)
db_taxon_higher <- c("Amphipoda", "Diplostraca", NA, "Notostraca", 
                     "Diplostraca", NA, NA, NA, NA)

# test1: Does the Clean_Taxon_names() function obtain the correct information?

# run the function on the test data
x <- Clean_Taxon_Names(data = df.test1, 
                       target_taxon = "taxon_name", 
                       life_stage = "Life_stage", database = "gbif")

# test whether all derived entries are correct i.e. TRUE
t1 <- x$clean_taxon_name == clean_taxon_name
t2 <- x$db == db
t3 <- x$acceptedNameUsageID == acceptedNameUsageID
t4 <- x$db_taxon_higher == db_taxon_higher

# combine these four tests
y <- c(t1, t2, t3, t4)

# all should either be true or NA
assertthat::assert_that(assertthat::see_if(all(y == TRUE | is.na(y))), 
                        msg = error_string)

# test2: Does the Clean_Taxon_Names() function output the correct additional identifier columns?

# add an identifier column to the df.test1 data
df.test2 <- dplyr::mutate(df.test1, 
                          site = 1:nrow(df.test1),
                          sex = c("male", "female", "female", "male", "female", "male", "male", "female", "male"))
head(df.test2)

# run the function
x <- Clean_Taxon_Names(data = df.test2, 
                       target_taxon = "taxon_name", 
                       life_stage = "Life_stage", database = "gbif")

# test if the columns are there and whether they are correct
assertthat::assert_that(assertthat::see_if( all( names(x) == c("taxon_name", "life_stage", "site", "sex", "clean_taxon_name", "db", "scientificName", "acceptedNameUsageID", "db_taxon_higher_rank", "db_taxon_higher")) ), 
                        msg = error_string)

# test if the identifier columns are correctly attached
assertthat::assert_that(assertthat::see_if(all( x[["site"]] == c(1, 2, 3, 4, 5, 6, 9, 7, 8)  ) ),
                        msg = error_string)

assertthat::assert_that(assertthat::see_if(all( x[["sex"]] == c("male", "female", "female", "male", "female", "male", "male", "male", "female")  ) ),
                        msg = error_string)

# test3: Does the Clean_Taxon_Names() function work when there are only special names?
df.test3 <- df.test1[c(7, 8), ]
head(df.test3)

# run the function
x <- Clean_Taxon_Names(data = df.test3, 
                       target_taxon = "taxon_name", 
                       life_stage = "Life_stage", database = "gbif")

# test if the output is correct
assertthat::assert_that(assertthat::see_if( all(x$db == "special") ), 
                        msg = error_string)

assertthat::assert_that(assertthat::see_if( all(x$scientificName == c("Turbellaria", "Nematoda")) ), 
                        msg = error_string)

# test4: Does the Clean_Taxon_Names() function work when there are no special names?
df.test4 <- df.test1[-c(7, 8), ]
head(df.test4)

# run the function
x <- Clean_Taxon_Names(data = df.test4, 
                       target_taxon = "taxon_name", 
                       life_stage = "Life_stage", database = "gbif")

# test if the output is correct
assertthat::assert_that(assertthat::see_if( all(x$db == "gbif" | is.na(x$db)) ), 
                        msg = error_string)


# Select_Traits_Tax_Dist() tests

# set an error message
error_string <- "Select_Traits_Tax_Dist() does not correctly find traits or equations"

# generate some test data and run through the two functions
df.test1 <- 
  data.frame(taxon_name = c("Gammarus", "Daphnia", "Triops granitica", "Triops", 
                            "Simocephalus vetulus", "Turbellaria", "Oligochaeta"),
             Life_stage = c("adult", "adult", "adult", "adult",
                            "adult", "none", "none"),
             lat = rep(50.55, 7 ),
             lon = rep(4.98, 7))

# run the Clean_Taxon_Names() function
x1 <- Clean_Taxon_Names(data = df.test1, 
                        target_taxon = "taxon_name", life_stage = "Life_stage",
                        database = "gbif"
)

# run the Get_Habitat_Data() function
y1 <- Get_Habitat_Data(data = x1, latitude_dd = "lat", longitude_dd = "lon")

# run the Select_Traits_Tax_Dist() function
z1 <- Select_Traits_Tax_Dist(data = y1, target_taxon = "taxon_name")
View(z1[[7]])

# test1: test if Select_Traits_Tax_Dist() the column names that are outputted are correct

t1 <- 
  lapply(z1, function(input) {
    
    all(names(input) == c("taxon_name", "Life_stage", "lat", "lon", "clean_taxon_name",
                          "db", "scientificName", "acceptedNameUsageID", "db_taxon_higher_rank",
                          "db_taxon_higher", "habitat_id", "realm", "major_habitat_type", "ecoregion",
                          "row_id", "db.scientificName", "trait_out", "id", "tax_distance"
    ))
    
  })

assertthat::assert_that(assertthat::see_if( all(unlist(t1)) ), 
                        msg = error_string)

# test2: test if Select_Traits_Tax_Dist() outputs entries that should have scientific names do have a non-missing scientificName column

t2 <- sapply(z1, function(input) unique(input[["scientificName"]])) == c("Gammarus", "Daphnia", NA, "Triops", "Simocephalus vetulus", NA, NA)

assertthat::assert_that(assertthat::see_if( all(is.na(t2) | (t2 == TRUE)) ), 
                        msg = error_string)

# test3: test if Select_Traits_Tax_Dist() outputs the taxonomic distances properly

t3 <- sapply(z1, function(input) { 
  
  all(is.numeric(input[["tax_distance"]]) | is.na(input[["tax_distance"]]))
  
}  )

assertthat::assert_that(assertthat::see_if( all( t3 ) ), 
                        msg = error_string)

# test4: test if Select_Traits_Tax_Dist() works correctly with only special names

# get the special names from the df.test1 data.frame
df.test2 <- df.test[c(6, 7), ]

# run the Clean_Taxon_Names() function
x2 <- Clean_Taxon_Names(data = df.test2, 
                        target_taxon = "taxon_name", life_stage = "Life_stage",
                        database = "gbif"
)

# run the Get_Habitat_Data() function
y2 <- Get_Habitat_Data(data = x2, latitude_dd = "lat", longitude_dd = "lon")

# run the Select_Traits_Tax_Dist() function
z2 <- Select_Traits_Tax_Dist(data = y2, target_taxon = "taxon_name")

t4 <- 
  mapply(function(spec, all) {
    
    spec <- spec[, names(spec) != "row_id"]
    all <- all[, names(all) != "row_id"]
    x <- (spec == all)
    all(is.na(x) | (x == TRUE))
    
  }, z2, z1[c(6, 7)])

assertthat::assert_that(assertthat::see_if( all( t4 ) ), 
                        msg = error_string)

# test5: test if Select_Traits_Tax_Dist() works correctly without any special names

# get the special names from the df.test1 data.frame
df.test3 <- df.test[-c(6, 7), ]

# run the Clean_Taxon_Names() function
x3 <- Clean_Taxon_Names(data = df.test3, 
                        target_taxon = "taxon_name", life_stage = "Life_stage",
                        database = "gbif"
)

# run the Get_Habitat_Data() function
y3 <- Get_Habitat_Data(data = x3, latitude_dd = "lat", longitude_dd = "lon")

# run the Select_Traits_Tax_Dist() function
z3 <- Select_Traits_Tax_Dist(data = y3, target_taxon = "taxon_name")

t5 <- 
  mapply(function(spec, all) {
    
    spec <- spec[, names(spec) != "row_id"]
    all <- all[, names(all) != "row_id"]
    x <- (spec == all)
    all(is.na(x) | (x == TRUE))
    
  }, z3, z1[-c(6, 7)])

assertthat::assert_that(assertthat::see_if( all( t4 ) ), 
                        msg = error_string)



