
# set-up test scenarios

# 1. various types of missing data 

# - length range
# - latitude-longitude data etc.

# 2. only special names


# load the search functions for testing

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

### END

