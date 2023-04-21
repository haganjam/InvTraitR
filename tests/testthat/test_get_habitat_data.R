test_that("given not a dataframe,
            when get_habitat_data,
            then error", {
    expect_error(
        get_habitat_data(c(), 1, 2),
        ".*data\\.frame.*"
    )
})

test_that("given a dataframe and mismatching lat/lon col names,
            when get_habitat_data,
            then error", {
    x <- c()
    lat <- c()
    lon <- c()
    df1 <- data.frame(x, lat)
    df2 <- data.frame(x, lon)

    expect_error(
        get_habitat_data(df1, "lat", "lon"),
        ".*not.*present.*"
    )
    expect_error(
        get_habitat_data(df2, "lat", "lon"),
        ".*not.*present.*"
    )
})

test_that("given an empty dataframe,
            when get_habitat_data,
            then error", {
    x <- c()
    lat <- c()
    lon <- c()
    df <- data.frame(x, lat, lon)

    expect_error(get_habitat_data(df, "lat", "lon"))
})

test_that("given an dataframe with bad lat/lon rows,
            when get_habitat_data,
            then error", {
    x <- c("a")
    lat <- c("nan")
    lon <- c(1)
    df <- data.frame(x, lat, lon)
    expect_error(get_habitat_data(df, "lat", "lon"), ".*not.*numeric.*")

    x <- c("a")
    lat <- c(1)
    lon <- c("nan")
    df <- data.frame(x, lat, lon)
    expect_error(get_habitat_data(df, "lat", "lon"), ".*not.*numeric.*")

    x <- c("a")
    lat <- c(1)
    lon <- c(NA)
    df <- data.frame(x, lat, lon)
    expect_error(get_habitat_data(df, "lat", "lon"), ".*not.*numeric.*")
})

test_that("given an dataframe with bad lat/lon rows,
            when get_habitat_data,
            then error", {
    x <- c("a")
    lat <- c(-999)
    lon <- c(1)
    df <- data.frame(x, lat, lon)
    expect_error(get_habitat_data(df, "lat", "lon"), ".*degree.*")

    x <- c("a")
    lat <- c(1)
    lon <- c(999)
    df <- data.frame(x, lat, lon)
    expect_error(get_habitat_data(df, "lat", "lon"), ".*degree.*")
})

setup_df <- function() {
    # set-up the test data
    df_test1 <- data.frame(
        latitude = c(
            49.76, # regular lat-lon point (Europe)
            47.19, # regular lat-lon point (North America)
            -33.56, # regular lat-lon point (Southern Africa)
            42.71, # lat-lon point in the ocean
            -21.55, # regular lat-lon point on an island (Madagascar)
            NA, # missing latitude data
            21.44, # missing longitude data
            NA
        ), # missing latitude and longitude data
        longitude = c(
            3.12,
            -105.79,
            20.859,
            -153.59,
            44.78,
            48.13,
            NA, NA
        )
    )

    # set-up the correct result
    df_test1_out <- data.frame(
        habitat_id = c(404, 142, 578, NA, 579, NA, NA, NA),
        major_habitat_type = c(
            "Temperate floodplain rivers and wetlands",
            "Temperate upland rivers",
            "Temperate coastal rivers",
            NA,
            "Xeric freshwaters and endorheic (closed) basins",
            NA,
            NA,
            NA
        ),
        ecoregion = c(
            "Central & Western Europe",
            "Upper Missouri",
            "Cape Fold",
            NA,
            "Western Madagascar",
            NA,
            NA,
            NA
        )
    )

    list(df_test1, df_test1_out)
}

test_that("count how many NAs there should be in each row after running
        get_habitat_data() on the df.test1 data.frame", {
    fixture <- setup_df()
    df_test1 <- fixture[[1]]
    df_test1_out <- fixture[[2]]

    # count how many NAs there should be in each row after running
    # get_habitat_data() on the df.test1 data.frame
    test1_na_output <- c(0, 0, 0, 4, 0, 5, 5, 6)

    # test1: Does the get_habitat_data() function
    # obtain the correct information?

    # run the function on the test data
    x <- get_habitat_data(
        data = df_test1,
        latitude_dd = "latitude",
        longitude_dd = "longitude"
    )

    # test whether all derived entries are correct
    y <- unlist(
        x[, names(df_test1_out)],
        use.names = FALSE
    ) == unlist(df_test1_out, use.names = FALSE)

    # all should either be true or NA
    expect_true(all(y == TRUE | is.na(y)))
    expect_true(all(apply(x, 1, function(x) sum(is.na(x))) == test1_na_output))
})

test_that("Are the output columns correct?", {
    fixture <- setup_df()
    df_test1 <- fixture[[1]]
    df_test1_out <- fixture[[2]]

    # add an identifier column to the df.test1 data
    df_test2 <- dplyr::mutate(df_test1, site = 1:nrow(df_test1))

    # run the function
    x <- get_habitat_data(
        data = df_test2,
        latitude_dd = "latitude",
        longitude_dd = "longitude"
    )

    # test if the columns are there and whether they are correct
    expect_true(all(names(x) == c(
        "latitude",
        "longitude",
        "site",
        "habitat_id",
        "realm",
        "major_habitat_type",
        "ecoregion"
    )), label = "unexpected column names")

    # test if the identifier column is correctly attached
    expect_true(all(x[["site"]] == df_test2[["site"]]))
})

