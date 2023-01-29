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

# TODO: add happy path tests
