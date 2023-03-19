# We cannot really test our get_db_file_path() function
# in all its branches but let's test what's easy.

test_that("reflects filename in error message", {
    expect_error(
        get_db_file_path("not there"),
        ".*not there.*"
    )
})

test_that("finds some file", {
    f <- "equation_database.rds"
    expect_gte(length(get_db_file_path(f)), length(f))
})
