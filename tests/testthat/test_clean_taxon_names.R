make_valid_df <- function() {
    a <- c("lorem ipsum", "dolor")
    b <- c("alpha beta", "gamma")
    data.frame(a, b)
}

test_that("given an unsupported database, when clean_taxon_names, then error", {
    expect_error(clean_taxon_names(NULL, NULL, NULL, "unsup"), ".*backbone.*")
})

test_that("given a data arg that's not a df,
            when clean_taxon_names,
            then error", {
    expect_error(clean_taxon_names(NULL, NULL, NULL, "gbif"))
    expect_error(clean_taxon_names(c("a", "b"), NULL, NULL, "gbif"))
})

test_that("given a dataframe and a column not found within,
            when clean_taxon_names,
            then error", {
    a <- c(1, 2)
    b <- c(3, 4)
    df <- data.frame(a, b)

    expect_error(clean_taxon_names(df, "nope", "larva", "gbif"))
    expect_error(clean_taxon_names(df, 1, "larva", "gbif"))
})

test_that("given a dataframe without string values in taxon col,
            when clean_taxon_names,
            then error", {
    a <- c(1, 2)
    b <- c(3, 4)
    df <- data.frame(a, b)

    expect_error(
        clean_taxon_names(df, "a", "larva", "gbif")
    )
})

test_that("given an unsupported life stage,
            when clean_taxon_names,
            then error", {
    expect_error(
        clean_taxon_names(make_valid_df(), "a", "unsup", "gbif"),
        regexp = ".*not.*appropriate life-stage.*"
    )
})

# TODO: add happy path tests
