make_valid_df <- function() {
    a <- c("lorem ipsum", "dolor")
    b <- c("alpha beta", "gamma")
    data.frame(a, b)
}

# set-up the test data
make_test_input <- function() {
    data.frame(
        taxon_name = c(
            "Gammarus_",
            "Daphnia",
            "Triops granitica",
            "Triops",
            "Simocephalus vetulus",
            NA,
            "Turbellaria",
            "Nematoda",
            "Bae.tidae"
        ),
        Life_stage = c(
            "adult",
            "adult",
            "adult",
            "adult",
            "adult",
            NA,
            "none",
            "none",
            NA
        )
    )
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

test_that("given a dataframe without string values in target_taxon col,
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
        clean_taxon_names(make_valid_df(), "a", "b", "gbif"),
        regexp = ".*not.*appropriate life-stage.*"
    )
})

test_that("Does the clean_taxon_names() function
            obtain the correct information?", {
    # set-up the correct result
    clean_taxon_name <- c(
        "Gammarus", "Daphnia", "Triops granitica", "Triops",
        "Simocephalus vetulus", NA, NA, "Turbellaria", "Nematoda"
    )
    db <- c(
        "gbif", "gbif", "gbif", "gbif", "gbif",
        "gbif", "gbif", "special", "special"
    )
    accepted_name_usage_id <- c(
        "GBIF:2218440", "GBIF:2234785", NA, "GBIF:2235057",
        "GBIF:2234807", NA, NA, NA, NA
    )
    taxon_order = c(
      "Amphipoda", "Diplostraca", NA, "Notostraca",
      "Diplostraca", NA, NA, NA, NA
    )
    taxon_family = c(
      "Gammaridae", "Daphniidae", NA, "Triopsidae",
      "Daphniidae", NA, NA, NA, NA
    )

    # when
    x <- clean_taxon_names(
        data = make_test_input(),
        target_taxon = "taxon_name",
        life_stage = "Life_stage",
        database = "gbif"
    )

    # then
    expect_equal(clean_taxon_name, x$clean_taxon_name)
    expect_equal(db, x$db)
    expect_equal(accepted_name_usage_id, x$acceptedNameUsageID)
    expect_equal(taxon_order, x$taxon_order)
    expect_equal(taxon_family, x$taxon_family)
})

test_that("Does the clean_taxon_names() function output
             the correct additional identifier columns?", {
    df_test1 <- make_test_input()
    df_test2 <- dplyr::mutate(df_test1,
        site = 1:nrow(df_test1),
        sex = c(
            "male", "female", "female",
            "male", "female", "male",
            "male", "female", "male"
        )
    )

    # when
    x <- clean_taxon_names(
        data = df_test2,
        target_taxon = "taxon_name",
        life_stage = "Life_stage",
        database = "gbif"
    )

    # test if the columns are there and whether they are correct
    expect_equal(names(x), c(
        "taxon_name",
        "Life_stage",
        "site",
        "sex",
        "clean_taxon_name",
        "db",
        "scientificName",
        "taxonRank",
        "acceptedNameUsageID",
        "taxon_order",
        "taxon_family"
    ))

    # test if the identifier columns are correctly attached
    expect_true(all(x[["site"]] == c(1, 2, 3, 4, 5, 6, 9, 7, 8)))

    expect_true(all(x[["sex"]] == c(
        "male",
        "female",
        "female",
        "male",
        "female",
        "male",
        "male",
        "male",
        "female"
    )))
})

test_that("Does the clean_taxon_names() function work
            when there are only special names?", {
    input_data <- make_test_input()[c(7, 8), ]

    # when
    x <- clean_taxon_names(
        data = input_data,
        target_taxon = "taxon_name",
        life_stage = "Life_stage",
        database = "gbif"
    )

    # test if the output is correct
    expect_true(all(x$db == "special"))
    expect_true(all(x$clean_taxon_name == c("Turbellaria", "Nematoda")))
})

test_that("Does the clean_taxon_names() function work
             when there are no special names?", {
    input_data <- make_test_input()[-c(7, 8), ]

    # when
    x <- clean_taxon_names(
        data = input_data,
        target_taxon = "taxon_name",
        life_stage = "Life_stage",
        database = "gbif"
    )

    # test if the output is correct
    expect_true(all(x$db == "gbif" | is.na(x$db)))
})
