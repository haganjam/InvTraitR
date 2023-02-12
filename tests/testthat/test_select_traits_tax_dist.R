make_test_input <- function() {
    data.frame(
        taxon_name = c(
            "Gammarus", "Daphnia", "Triops granitica", "Triops",
            "Simocephalus vetulus", "Turbellaria", "Oligochaeta"
        ),
        Life_stage = c(
            "adult", "adult", "adult", "adult",
            "adult", "none", "none"
        ),
        lat = rep(50.55, 7),
        lon = rep(4.98, 7)
    )
}

# extract_genus ---

test_that("given a non-string input, when extract_genus, then error out", {
    expect_error(extract_genus(1))
})

test_that("given non-letter inputs, when extract_genus, then error out", {
    expect_error(extract_genus("some sp3(!al input"))
    expect_error(extract_genus("some.interpunct"))
})

test_that(
    "given a single word, when extract_genus, then return it with attr n=1",
    {
        word <- "word"

        extracted <- extract_genus(word)

        expect_equal(extracted[1], word)
        expect_length(attributes(extracted), 1)
        expect_equal(1, attributes(extracted)$n)
    }
)

test_that(
    "given N words, when extract_genus, then return it with attr n=N",
    {
        word <- "lorem"
        words <- paste(word, "ipsum dolor amet")

        extracted <- extract_genus(words)

        expect_equal(extracted[1], word)
        expect_length(attributes(extracted), 1)
        expect_equal(4, attributes(extracted)$n)
    }
)

# select_traits_tax_dist ---

test_that("given not a number for dists,
             when select_traits_tax_dist,
             then error", {
    expect_error(select_traits_tax_dist(NULL, "", "not a number"))
    expect_error(select_traits_tax_dist(NULL, "", 1, "equation", "not a number"))
})

test_that("given a number < 0 for dists,
             when select_traits_tax_dist,
             then error", {
    expect_error(select_traits_tax_dist(NULL, "", 0))
    expect_error(select_traits_tax_dist(NULL, "", -1))
    expect_error(select_traits_tax_dist(NULL, "", 1, "equation", 0))
    expect_error(select_traits_tax_dist(NULL, "", 1, "equation", -1))
})

test_that("given an unsupported trait,
            when select_traits_tax_dist
            then error", {
    expect_error(
        select_traits_tax_dist(data.frame(), "tax", 3, "nope"),
        ".*trait.*"
    )
})


test_that("test if select_traits_tax_dist() the column
            names that are outputted are correct", {
    # TODO: this should be static and not rely on other units (functions)
    x1 <- clean_taxon_names(
        data = make_test_input(),
        target_taxon = "taxon_name", life_stage = "Life_stage",
        database = "gbif"
    )

    # TODO: this should be static and not rely on other units (functions)
    y1 <- get_habitat_data(data = x1, latitude_dd = "lat", longitude_dd = "lon")

    # when
    z1 <- select_traits_tax_dist(data = y1, target_taxon = "taxon_name")

    t1 <- lapply(z1, function(input) {
        all(names(input) == c(
            "taxon_name", "Life_stage", "lat", "lon", "clean_taxon_name",
            "db", "scientificName", "acceptedNameUsageID",
            "db_taxon_higher_rank", "db_taxon_higher", "habitat_id",
            "realm", "major_habitat_type", "ecoregion",
            "db.scientificName", "trait_out", "id", "tax_distance"
        ))
    })

    expect_true(all(unlist(t1)))
})

test_that("test if select_traits_tax_dist() outputs entries that should
            have scientific names do have a non-missing
            scientificName column", {
    # TODO: this should be static and not rely on other units (functions)
    x1 <- clean_taxon_names(
        data = make_test_input(),
        target_taxon = "taxon_name", life_stage = "Life_stage",
        database = "gbif"
    )

    # TODO: this should be static and not rely on other units (functions)
    y1 <- get_habitat_data(data = x1, latitude_dd = "lat", longitude_dd = "lon")

    # when
    z1 <- select_traits_tax_dist(data = y1, target_taxon = "taxon_name")

    t2 <- sapply(z1, function(input) unique(input[["scientificName"]])) == c(
        "Gammarus",
        "Daphnia",
        NA,
        "Triops",
        "Simocephalus vetulus",
        NA,
        NA
    )
    expect_true(all(is.na(t2) | (t2 == TRUE)))
})

test_that("test if select_traits_tax_dist() outputs
            the taxonomic distances properly", {
    # TODO: this should be static and not rely on other units (functions)
    x1 <- clean_taxon_names(
        data = make_test_input(),
        target_taxon = "taxon_name", life_stage = "Life_stage",
        database = "gbif"
    )

    # TODO: this should be static and not rely on other units (functions)
    y1 <- get_habitat_data(data = x1, latitude_dd = "lat", longitude_dd = "lon")

    # when
    z1 <- select_traits_tax_dist(data = y1, target_taxon = "taxon_name")

    t3 <- sapply(z1, function(input) {
        all(is.numeric(input[["tax_distance"]]) | is.na(input[["tax_distance"]]))
    })
    expect_true(all(t3))
})

test_that("test if select_traits_tax_dist() works
            correctly with only special names", {
    # get the special names from the make_test_input() data.frame
    df.test2 <- make_test_input()[c(6, 7), ]

    # TODO: this should be static and not rely on other units (functions)
    x2 <- clean_taxon_names(
        data = df.test2,
        target_taxon = "taxon_name", life_stage = "Life_stage",
        database = "gbif"
    )

    # TODO: this should be static and not rely on other units (functions)
    y2 <- get_habitat_data(data = x2, latitude_dd = "lat", longitude_dd = "lon")

    # run the select_traits_tax_dist() function
    z2 <- select_traits_tax_dist(data = y2, target_taxon = "taxon_name")

    t4 <- mapply(function(spec, all) {
        spec <- spec[, names(spec) != "row_id"]
        all <- all[, names(all) != "row_id"]
        x <- (spec == all)
        all(is.na(x) | (x == TRUE))
    }, z2, z1[c(6, 7)])
    expect_true(all(t4))
})

test_that("test if select_traits_tax_dist() works
            correctly without any special names", {
    # get the special names from the make_test_input() data.frame
    df.test3 <- make_test_input()[-c(6, 7), ]

    # TODO: this should be static and not rely on other units (functions)
    x3 <- clean_taxon_names(
        data = df.test3,
        target_taxon = "taxon_name", life_stage = "Life_stage",
        database = "gbif"
    )

    # TODO: this should be static and not rely on other units (functions)
    y3 <- get_habitat_data(data = x3, latitude_dd = "lat", longitude_dd = "lon")

    # run the select_traits_tax_dist() function
    z3 <- select_traits_tax_dist(data = y3, target_taxon = "taxon_name")

    t5 <- mapply(function(spec, all) {
        spec <- spec[, names(spec) != "row_id"]
        all <- all[, names(all) != "row_id"]
        x <- (spec == all)
        all(is.na(x) | (x == TRUE))
    }, z3, z1[-c(6, 7)])
    expect_true(all(t4))
})
