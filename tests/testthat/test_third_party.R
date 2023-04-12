# we want to verify these once, so that we notice breaking changes
# in the dependencies, but then replace the actual calls with test
# doubles to speed things up

# -- tests

test_that("bdc::bdc_clean_names actually returns expected results", {
    expected <- data.frame(
        scientificName = c(
            "Gammarus_", "Daphnia", "Triops granitica",
            "Triops", "Simocephalus vetulus", NA,
            "Turbellaria", "Nematoda", "Bae.tidae"
        ),
        names_clean = c(
            "Gammarus", "Daphnia", "Triops granitica",
            "Triops", "Simocephalus vetulus", NA,
            "Turbellaria", "Nematoda", NA
        )
    )

    actual <- bdc::bdc_clean_names(
        sci_names = c(
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
        save_outputs = FALSE
    )
    expect_equal(actual[["scientificName"]], expected[["scientificName"]])
    expect_equal(actual[["names_clean"]], expected[["names_clean"]])
})

test_that("bdc::bdc_query_names_taxadb returns expected result", {
    expected <- data.frame(
        scientificName = c(
            "Gammarus", "Daphnia", NA,
            "Triops", "Simocephalus vetulus", NA, NA
        ),
        taxonRank = c("genus", "genus", NA, "genus", "species", NA, NA),
        acceptedNameUsageID = c(
            "GBIF:2218440", "GBIF:2234785", NA,
            "GBIF:2235057", "GBIF:2234807", NA, NA
        )
    )

    actual <- bdc::bdc_query_names_taxadb(
        sci_name = c(
            "Gammarus", "Daphnia", "Triops granitica",
            "Triops", "Simocephalus vetulus", NA, NA
        ),
        db = "gbif",
        rank_name = "Animalia",
        rank = "kingdom",
        export_accepted = FALSE
    )

    expect_equal(actual[["scientificName"]], expected[["scientificName"]])
    expect_equal(actual[["taxonRank"]], expected[["taxonRank"]])
    expect_equal(actual[["acceptedNameUsageID"]], expected[["acceptedNameUsageID"]])
})

test_that("taxadb::td_create works", {
    taxadb::td_create(provider = "gbif")
    expect_equal(TRUE, TRUE) # dummy
})
