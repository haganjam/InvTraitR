test_that("given a bad workflow,
            when get_trait_from_taxon,
            then error", {
    expect_error(
        get_trait_from_taxon(data.frame(), "tax", "larva", "lat", "lon", "size", "unsupported"),
        ".*workflow.*"
    )
})

# TODO: add tests
