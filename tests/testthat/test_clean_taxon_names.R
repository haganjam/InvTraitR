test_that("given an unsupported database, when clean_taxon_names, then error", {
    expect_error(clean_taxon_names(NULL, NULL, NULL, "unsup"))
})
