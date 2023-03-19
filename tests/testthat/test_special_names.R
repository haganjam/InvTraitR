# This test is a bit stupid but may become more useful
# later in case the function evolves, and following TDD
# it's been here before its implementation (right?)
test_that("special_taxon_names yields correct names", {
    expect_equal(
        c(
            "Rotifera",
            "Tardigrada",
            "Nematoda",
            "Platyhelminthes",
            "Turbellaria",
            "Annelida",
            "Oligochaeta"
        ),
        special_taxon_names()
    )
})
