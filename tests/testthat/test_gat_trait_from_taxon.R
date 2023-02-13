make_test_input <- function() {
    df_test1 <-
        data.frame(
            taxon_name = c(
                "Gammarus", "Gammarus", "Gammarus", "Daphnia",
                "Triops granitica", "Triops",
                "Simocephalus vetulus", "Simocephalus vetulus",
                "Turbellaria", "Oligochaeta", "Oligochaeta", "Daphnia gessneri"
            ),
            Life_stage = c(
                "adult", "adult", "adult", "adult", "adult", "adult",
                "adult", "adult", "none", "none", "none",
                "adult"
            ),
            lat = c(rep(50.5, 6), rep(47.5, 6)),
            lon = c(rep(4.98, 6), rep(-105.4, 6)),
            body_size_mm = rnorm(12, 10, 2)
        )

    # add an NA to see how the functions react
    df_test1[9, ]$body_size_mm <- NA

    df_test1
}

test_that("given a bad workflow,
            when get_trait_from_taxon,
            then error", {
    expect_error(
        get_trait_from_taxon(
            data.frame(),
            "tax", "larva", "lat", "lon", "size", "unsupported"
        ),
        ".*workflow.*"
    )
})

test_that("a lot of things", {
    # TODO: break this apart

    x <- get_trait_from_taxon(
        data = df.test1,
        target_taxon = "taxon_name",
        life_stage = "Life_stage",
        latitude_dd = "lat",
        longitude_dd = "lon",
        body_size = "body_size_mm",
        max_tax_dist = 3,
        trait = "equation",
        gen_sp_dist = 0.5
    )

    # test if all the relevant taxon names are present
    expect_true(all(x[["taxon_name"]] == df.test1[["taxon_name"]]))

    # body size must be within the min and max
    t1 <- (x[["body_size_mm"]] >= x[["min_body_size_mm"]]) & (x[["body_size_mm"]] <= x[["max_body_size_mm"]])
    expect_true(all(t1 | is.na(t1)))

    # if the equation is NA then the dry_biomass_mg should be NA
    expect_true(all(is.na(x[["equation"]]) == is.na(x[["dry_biomass_mg"]])))

    # if the body_length_mm column is NA then the dry_biomass_mg column should be NA
    t2 <- is.na(x[["body_size_mm"]])
    expect_true(all(is.na(x[["body_size_mm"]][t2]) == is.na(x[["dry_biomass_mg"]][t2])))

    # if the id column is NA, then the dry_biomass_mg
    # and equation columns should be NA
    t3 <- is.na(x[["id"]])
    expect_true(all(is.na(c(x[["dry_biomass_mg"]][t3], x[["equation"]][t3]))))

    # check if the outputted taxonomic distances are less than the max tax distance
    expect_true(all(x[["tax_distance"]][!is.na(x[["tax_distance"]])] < 3))
})

test_that("test a highly marginal case where
            there are no matches for the life-stages", {
    x <- get_trait_from_taxon(
        data = data.frame(
            taxon_name = c("Gammarus", "Daphnia"),
            Life_stage = c("larva", "none"),
            lat = rep(50.5, 1),
            lon = rep(4.98, 1),
            body_size_mm = rnorm(2, 10, 2)
        ),
        target_taxon = "taxon_name",
        life_stage = "Life_stage",
        latitude_dd = "lat",
        longitude_dd = "lon",
        body_size = "body_size_mm",
        max_tax_dist = 3,
        trait = "equation",
        gen_sp_dist = 0.5
    )

    expect_true(all(is.na(x[, c("id", "tax_distance", names(x)[grepl("_match", x = names(x))])])))
})

test_that("test the case where there are only special names", {
    input_data <- data.frame(
        taxon_name = c("Oligochaeta", "Oligochaeta", "Turbellaria"),
        Life_stage = c("none", "none", "none"),
        lat = rep(50.5, 1),
        lon = rep(4.98, 1),
        body_size_mm = rnorm(3, 10, 2)
    )

    x <- get_trait_from_taxon(
        data = input_data,
        target_taxon = "taxon_name",
        life_stage = "Life_stage",
        latitude_dd = "lat",
        longitude_dd = "lon",
        body_size = "body_size_mm",
        max_tax_dist = 3,
        trait = "equation",
        gen_sp_dist = 0.5
    )

    expect_equal(input_data[["taxon_name"]], x[["taxon_name"]])
})
