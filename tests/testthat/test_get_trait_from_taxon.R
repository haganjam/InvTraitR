make_test_input <- function() {
  df <-
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
  df[9, ]$body_size_mm <- NA

  df
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

test_that("if all taxon_names are not present in output,
           when get_trait_from_taxon,
           then error", {
  input <- make_test_input()

  # when
  output <- get_trait_from_taxon(
    data = input,
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    latitude_dd = "lat",
    longitude_dd = "lon",
    body_size = "body_size_mm",
    max_tax_dist = 3,
    trait = "equation",
    gen_sp_dist = 0.5
  )

  # test if all the relevant taxon names are present in the output
  expect_true(all(output[["taxon_name"]] == input[["taxon_name"]]))
})

test_that("given the equation, body_size_mm or id columns are NA,
           when get_trait_from_taxon,
           then dry_biomass_mg should be NA", {
  input <- make_test_input()

  # when
  output <- get_trait_from_taxon(
    data = input,
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    latitude_dd = "lat",
    longitude_dd = "lon",
    body_size = "body_size_mm",
    max_tax_dist = 3,
    trait = "equation",
    gen_sp_dist = 0.5
  )

  # if the equation and body_size_mm or id colums are NA,
  # then dry_biomass_mg should be NA

  # equation column
  expect_true(all(is.na(output[["equation"]]) == is.na(output[["dry_biomass_mg"]])))

  # body_size_mm
  x <- is.na(output[["body_size_mm"]])
  expect_true(all(is.na(output[["body_size_mm"]][x]) == is.na(output[["dry_biomass_mg"]][x])))

  # id
  y <- is.na(output[["id"]])
  expect_true(all(is.na(c(output[["dry_biomass_mg"]][y]))))
})

test_that("given some of the outputted taxonomic distance values
             are greater than the max_tax_distance argument,
           when get_trait_from_taxon,
           then error", {
  input <- make_test_input()

  # when
  output <- get_trait_from_taxon(
    data = input,
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    latitude_dd = "lat",
    longitude_dd = "lon",
    body_size = "body_size_mm",
    max_tax_dist = 3,
    trait = "equation",
    gen_sp_dist = 0.5
  )

  # test if the outputted taxonomic distances are less than the max tax distance
  expect_true(all(output[["tax_distance"]][!is.na(output[["tax_distance"]])] <= 3))
})

test_that("does the get_trait_from_taxon() function output the correct
          additional identifier columns?", {
  input <- make_test_input()

  # add additional columns
  input$sex <- "male"

  # when
  output <- get_trait_from_taxon(
    data = input,
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    latitude_dd = "lat",
    longitude_dd = "lon",
    body_size = "body_size_mm",
    max_tax_dist = 3,
    trait = "equation",
    gen_sp_dist = 0.5
  )

  # test if the sex column is in the output
  expect_true("sex" %in% names(output))
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
