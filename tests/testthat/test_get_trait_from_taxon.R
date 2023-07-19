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
  x <- make_test_input()

  # when
  y <- get_trait_from_taxon(
    data = x,
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
  expect_true(all(y[["data"]][["taxon_name"]] == x[["taxon_name"]]))
})

test_that("given the equation, body_size_mm or id columns are NA,
           when get_trait_from_taxon,
           then dry_biomass_mg should be NA", {
  x <- make_test_input()

  # when
  y <- get_trait_from_taxon(
    data = x,
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
  expect_true(all(is.na(y[["data"]][["equation"]]) == is.na(y[["data"]][["dry_biomass_mg"]])))

  # body_size_mm
  z <- is.na(y[["data"]][["body_size_mm"]])
  expect_true(all(is.na(y[["data"]][["body_size_mm"]][z]) == is.na(y[["data"]][["dry_biomass_mg"]][z])))

  # id
  u <- is.na(y[["data"]][["id"]])
  expect_true(all(is.na(c(y[["data"]][["dry_biomass_mg"]][u]))))
})

test_that("given some of the outputted taxonomic distance values
             are greater than the max_tax_distance argument,
           when get_trait_from_taxon,
           then error", {
  x <- make_test_input()

  # when
  y <- get_trait_from_taxon(
    data = x,
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
  expect_true(all(y[["data"]][["tax_distance"]][!is.na(y[["data"]][["tax_distance"]])] <= 3))
})

test_that("does the get_trait_from_taxon() function output the correct
          additional identifier columns?", {
  x <- make_test_input()

  # add additional columns
  x$sex <- "male"

  # when
  y <- get_trait_from_taxon(
    data = x,
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
  expect_true("sex" %in% names(y[["data"]]))
})

test_that("test a highly marginal case where
            there are no matches for the life-stages", {
              
    x <- 
      data.frame(
        taxon_name = c("Gammarus", "Daphnia"),
        Life_stage = c("larva", "none"),
        lat = rep(50.5, 1),
        lon = rep(4.98, 1),
        body_size_mm = rnorm(2, 10, 2)
      )
                
  y <- 
    get_trait_from_taxon(
      data = x,
      target_taxon = "taxon_name",
      life_stage = "Life_stage",
      latitude_dd = "lat",
      longitude_dd = "lon",
      body_size = "body_size_mm",
      max_tax_dist = 3,
      trait = "equation",
      gen_sp_dist = 0.5
    )
  
  expect_true(all(y[["decision_data"]][["workflow2_choice"]] == FALSE))
  
})

test_that("test the case where there are only special names", {
  
  x <- data.frame(
    taxon_name = c("Oligochaeta", "Oligochaeta", "Turbellaria"),
    Life_stage = c("none", "none", "none"),
    lat = rep(50.5, 1),
    lon = rep(4.98, 1),
    body_size_mm = rnorm(3, 10, 2)
  )

  y <- get_trait_from_taxon(
    data = x,
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    latitude_dd = "lat",
    longitude_dd = "lon",
    body_size = "body_size_mm",
    max_tax_dist = 3,
    trait = "equation",
    gen_sp_dist = 0.5
  )

  expect_equal(x[["taxon_name"]], y[["data"]][["taxon_name"]])
  
})

# test if the decision df matches the output df
test_that("does the decision data match the output data", {
            
  x <- make_test_input()
            
  y <- get_trait_from_taxon(
    data = x,
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    latitude_dd = "lat",
    longitude_dd = "lon",
    body_size = "body_size_mm",
    max_tax_dist = 3,
    trait = "equation",
    gen_sp_dist = 0.5
  )
  
  # get the samples for which equations were given
  u <- dplyr::filter(y[["data"]], !is.na(dry_biomass_mg))
  z <- dplyr::filter(y[["decision_data"]], workflow2_choice == TRUE)
  
  expect_true(all(unique(u$taxon_name) == unique(z$taxon_name)))
            
})

