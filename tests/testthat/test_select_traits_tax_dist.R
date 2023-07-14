make_test_input <- function() {
  data.frame(
    row = 1:7,
    taxon_name = c(
      "Gammarus", "Daphnia", "Triops granitica", "Triops",
      "Simocephalus vetulus", "Turbellaria", "Oligochaeta"
    ),
    Life_stage = c(
      "adult", "adult", "adult", "adult",
      "adult", "none", "none"
    ),
    lat = rep(50.55, 7),
    lon = rep(4.98, 7),
    length_mm = rnorm(7, 10, 2),
    clean_taxon_name = c(
      "Gammarus", "Daphnia", "Triops granitica", "Triops",
      "Simocephalus vetulus", "Turbellaria", "Oligochaeta"
    ),
    db = c(
      "gbif", "gbif", "gbif", "gbif",
      "gbif", "special", "special"
    ),
    scientificName = c(
      "Gammarus", "Daphnia", NA, "Triops",
      "Simocephalus vetulus", NA, NA
    ),
    taxonRank = c("Genus", "Genus", NA, "Genus", 
                  "Species", NA, NA),
    acceptedNameUsageID = c(
      "GBIF:2218440", "GBIF:2234785", NA,
      "GBIF:2235057", "GBIF:2234807", NA, NA
    ),
    taxon_order = c(
      "Amphipoda", "Diplostraca", NA, "Notostraca",
      "Diplostraca", NA, NA
    ),
    taxon_family = c(
      "Gammaridae", "Daphniidae", NA, "Triopsidae",
      "Daphniidae", NA, NA
    ),
    habitat_id = c(
      404, 404, 404, 404, 404, 404, 404
    ),
    realm = c(
      "Palearctic", "Palearctic", "Palearctic",
      "Palearctic", "Palearctic", "Palearctic",
      "Palearctic"
    ),
    major_habitat_type = c(
      "Temperate floodplain rivers and wetlands",
      "Temperate floodplain rivers and wetlands",
      "Temperate floodplain rivers and wetlands",
      "Temperate floodplain rivers and wetlands",
      "Temperate floodplain rivers and wetlands",
      "Temperate floodplain rivers and wetlands",
      "Temperate floodplain rivers and wetlands"
    ),
    ecoregion = c(
      "Central & Western Europe",
      "Central & Western Europe",
      "Central & Western Europe",
      "Central & Western Europe",
      "Central & Western Europe",
      "Central & Western Europe",
      "Central & Western Europe"
    )
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
    select_traits_tax_dist(data.frame(), "tax", 3, "nope"))
})

test_that("test if select_traits_tax_dist() the column
            names that are outputted are correct", {
  # when
  output <- select_traits_tax_dist(
    data = make_test_input(),
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    body_size = "length_mm"
  )

  # extract names from each element of the output list
  expect_true( all(names(output) == c(
      c("row", "taxon_name", "Life_stage", "lat", "lon", "length_mm", 
        "clean_taxon_name", "db", "scientificName", "taxonRank", 
        "acceptedNameUsageID", "taxon_order", "taxon_family", 
        "habitat_id", "realm", "major_habitat_type", "ecoregion", 
        "trait_out", "db_scientificName", "id", "tax_distance", 
        "body_size_range_match", "life_stage_match", "r2_match", "n", 
        "db_min_body_size_mm", "db_max_body_size_mm", "realm_match", 
        "major_habitat_type_match", "ecoregion_match", "recommend", 
        "explanation", "workflow2_choice")
    )) 
    )
  
})

test_that("test if select_traits_tax_dist() outputs entries that should
            have scientific names do have a non-missing
            scientificName column", {
  # when
  output <- select_traits_tax_dist(
    data = make_test_input(),
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    body_size = "length_mm"
  )
  
  # extract unique expected names
  x <- 
    output |>
    dplyr::group_by(row) |>
    dplyr::summarise(scientificName = unique(scientificName)) |>
    dplyr::pull(scientificName)
  
  # set the correct answers
  y <- c(
    "Gammarus",
    "Daphnia",
    NA,
    "Triops",
    "Simocephalus vetulus",
    NA,
    NA
  )
  
  # test if these are equal
  z <- mapply(function(x, y){
    (x == y) | (is.na(x) && is.na(y))
  }, x, y, SIMPLIFY = TRUE, USE.NAMES = FALSE)

  # make sure the outputted scientific names are correct
  expect_true(all(z))
  
})

test_that("test if select_traits_tax_dist() outputs
            the taxonomic distances properly", {
  # when
  output <- select_traits_tax_dist(
    data = make_test_input(),
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    body_size = "length_mm"
  )

  # all taxonomic distances should be numeric or NA
  expect_true( 
    all( is.numeric(output[["tax_distance"]]) | is.na(output[["tax_distance"]]) )
    )
  
})

test_that("test if select_traits_tax_dist() works
            correctly with only special names", {
              
  test_input <- make_test_input()   
  
  # run the select_traits_tax_dist() function with only special names
  output1 <- select_traits_tax_dist(
    data = test_input[c(6, 7), ],
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    body_size = "length_mm"
  )

  # run the select_traits_tax_dist() function with all names
  output2 <- select_traits_tax_dist(
    data = test_input,
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    body_size = "length_mm"
  )

  # make sure output is correct
  spec <- output1[, names(output1) != "row"]
  all <- output2[c(27:30),][, names(output2[c(27:30),]) != "row"]
  x <- (spec == all)
  
  expect_true(all( all(is.na(x) | (x == TRUE)) ))
  
})

test_that("test if select_traits_tax_dist() works
            correctly without any special names", {
              
  test_input <- make_test_input()             
  
  # run the select_traits_tax_dist() function without special names
  output1 <- select_traits_tax_dist(
    data = test_input[-c(6, 7), ],
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    body_size = "length_mm"
  )

  # run the select_traits_tax_dist() function with all names
  output2 <- select_traits_tax_dist(
    data = test_input,
    target_taxon = "taxon_name",
    life_stage = "Life_stage",
    body_size = "length_mm"
  )

  # make sure output is correct
  spec <- output1[, names(output1) != "row"]
  all <- output2[-c(27:30),][, names(output2[-c(27:30),]) != "row"]
  x <- (spec == all)
  
  expect_true(all( all(is.na(x) | (x == TRUE)) ))
  
})

