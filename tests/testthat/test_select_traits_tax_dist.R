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
    lon = rep(4.98, 7),
    clean_taxon = c(
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
    db_taxon_higher_rank = c(
      "order", "order", NA, "order", "order",
      NA, NA
    ),
    db_taxon_higher = c(
      "Amphipoda", "Diplostraca", NA, "Notostraca",
      "Diplostraca", NA, NA
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
    select_traits_tax_dist(data.frame(), "tax", 3, "nope"),
    ".*trait.*"
  )
})

test_that("test if select_traits_tax_dist() the column
            names that are outputted are correct", {
  # when
  output <- select_traits_tax_dist(
    data = make_test_input(),
    target_taxon = "taxon_name"
  )

  # extract names from each element of the output list
  x <- lapply(output, function(input) {
    all(names(input) == c(
      "taxon_name", "Life_stage", "lat", "lon", "clean_taxon",
      "db", "scientificName", "taxonRank", "acceptedNameUsageID",
      "db_taxon_higher_rank", "db_taxon_higher", "habitat_id",
      "realm", "major_habitat_type", "ecoregion",
      "db.scientificName", "trait_out", "id", "tax_distance"
    ))
  })

  # test if names in all the different list elements were correct
  expect_true(all(unlist(x)))
})

test_that("test if select_traits_tax_dist() outputs entries that should
            have scientific names do have a non-missing
            scientificName column", {
  # when
  output <- select_traits_tax_dist(
    data = make_test_input(),
    target_taxon = "taxon_name"
  )

  # make sure the outputted scientific names are correct
  x <- sapply(output, function(input) unique(input[["scientificName"]])) == c(
    "Gammarus",
    "Daphnia",
    NA,
    "Triops",
    "Simocephalus vetulus",
    NA,
    NA
  )

  # all should be true or NA
  expect_true(all(is.na(x) | (x == TRUE)))
})

test_that("test if select_traits_tax_dist() outputs
            the taxonomic distances properly", {
  # when
  output <- select_traits_tax_dist(
    data = make_test_input(),
    target_taxon = "taxon_name"
  )

  # get the taxonomic distance values for all names
  x <- sapply(output, function(input) {
    all(is.numeric(input[["tax_distance"]]) | is.na(input[["tax_distance"]]))
  })

  # all taxonomic distances should be numeric or NA
  expect_true(all(x))
})

test_that("test if select_traits_tax_dist() works
            correctly with only special names", {
  # run the select_traits_tax_dist() function with only special names
  output1 <- select_traits_tax_dist(
    data = make_test_input()[c(6, 7), ],
    target_taxon = "taxon_name"
  )

  # run the select_traits_tax_dist() function with all names
  output2 <- select_traits_tax_dist(
    data = make_test_input(),
    target_taxon = "taxon_name"
  )

  # make sure output is correct
  x <- mapply(function(spec, all) {
    spec <- spec[, names(spec) != "row_id"]
    all <- all[, names(all) != "row_id"]
    x <- (spec == all)
    all(is.na(x) | (x == TRUE))
  }, output1, output2[c(6, 7)])

  expect_true(all(x))
})

test_that("test if select_traits_tax_dist() works
            correctly without any special names", {
  # run the select_traits_tax_dist() function without special names
  output1 <- select_traits_tax_dist(
    data = make_test_input()[-c(6, 7), ],
    target_taxon = "taxon_name"
  )

  # run the select_traits_tax_dist() function with all names
  output2 <- select_traits_tax_dist(
    data = make_test_input(),
    target_taxon = "taxon_name"
  )

  x <- mapply(
    function(spec, all) {
      spec <- spec[, names(spec) != "row_id"]
      all <- all[, names(all) != "row_id"]
      x <- (spec == all)
      all(is.na(x) | (x == TRUE))
    },
    output1, output2[-c(6, 7)]
  )

  expect_true(all(x))
})
