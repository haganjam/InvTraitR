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

# TODO: add happy path tests
