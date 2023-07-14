#' @title get_habitat_data()
#' @description Get habitat data from Abell et al.'s (2008) freshwater
#'  ecoregion map
#' @details This function takes a data.frame with latitude and longitude data
#'  and gets data on biogeographic realm, major habitat type and ecoregion
#'  based on Abell et al.'s (2008) freshwater ecoregion map.
#' @author James G. Hagan (james_hagan(at)outlook.com) and Ronald Bergmann
#' @param data - data.frame with a column containing latitude and longitude in
#'  decimal degrees
#' @param latitude_dd - character string with the column name containing the
#'  latitude in decimal degrees
#' @param longitude_dd - character string with the column name containing the
#'  longitude in decimal degrees
#' @return tibble of the input data with habitat data from Abell et al.s (2008)
#'  ecoregion map attached as additional columns function to get habitat data
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.string
get_habitat_data <- function(data, latitude_dd, longitude_dd) {
    # check if the data input is a data.frame or a tibble
    assert_that(
        is.data.frame(data) || dplyr::is.tbl(data),
        msg = paste(data, "is not a data.frame or tibble object")
    )

    # check if the latitude and longitude columns are present in the data object
    assert_that(
        is.string(latitude_dd) &&
            is.string(longitude_dd) &&
            all(c(latitude_dd, longitude_dd) %in% names(data)),
        msg = paste(
            latitude_dd,
            "or",
            longitude_dd,
            "are not strings and/or are not present in the data object"
        )
    )

    # check if the data object has more than 0 rows
    assert_that(
        nrow(data) > 0,
        msg = paste(data, "object has zero rows and therefore no data")
    )

    # check if all lat/lon rows are numeric or all NA
    data_lat <- data[[latitude_dd]]
    data_lon <- data[[longitude_dd]]
    assert_that(
        (is.numeric(data_lat) && is.numeric(data_lon)) ||
            (all(is.na(data_lat)) && all(is.na(data_lon))), # TODO: this is intended?
        msg = paste(data_lat, "or", data_lon, "are not numeric variables")
    )

    # make sure the decimal degree variables are within the ranges of
    # the variables
    assert_that(
        all(data_lat[!is.na(data_lat)] <= 90) &&
            all(data_lat[!is.na(data_lat)] >= -90) &&
            all(data_lon[!is.na(data_lon)] <= 180) &&
            all(data_lon[!is.na(data_lon)] >= -180),
        msg = paste(
            data_lat,
            "or",
            data_lon,
            "are too big/small to be valid decimal degrees"
        )
    )

    # load the freshwater habitat map
    if (!exists("fw_map")) {
        fw_map <- readRDS(file = get_db_file_path("freshwater_ecoregion_map.rds"))
    }

    # load the freshwater habitat map
    if (!exists("fw_meta")) {
        fw_meta <- readRDS(file = get_db_file_path("freshwater_ecoregion_metadata.rds"))
    }

    # make sure the latitude-longitude data are numeric variables
    data[[latitude_dd]] <- as.numeric(data[[latitude_dd]])
    data[[longitude_dd]] <- as.numeric(data[[longitude_dd]])

    # add a row id column
    data$row_id <- 1:nrow(data)

    # remove the NA values
    lat.lon <- data[!(is.na(data[[longitude_dd]]) | is.na(data[[latitude_dd]])), ]

    # if there are no non-NA values then we don't get habitat data
    if (nrow(lat.lon) > 0) {
        # convert this to a spatial points object
        sp.pts <- sp::SpatialPoints(
            lat.lon[, c(longitude_dd, latitude_dd)],
            proj4string = raster::crs(fw_map)
        )

        # check where pts overlap with freshwater ecoregion
        hab.dat <- sp::over(sp.pts, fw_map)

        # add row id to these data
        hab.dat$row_id <- lat.lon$row_id
        names(hab.dat) <- c("habitat_id", "area_km2", "row_id")

        # join the habitat data to the original data
        data <- dplyr::full_join(data, hab.dat, by = "row_id")
    } else {
        data[["habitat_id"]] <- NA
        data[["area_km2"]] <- NA
        data[["row_id"]] <- NA
    }

    # arrange by row_id
    data <- dplyr::arrange(data, row_id)

    # remove the row_id column
    data <- dplyr::select(data, -row_id, -area_km2)

    # join these data to the metadata
    data <- dplyr::left_join(data, fw_meta, by = "habitat_id")

    # convert to a tibble
    data <- dplyr::as_tibble(data)

    data
}
