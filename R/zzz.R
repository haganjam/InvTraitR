# check and load DB files from GH for the installed version
.onLoad <- function(libname, pkgname) {
    pkg_version <- utils::packageVersion(pkgname)
    local_db_dir <- file.path(
        fs::path_package(pkgname),
        .Platform$file.sep,
        file.path("database"),
        .Platform$file.sep
    )
    remote_db_dir <- "database"
    dir.create(local_db_dir, showWarnings = FALSE)
    base_url <- "https://github.com/haganjam/FW_invert_biomass_allometry"

    db_files <- list(
        "col_higher_taxon_matrices.rds",
        "col_taxon_database.rds",
        "equation_database.rds",
        "freshwater_ecoregion_data.rds",
        "freshwater_ecoregion_map.rds",
        "freshwater_ecoregion_metadata.rds",
        "gbif_higher_taxon_matrices.rds",
        "gbif_taxon_database.rds",
        "itis_higher_taxon_matrices.rds",
        "itis_taxon_database.rds",
        "taxon_database.rds"
    )

    versioned_base_url <- NULL
    for (file in db_files) {
        file_path <- file.path(local_db_dir, file)
        if (!file.exists(file_path)) {
            # check if versioned files are available
            # - pull mainline otherwise
            if (is.null(versioned_base_url)) {
                tag_available <- httr::GET(paste(
                    base_url,
                    "blob",
                    pkg_version,
                    remote_db_dir,
                    sep = "/"
                ))$status_code == 200

                if (tag_available) {
                    versioned_base_url <- paste(
                        base_url,
                        "blob",
                        pkg_version,
                        "database",
                        sep = "/"
                    )
                } else {
                    versioned_base_url <- paste(
                        base_url,
                        "raw/main/database",
                        sep = "/"
                    )
                }
            }

            packageStartupMessage(sprintf(
                "required DB %s not found locally; downloading ...",
                file_path
            ))
            download.file(paste(versioned_base_url, file, sep = "/"), file_path)
            packageStartupMessage("done")
        }
    }
    invisible()
}

.onAttach <- function(libname, pkgname) {
    # let's greet our users :)
    pkg_version <- utils::packageVersion(pkgname)
    packageStartupMessage(sprintf(
        "Welcome to %s in version %s",
        pkgname,
        pkg_version
    ))
    invisible()
}
