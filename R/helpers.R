appdata_subdir <- "InvTraitR"
db_dir <- "database"

#' @title get_db_file_path
#' @description Will try to fetch the path of given DB file from:
#'  * platform-specific appdata dir
#'  * project dir/database (during development)
#'  * testdata dir
#' @param filename DB file to look up
#' @return path to the requested file
get_db_file_path <- function(filename) {
    ## appdir
    appdir_path <- file.path(
        rappdirs::user_data_dir(),
        appdata_subdir,
        db_dir,
        filename
    )
    if (file.exists(appdir_path)) {
        return(appdir_path)
    }

    # project dir
    project_path <- file.path(
        db_dir,
        filename
    )
    if (file.exists(project_path)) {
        return(project_path)
    }

    # test dir
    test_path <- system.file(
        "testdata",
        filename
    )
    if (test_path != "") {
        return(test_path)
    }

    # fail
    stop(paste("DB file '", filename, "' could not be found"))
}

#' @title update_user_db
#' @description Copy over DB files to appdata dir. Comes in
#'  handy in case they got stale and you want to update them.
update_user_db <- function() {
    appdata_db_path <- file.path(
        rappdirs::user_data_dir(),
        appdata_subdir,
        db_dir
    )
    if (!dir.exists(appdata_db_path)) {
        dir.create(appdata_db_path, recursive = TRUE)
    }

    files <- list.files(db_dir)
    for (file in files) {
        file.copy(
            from = file.path(db_dir, file),
            to = file.path(
                appdata_db_path,
                file
            ),
            overwrite = TRUE
        )
    }
}

#' @title download_database
#' @description Will setup the local database in a directory
#'  chosen by rappdirs. Will download DB files for the current
#'  package version (tag) if it can be determined and is available
#'  or from mainline otherwise.
download_database <- function(pkgname) {
    pkg_version <- utils::packageVersion(pkgname)
    local_db_dir <- file.path(
        rappdirs::user_data_dir(),
        appdata_subdir,
        db_dir
    )
    remote_db_dir <- "database"

    base_url <- "https://github.com/haganjam/InvTraitR"
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
    local_db_dir_exists <- file.exists(local_db_dir)
    consent <- !interactive()
    for (file in db_files) {
        file_path <- file.path(local_db_dir, file)
        if (!file.exists(file_path)) {
            # ask for consent when run interactively
            if (!consent) {
                consent <- utils::menu(c("Yes", "No"), title = paste(
                    "Will create",
                    local_db_dir,
                    "and download DB files from GitHub - okay?"
                ))
                if (consent != "1") {
                    stop("DB setup aborted.")
                }
            }

            # create local db dir
            if (!local_db_dir_exists) {
                local_db_dir_exists <- dir.create(
                    local_db_dir,
                    recursive = TRUE,
                    showWarnings = FALSE
                )
            }

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
            utils::download.file(
                paste(versioned_base_url, file, sep = "/"),
                file_path
            )
            packageStartupMessage("done")
        }
    }
}
