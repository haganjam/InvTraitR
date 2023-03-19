# check and load DB files from GH for the installed version
.onLoad <- function(libname, pkgname) {
    download_database(pkgname)
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
