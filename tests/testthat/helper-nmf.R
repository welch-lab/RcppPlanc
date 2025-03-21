skip_on_winbuilder <- function () {
  skip_if(!interactive() && !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false"))) && tolower(Sys.info()[["sysname"]]) == "windows")
}
