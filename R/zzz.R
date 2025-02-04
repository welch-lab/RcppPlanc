.onLoad <- function(libname, pkgname) {
  .openblaspthreadoff(getLoadedDLLs()[["RcppPlanc"]][[4]])
}

.onUnload <- function(libname, pkgname) {
  .openblaspthreadon(getLoadedDLLs()[["RcppPlanc"]][[4]])
}