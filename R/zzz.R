.onLoad <- function(libname, pkgname) {
  .openblaspthreadoff(getLoadedDLLs()[["RcppPlanc"]][[4]])
  return(invisible(NULL))
}

.onUnload <- function(libpath) {
  .openblaspthreadon(getLoadedDLLs()[["RcppPlanc"]][[4]])
  library.dynam.unload("RcppPlanc", libpath)
  return(invisible(NULL))
}
