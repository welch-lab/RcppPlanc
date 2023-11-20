.onLoad <- function(libname, pkgname) {
    .openblaspthreadoff()
}

.onUnload <- function(libname, pkgname) {
  .openblaspthreadon()
}