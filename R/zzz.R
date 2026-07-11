#' Retrieve List of NIMBLE Distributions Defined in faaime
#'
#' Function to rertrieve a list of the NIMBLE distributions defined in the
#' faaime package.
#'
#' @return A named list with each element containing the distribution
#' specification required for the \code{\link[nimble]{registerDistributions}}
#' function.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{registerDistributions}}
#' @export
getfaaimeDistributions <- function() {
  return(list(
    # Poisson Binomial distribution
    dpoisbinom = list(
      BUGSdist = "dpoisbinom(probs)",
      types = c("value = integer(0)", "probs = double(1)"),
      discrete = TRUE,
      pqAvail = TRUE,
      range = c(0, Inf)
    )
  ))
}

# Perform package setup upon loading of package namespace
.onLoad <- function(libname, pkgname) {
  # The location where the RcppArmadillo headers are held
  armadilloheaderloc <- file.path(find.package("RcppArmadillo"), "include")
  # The location where the Rcpp headers are held
  rcppheaderloc <- file.path(find.package("Rcpp"), "include")
  # The location where the NIMBLE-related source code is held
  nimblesrcloc <- file.path(find.package(pkgname), "nimble_ext")
  if(!dir.exists(nimblesrcloc)) {
    nimblesrcloc <- file.path(find.package(pkgname), "inst", "nimble_ext")
  }
  # ---- 1.1. Compile the NIMBLE FFT C++ bindings ----
  # Setup a temporary compilation location
  nimble_fft_compileloc <- gsub("\\", "/", tempfile(pattern = "nimble_fft_", fileext = ".o"), fixed = TRUE)
  # Compile the nimble_fft C++ code
  system(paste("g++ -std=c++14 -I\"", R.home("include"), "\" -I\"", armadilloheaderloc, "\" -I\"", rcppheaderloc, "\" \"", file.path(nimblesrcloc, "nimble_fft.cpp"), "\" -c -o \"", nimble_fft_compileloc, "\"", sep = ""))
  nimble_fft <- nimble::nimbleExternalCall(
    function(pReal = double(1), pImag = double(1), iLength = integer(0), bInverse = logical(0)){},
    Cfun = "nimble_fft",
    returnType = void(),
    headerFile = file.path(nimblesrcloc, "nimble_fft.h"),
    oFile = nimble_fft_compileloc
  )
  # Assign the re-defined NIMBLE function
  assignInMyNamespace("nimble_fft", nimble_fft)
  # ---- 1.2. Register the NIMBLE distributions ----
  # De-register the NIMBLE distributions that are currently registered
  faaimeDistributions <- getfaaimeDistributions()
  sapply(X = names(faaimeDistributions), FUN = nimble::deregisterDistributions, warn = FALSE)
  # Register the NIMBLE distributions
  nimble::registerDistributions(faaimeDistributions, verbose = FALSE)
}

# Perform package tidying after unloading the package namespace
.onUnload <- function(libpath) {
  # De-register the NIMBLE distributions that have been registered
  sapply(X = names(getfaaimeDistributions()), FUN = nimble::deregisterDistributions, warn = FALSE)
}
