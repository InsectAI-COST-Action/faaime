// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <complex>
#include "nimble_fft.h"

// ---- 1. External Interface for Armadillo's Fast Fourier Transform Function ----
// A function the encapsulates the fft/ifft functions of the Armadillo library in a
// way that is accessible to be called from the nimbleExternalCall function from the
// nimble package.
extern "C" void nimble_fft(double *pReal, double *pImag, int iLength, bool bInverse) {
  try {
    const std::size_t tLength = std::size_t(iLength);
    double dTempReal, dTempImag;
    arma::cx_vec vCmplx(tLength);
    std::complex<double> cVal;
    if(iLength <= 0) {
      Rcpp::stop("length must be strictly positive");
    }
    // Copy across the input elements across to the complex vector required by
    // the armadillo interface
    for(size_t tIter = 0; tIter < tLength; ++tIter) {
      dTempReal = dTempImag = 0.0;
      if(pReal != NULL) {
        dTempReal = pReal[tIter];
      }
      if(pImag != NULL) {
        dTempImag = pImag[tIter];
      }
      vCmplx(tIter) = std::complex<double>(dTempReal, dTempImag);
    }
    // Call the forward or inverse Fast Fourier Transform from armadillo
    // depending on the value of the bInverse argument
    if(bInverse) {
      vCmplx = arma::ifft(vCmplx);
    } else {
      vCmplx = arma::fft(vCmplx);
    }
    // Copy across the outputs from the FFT operation into the relevant arrays
    for(size_t tIter = 0; tIter < tLength; ++tIter) {
      cVal = vCmplx(tIter);
      pReal[tIter] = cVal.real();
      pImag[tIter] = cVal.imag();
    }
  } catch(...) {
    // Avoid any exceptions being thrown outside of the external interface
    Rcpp::stop("exception encountered during processing");
  }
  return;
}
