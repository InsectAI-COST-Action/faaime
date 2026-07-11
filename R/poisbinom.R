# 1.1 ---- atan2 function for NIMBLE ----
#' NIMBLE Implementation of arctan2 function
#'
#' A function for NIMBLE to allow the calculation of the arctan2 function. This
#' function has not been vectorised.
#'
#' @param iny Numeric scalar
#' @param inx Numeric scalar
#'
#' @return A numeric scalar returning the angle (in radians) between the x-axis
#' and the vector from the origin to (\code{inx}, \code{iny}).
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleExternalCall}}
#' @export
nimatan2 <- nimble::nimbleFunction(run = function(iny = double(0), inx = double(0)) {
  returnType(double(0))
  outVal <- 0.0
  if(inx > 0.0) {
    outVal <- atan(iny / inx)
  } else if(inx < 0.0) {
    if(iny >= 0.0) {
      outVal <- pi + atan(iny / inx)
    } else {
      outVal <- -pi + atan(iny / inx)
    }
  } else if(iny > 0.0) {
    outVal <- pi / 2.0
  } else if(iny < 0.0) {
    outVal <- -pi / 2.0
  }
  return(outVal)
})

# 1.2 ---- Calculate Vector of Probability Mass Values for Poisson Binomial Distribution ----
#' Calculate Vector of Poisson Binomial Distribution Values
#'
#' Function that provides an interface to the NIMBLE modelling framework to
#' create a vector of probability mass values resulting from a Poisson Binomial
#' distribution.
#'
#' @param probVec A numeric vector of probability values for each of the
#' independant trials that comprise the Poisson binomial distribution.
#'
#' @return A numeric vector with length \code{length(probVec) + 1} containing
#' the probability masses from the corresponding Poisson binomial distribution.
#' The probability of exactly n successful trails can be found in the (n + 1)th
#' element of the output vector.
#'
#' @details The calculation of the probability masses is performed according to
#' the discrete Fourier transform of the Poisson binomial characteristic
#' function as described by Hong (2013).
#'
#' @references Hong (2013) On computing the distribution function for the
#' Poisson binomial distribution *Computational Statistics and Data Analysis*
#' 59, 41-51.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleExternalCall}}
#' @export
poisbinom_massvec <- nimble::nimbleFunction(run = function(probVec = double(1)) {
  returnType(double(1))
  vecLength <- length(probVec)
  realPart <- numeric(vecLength + 1)
  imagPart <- numeric(vecLength + 1)
  realPart[1] <- 1.0
  imagPart[1] <- 0.0
  # Sanity check the inputs
  if(any(probVec < 0.0 | probVec > 1.0)) {
    stop("invalid probabilitiy value input")
  }
  for(iIter in 2:(ceiling(vecLength / 2.0) + 1)) {
    c1 <- 0.0
    c2 <- 0.0
    tt <- (iIter - 1) * 2.0 * pi / (vecLength + 1.0)
    for(jIter in 1:vecLength) {
      pj <- probVec[jIter]
      ax <- 1.0 - pj + pj * cos(tt)
      bx <- pj * sin(tt)
      tmpOne <- sqrt(ax * ax + bx * bx)
      tmpTwo <- nimatan2(bx, ax)
      c1 <- c1 + log(tmpOne)
      c2 <- c2 + tmpTwo
    }
    a2 <- exp(c1) * cos(c2)
    b2 <- exp(c1) * sin(c2)
    realPart[iIter] <- a2
    imagPart[iIter] <- b2
    realPart[vecLength + 3 - iIter] <- a2
    imagPart[vecLength + 3 - iIter] <- -b2
  }
  # Normalise the probabilities and set any negative values to zero (negative values will only come about from underflow errors during the FFT)
  outProbs <- pmax(nimfft(realPart, imagPart, FALSE)[1, ] / (vecLength + 1), rep(0.0, vecLength + 1))
  return(outProbs)
})

# 1.3 ---- Probability Mass Function for Poisson Binomial Distribution ----
#' Probability Mass Function for the Poisson Binomial Distribution
#'
#' Function that provides the probability mass function of the Poisson binomial
#' distribution.  This function is encapsulated in a way that it can extend the
#' functionality of the NIMBLE modelling framework.
#'
#' @param x Integer scalar that contains the number of successful trials from
#' the Poisson Binomial process.
#' @param probs Numeric vector containing the probability of success of each of
#' the independent trails that comprise the Poisson Binomial distribution.
#' @param log Logical scalar indicating whether the log of the probability mass
#' is required.
#'
#' @return A numeric scalar containing the probability mass (or the log of the
#' probability mass if \code{log = TRUE}) of obtaining exactly \code{x}
#' successful trials.
#'
#' @details The calculation of the probability masses is performed according to
#' the discrete Fourier transform of the Poisson binomial characteristic
#' function as described by Hong (2013).
#'
#' @references Hong (2013) On computing the distribution function for the
#' Poisson binomial distribution *Computational Statistics and Data Analysis*
#' 59, 41-51.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleExternalCall}} \code{\link{poisbinom_massvec}}
#' @export
dpoisbinom <- nimble::nimbleFunction(run = function(x = integer(0), probs = double(1), log = integer(0, default = 0)) {
  returnType(double(0))
  # Retrieve the number of trials
  numTrials <- length(probs)
  tarProb <- -Inf
  if(x >= 0 & x <= numTrials) {
    # Retrieve the vector of mass probabilities from convolution of the characteristic function
    outProbs <- poisbinom_massvec(probs)
    if(log) {
      tarProb <- log(outProbs[x + 1])
    } else {
      tarProb <- outProbs[x + 1]
    }
  } else if(!log) {
    tarProb <- 0.0
  }
  return(tarProb)
})

# 1.4 ---- Random Number Generator for the Poisson Binomial Distribution ----
#' Random Number Generator for the Poisson Binomial Distribution
#'
#' Generate a random number from the Poisson Binomial Distribution. This
#' function is encapsulated in a way that it can extend the functionality of the
#' NIMBLE modelling framework.
#'
#' @param n Integer scalar indicating the number of random samples to generate
#' (currently this value must be 1).
#' @param probs Numeric vector containing the probability of success of each of
#' the independent trails that comprise the Poisson Binomial distribution.
#'
#' @return An integer scalar containing a random draw of the number of
#' successful trials.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleExternalCall}}
#' @export
rpoisbinom <- nimble::nimbleFunction(run = function(n = integer(0), probs = double(1)) {
  returnType(integer(0))
  if(n != 1) {
    print("rpoisbinom only allows n = 1; using n = 1.")
  }
  numSuccess <- 0
  for(probIter in 1:length(probs)) {
    numSuccess <- numSuccess + rbinom(1, 1, probs[probIter])
  }
  return(numSuccess)
})

# 1.5 ---- Cumulative Probability Mass Function for the Poisson Binomial Distribution ---
#' Cumulative Probability Mass Function for the Poisson Binomial Distribution
#'
#' Calculate the cumulative probability mass for the Poisson binomial
#' distribution. This function is encapsulated in a way that it can extend the
#' functionality of the NIMBLE modelling framework.
#'
#' @param q Integer scalar containing the number of successful trials to
#' calculate the cumulative probability mass for.
#' @param probs Numeric vector containing the probability of success of each of
#' the independent trails that comprise the Poisson Binomial distribution.
#' @param lower.tail Logical scalar. If \code{TRUE} calculate the probability
#' that there \code{q} or fewer successful trials.  If \code{FALSE} calculate
#' the complementary cumulative distribution function.
#' @param log.p Logical scalar indicating whether the log of the cumulative
#' probability mass is required.
#'
#' @return A numeric scalar containing the cumulative probability.
#'
#' @details The calculation of the probability masses is performed according to
#' the discrete Fourier transform of the Poisson binomial characteristic
#' function as described by Hong (2013).
#'
#' @references Hong (2013) On computing the distribution function for the
#' Poisson binomial distribution *Computational Statistics and Data Analysis*
#' 59, 41-51.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleExternalCall}} \code{\link{poisbinom_massvec}}
#' @export
ppoisbinom <- nimble::nimbleFunction(run = function(q = integer(0), probs = double(1), lower.tail = integer(0, default = 1), log.p = integer(0, default = 0)) {
  returnType(double(0))
  numTrials <- length(probs)
  tarProb <- 0.0
  if(q >= 0 & q < numTrials) {
    # Retrieve the vector of mass probabilities from convolution of the characteristic function
    outProbs <- poisbinom_massvec(probs)
    tarProb <- sum(outProbs[1:(q + 1)])
  } else if(q >= numTrials) {
    tarProb <- 1.0
  }
  if(!lower.tail) {
    tarProb <- 1.0 - tarProb
  }
  if(log.p) {
    tarProb <- log(tarProb)
  }
  return(tarProb)
})

# 1.6 ---- Quantile Function for the Poisson Binomial Distribution ----
#' Quantile Function for the Poisson Binomial Distribution
#'
#' Calculate the quantiles for the Poisson binomial distribution.  This function
#' is encapsulated in a way that it can extend the functionality of the NIMBLE
#' modelling framework.
#'
#' @param p A numeric scalar representing the cumulative probability mass for
#' which the corresponding quantile is going to be calculated.
#' @param probs Numeric vector containing the probability of success of each of
#' the independent trails that comprise the Poisson Binomial distribution.
#' @param lower.tail Logical scalar. If \code{TRUE} then the quantile is being
#' calculated from the lower tail of the probability distribution.  If
#' \code{FALSE} then the quantile is being calculated from the upper tail of the
#' probability distribution.
#' @param log.p Logical scalar indicating whether the cumulative probability
#' mass is given on a log scale.
#'
#' @return The number of successful trials associated with the relevant
#' quantile.
#'
#' @details The calculation of the probability masses is performed according to
#' the discrete Fourier transform of the Poisson binomial characteristic
#' function as described by Hong (2013).
#'
#' @references Hong (2013) On computing the distribution function for the
#' Poisson binomial distribution *Computational Statistics and Data Analysis*
#' 59, 41-51.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[nimble]{nimbleExternalCall}} \code{\link{poisbinom_massvec}}
#' @export
qpoisbinom <- nimble::nimbleFunction(run = function(p = double(0), probs = double(1), lower.tail = integer(0, default = 1), log.p = integer(0, default = 0)) {
  returnType(integer(0))
  numTrials <- length(probs)
  nOut <- 0
  tarp <- p
  if(log.p) {
    tarp <- exp(p)
  }
  if(!lower.tail) {
    tarp <- 1.0 - tarp
  }
  if(tarp > 1.0) {
    nOut <- numTrials
  } else if(tarp > 0.0) {
    outProbs <- poisbinom_massvec(probs)
    startp <- outProbs[1]
    while(tarp > startp) {
      nOut <- nOut + 1
      startp <- startp + outProbs[nOut + 1]
    }
  }
  return(nOut)
})
