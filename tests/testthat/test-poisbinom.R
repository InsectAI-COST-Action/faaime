# Create a vector of probabilities of independent trials for use in testing
testprobvec <- seq(0.0, 1.0, by = 0.05)
# Set the maximum tolerance for differences in probabilities from the benchmark
baseeps <- 0.00001
probeps <- (1.0 / length(testprobvec)) * baseeps

# Compile the poisbinom_massvec function
compiled_massvec <- nimble::compileNimble(poisbinom_massvec)
# Calculate the Poisson binomial mass probabilities according to the faaime package
faaime_calc <- compiled_massvec(testprobvec)

test_that("calculated Poisson binomial mass probabilities sum to one", {
  expect_lt(abs(sum(faaime_calc) - 1.0), baseeps)
})

test_that("calculated Poisson binomial mass probabilities are consistent with the poibin package implementation", {
  # Calculate the Poisson binomial mass probabilities according to the poibin package
  poibin_calc <- poibin::dpoibin(0:length(testprobvec), testprobvec)
  # Test to see if the calculates probabilities are close enough (within a small error due to different FFT algorithms)
  expect_all_true(abs(poibin - faaime_calc) < probeps)
})
