#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List identificationROC(NumericVector vprobs, LogicalVector vcalli) {
  NumericMatrix curveCoords;
  NumericVector uniqueProbs, allProbs(vprobs.size() + 2), aucVal(1, 0.0);
  std::copy(vprobs.begin(), vprobs.end(), allProbs.begin());
  allProbs[vprobs.size() + 1] = 0.0;
  allProbs[vprobs.size() + 2] = 1.0;
  uniqueProbs = unique(allProbs);
  return List::create(
    Named("curveCoords") = curveCoords,
    Named("AUC") = aucVal
  );
}
