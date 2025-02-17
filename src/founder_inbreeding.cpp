#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector v_prior_with_f_cpp(IntegerVector founder_masks,
                                 NumericVector founder_f,
                                 int number_of_transmissions,
                                 int number_of_fixed_transmissions) {
  if (founder_masks.size() != founder_f.size()){
    Rcpp::stop("founder_masks and founder_f need to have the same length");
  }

  int number_of_non_fixed_transmissions = number_of_transmissions -
    number_of_fixed_transmissions;

  // start from uniform prior on transmission vectors
  int number_of_v = 1 << number_of_non_fixed_transmissions;
  NumericVector v_prior(number_of_v, 1.0 / number_of_v);

  // bias prior for every inbred founder
  for (int i_founder = 0; i_founder < founder_f.size(); i_founder++){
    double f = founder_f[i_founder];

    if (f > 0){
      // bias the probability of the other transmissions
      // from this founder such that those are more likely
      // to be ibd
      int mask = founder_masks[i_founder];
      double original = std::pow(0.5, __builtin_popcount(mask));

      for (int v = 0; v < number_of_v; v++){
        if ((v & mask) == 0){
          v_prior[v] = f * (v_prior[v] / original) +
            (1-f) * v_prior[v];
        }else{
          v_prior[v] = (1-f) * v_prior[v];
        }
      }
    }
  }

  return v_prior;
}
