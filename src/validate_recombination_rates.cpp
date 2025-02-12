#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void validate_recombination_rates_cpp(NumericVector x) {

  for (int i_locus = 0; i_locus < x.size(); i_locus++){
    if (x[i_locus] < 0){
      Rcpp::stop("recombination rates should not be negative");
    }
    if (x[i_locus] > 0.5){
      Rcpp::stop("recombination rates should not exceed 0.5");
    }
  }
}
