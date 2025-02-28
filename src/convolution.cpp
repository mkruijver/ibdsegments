#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int find_index_of_last_non_zero(NumericVector x, double eps) {

  for (int i = x.size() - 1; i >=0; i--){
    if (x[i] > eps){
      return i + 1;
    }
  }

  return 1;
}

// [[Rcpp::export]]
int find_index_of_first_non_zero(NumericVector x, double eps) {

  int n = x.size();
  for (int i = 0; i < n; i++){
    if (x[i] > eps){
      return i + 1;
    }
  }

  return n + 1;
}

// [[Rcpp::export]]
NumericVector cumulative_simpson_cpp(NumericVector fx) {

  // this is ported from distr:::.csimpsum

  int n = fx.size();

  if (n % 2 != 1) {
    Rcpp::stop("fx needs to have an odd number of points");
  }

  int n_half = n/2; // integer division
  NumericVector ff(n_half + 1);

  double sum_odd = 0;
  double sum_even = fx[0];

  for (int i = 1; i <= n_half; ++i) {
    double f_even = fx[2 * i];
    double f_odd = fx[2 * i - 1];

    sum_odd += f_odd;
    sum_even += f_even;

    ff[i] = (2 * sum_even - f_even - fx[0] + 4 * sum_odd) / 3.0;
  }

  return ff;
}

// [[Rcpp::export]]
List pmf_of_sum(NumericVector x1, NumericVector p1,
                NumericVector x2, NumericVector p2,
                double eps) {

  int n1 = x1.size();
  if (p1.size() != n1) Rcpp::stop("p1 does not have the same length as x1");

  int n2 = x2.size();
  if (p2.size() != n2) Rcpp::stop("p2 does not have the same length as x2");

  std::map<double, double> sum_pmf;

  for (int i = 0 ; i < n1; i++){
    for (int j = 0; j < n2; j++){
      double sum_value = x1[i] + x2[j];
      double prob = p1[i] * p2[j];

      if (prob >= eps){
        sum_pmf[sum_value] += prob;
      }
    }
  }

  NumericVector sum_x;
  NumericVector sum_p;

  for (const auto& x_px_pair : sum_pmf) {
    sum_x.push_back(x_px_pair.first);
    sum_p.push_back(x_px_pair.second);
  }

  return DataFrame::create(Named("x") = sum_x,
                           Named("px") = sum_p);
}
