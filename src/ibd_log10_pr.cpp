#include <Rcpp.h>
using namespace Rcpp;

int get_number_of_flipped_bits(int x, int y){
  unsigned int flipped = x ^ y;

  // count bit by bit
  int count = 0;

  while (flipped) {
    // check last bit
    count += flipped & 1;

    // shift right
    flipped >>= 1;
  }

  return count;
}

int count_number_of_set_bits(int v){
  return get_number_of_flipped_bits(0, v);
}

// [[Rcpp::export]]
NumericVector FFT_p(NumericVector p, int number_of_bits) {
  // computes the FFT of p0 as described by Kruglyak

  int size = p.size();

  NumericVector p0 = Rcpp::clone(p);
  NumericVector p1(size);

  bool toggle = false;

  for(int i = 0; i < number_of_bits; i++){

    for(int v = 0; v < size; v++){

      int shift = 1 << i;
      NumericVector& src = toggle ? p1 : p0;
      NumericVector& dst = toggle ? p0 : p1;

      dst[v] = ((v & shift) ? -src[v] : src[v]) + src[v ^ shift];
    }

    toggle = !toggle;
  }

  if (toggle) {
    return(p1);
  }
  else{
    return(p0);
  }
}

// [[Rcpp::export]]
NumericVector FFT_T(int number_of_v, double theta, IntegerVector fixed_transmission_masks,
                    int number_of_transmissions, int number_of_fixed_transmissions){

  // computes the FFT of the transition probabilities as described Kruglyak
  NumericVector ret(number_of_v);

  double x = 1.0 - 2.0 * theta;

  NumericVector pows(number_of_transmissions + 1);
  for(int i = 0; i < number_of_transmissions + 1; i++) pows[i] = pow(x,i);

  for(int v = 0; v < number_of_v; v++){
    int pi_v = v;

    for(int i = 0; i < number_of_fixed_transmissions; i++){
      int bits_from_i = v & fixed_transmission_masks[i];

      pi_v = pi_v ^ (count_number_of_set_bits(bits_from_i) % 2 <<
        (number_of_transmissions - number_of_fixed_transmissions + i));
    }

    ret[v] = pows[count_number_of_set_bits(pi_v)];
  }

  return(ret);
}

// [[Rcpp::export]]
NumericVector ibd_log10_pr_cpp(IntegerVector ibd_state_by_v,
                               IntegerVector ibd_by_locus,
                               NumericVector recombination_rate_by_locus,
                               int number_of_transmissions,
                               IntegerVector fixed_transmission_masks,
                               double pr_v_constant,
                               NumericVector pr_v_biased){

  int number_of_fixed_transmissions = fixed_transmission_masks.size();
  int number_of_bits = number_of_transmissions - number_of_fixed_transmissions;

  double log10_pr = 0;
  NumericVector log10_pr_by_locus(ibd_by_locus.size(), -INFINITY);

  // assign uniform prior on transmission vectors
  NumericVector v_prior;
  if (pr_v_constant >= 0){
    v_prior = NumericVector(ibd_state_by_v.size(), pr_v_constant);
  }
  else{
    v_prior = pr_v_biased;
  }

  // reserve memory for a posterior
  NumericVector v_posterior(ibd_state_by_v.size());

  NumericVector v_pr(ibd_state_by_v.size());

  int number_of_loci = ibd_by_locus.size();

  for (int i_locus = 0; i_locus < number_of_loci; i_locus++){
    double locus_pr = 0;

    int ibd_locus = ibd_by_locus[i_locus];

    if (i_locus > 0){

      double theta = recombination_rate_by_locus[i_locus - 1];
      bool recombines = theta > 0;

      if (recombines){
        NumericVector v_fft = FFT_p(v_posterior, number_of_bits);

        double v_fft_sum = Rcpp::sum(Rcpp::abs(v_fft));

        NumericVector T_fft = FFT_T(ibd_state_by_v.size(),
                                    theta,
                                    fixed_transmission_masks,
                                    number_of_transmissions,
                                    number_of_fixed_transmissions);

        // multiply the transforms
        for (int i = 0; i < v_fft.size(); i++){
          v_fft[i] *= T_fft[i] / v_fft_sum;
        }

        // transform back
        NumericVector v = FFT_p(v_fft, number_of_bits);

        // normalise
        double v_sum = Rcpp::sum(Rcpp::abs(v));

        for (int i_v = 0; i_v < v_prior.size(); i_v++){
          v_prior[i_v] = v[i_v] / v_sum;
        }

      }else{
        for (int i_v = 0; i_v < v_prior.size(); i_v++){
          v_prior[i_v] = v_posterior[i_v];
        }
      }
    }

    // add to the locus pr
    for (int i_v = 0; i_v < v_prior.size(); i_v++){
      double pr = ibd_locus == ibd_state_by_v[i_v] ?  v_prior[i_v] : 0;

      v_pr[i_v] = pr;

      locus_pr += pr;
    }

    // compute a posterior on transmission vectors
    for (int i_v = 0; i_v < v_prior.size(); i_v++){
      v_posterior[i_v] = v_pr[i_v] / locus_pr;
    }

    double log10_pr_locus = std::log10(locus_pr);

    log10_pr += log10_pr_locus;
    log10_pr_by_locus[i_locus] = log10_pr_locus;

    // a zero-likelihood knocks out all other locus contributions
    if (locus_pr == 0) return(log10_pr_by_locus);
  }

  return log10_pr_by_locus;
}
