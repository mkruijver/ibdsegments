#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void one_p_step(const NumericVector& a_prev,
                NumericVector& a_next,
                const int stay_in_ibd,
                const IntegerVector &ibd_state_by_v,
                const int number_of_transmissions,
                const IntegerVector& masks){

  a_next.fill(0.0);
  double pr = 1.0 / number_of_transmissions;

  int number_of_fixed_transmissions = masks.size();

  int number_of_non_fixed_transmissions = number_of_transmissions -
    number_of_fixed_transmissions;

  int v_n = a_prev.size();

  for (int v = 0; v < v_n; v++){

    if (a_prev[v] > 0){
      // flip all bits
      for (int i = 0; i < number_of_non_fixed_transmissions; i++){

        int w = v ^ (1 << i);

        if (ibd_state_by_v[w] == stay_in_ibd){
          a_next[w] += a_prev[v] * pr;
        }
      }

      // founder symmetry
      for (int i = 0; i < number_of_fixed_transmissions; i++){
        int w = v ^ masks[i];

        if (ibd_state_by_v[w] == stay_in_ibd){
          a_next[w] += a_prev[v] * pr;
        }
      }
    }
  }
}

NumericVector p_step(const NumericVector& a0,
                     const double lambda,
                     const IntegerVector& ibd_state_by_v,
                     const int stay_in_ibd,
                     const IntegerVector& masks,
                     const int number_of_transmissions,
                     const int number_of_fixed_transmissions) {

  int number_of_non_fixed_transmissions = number_of_transmissions - number_of_fixed_transmissions;
  double pr = 1.0 / number_of_transmissions;

  int min_number_of_crossovers = R::qpois(1e-16 / 2, lambda, 1, 0);
  int max_number_of_crossovers = R::qpois(1.0 - 1e-16, lambda, 1, 0);

  if (!std::isfinite(max_number_of_crossovers)){
    Rcpp::stop("max_number_of_crossovers is not finite");
  }

  int a_n = a0.size();
  NumericVector a_final(a_n);
  NumericVector a_prev(a_n);
  NumericVector a_next(a_n);

  // 0 steps
  for (int i = 0; i < a_n; i++){
    if (ibd_state_by_v[i] == stay_in_ibd){
      double pr_step =  R::dpois(0, lambda, false); //step_prs[0];
      a_final[i] = pr_step * a0[i];

      a_prev[i] = a0[i];
    }
  }

  double scale = 1.0;

  for (int number_of_steps = 1; number_of_steps < max_number_of_crossovers; number_of_steps++){

    one_p_step(a_prev, a_next, stay_in_ibd, ibd_state_by_v,
               number_of_transmissions, masks);

    // rescale a_next to avoid underflow in a_next
    double a_next_sum = 0.0;
    for (int i = 0; i < a_n; i++){
      a_next_sum += a_next[i];
    }
    scale *= a_next_sum;

    // if it is extremely improbable to remain in the state after a large
    // number of steps then we terminate the sum early
    if (scale < 1e-32){
      return(a_final);
    }

    for (int i = 0; i < a_n; i++){
      a_next[i] /= a_next_sum;
    }

    // add to the total
    if (number_of_steps >= min_number_of_crossovers){
      double const_term = R::dpois(number_of_steps, lambda, false) * scale;
      if (const_term > 0){
        for (int i = 0; i < a_n; i++){
          a_final[i] += a_next[i] * const_term;
        }
      }
    }

    std::swap(a_next, a_prev);
  }

  return a_final;
}

// [[Rcpp::export]]
NumericVector q_step(const NumericVector &a0, const int number_of_transmissions,
                     const int number_of_fixed_transmissions, const IntegerVector &masks) {

  int number_of_non_fixed_transmissions = number_of_transmissions - number_of_fixed_transmissions;

  double rate = 0.01;

  int number_of_v = a0.size();
  NumericVector a_next(number_of_v);

  for (int v = 0; v < number_of_v; v++){
    // flip all bits
    for (int i = 0; i < number_of_non_fixed_transmissions; i++){
      int w = v ^ (1 << i);

      a_next[w] += a0[v] * rate;
    }

    // founder symmetry
    for (int i = 0; i < number_of_fixed_transmissions; i++){
      int w = v ^ masks[i];

      a_next[w] += a0[v] * rate;
    }
  }

  return a_next;
}

double rate_instant_move(NumericVector alpha,
                         int from_ibd_state,
                         int number_of_transmissions,
                         IntegerVector masks,
                         IntegerVector ibd_state_by_v){

  int number_of_fixed_transmissions = masks.size();
  int number_of_non_fixed_transmissions = number_of_transmissions - number_of_fixed_transmissions;

  double q = 0.0;

  double rate = 0.01;

  double alpha_sum = 0;
  for (int i = 0; i < alpha.size(); i++){
    alpha_sum += alpha[i];
  }
  rate /= alpha_sum;

  for (int v = 0; v < alpha.size(); v++){
    // flip all bits
    for (int i = 0; i < number_of_non_fixed_transmissions; i++){
      int w = v ^ (1 << i);

      if (ibd_state_by_v[w] != from_ibd_state){
        q += alpha[v] * rate;
      }
    }

    // founder symmetry
    for (int i = 0; i < number_of_fixed_transmissions; i++){
      int w = v ^ masks[i];

      if (ibd_state_by_v[w] != from_ibd_state){
        q += alpha[v] * rate;
      }

    }
  }

  return q;
}

// [[Rcpp::export]]
void one_F_step(const NumericMatrix F,
                NumericMatrix F_next,
                const int stay_in_ibd,
                const IntegerVector ibd_state_by_v,
                const int number_of_transmissions,
                const int number_of_fixed_transmissions,
                const IntegerVector unique_masks,
                const IntegerVector unique_masks_count){

  F_next.fill(0.0);
  double pr = 1.0 / number_of_transmissions;

  int number_of_non_fixed_transmissions = number_of_transmissions -
    number_of_fixed_transmissions;

  int n_row = F.nrow();
  int n_col = F.ncol();

  int number_of_unique_masks = unique_masks.size();

  for (int v = 0; v < n_row; v++){
    for(int k = 0; k < n_col - 1; k++){
      double F_v_k = F(v, k);
      if (F_v_k == 0) continue;

      double F_v_k_times_pr = F_v_k * pr;

      // flip all bits
      for (int i = 0; i < number_of_non_fixed_transmissions; i++){

        int w = v ^ (1 << i);

        bool goes_to_ibd = ibd_state_by_v[w] == stay_in_ibd;
        int delta = (goes_to_ibd) ? 1 : 0;

        F_next(w, k + delta) += F_v_k_times_pr;
        // if (k ==0){
        //   Rcpp::Rcout << "v = " << v << " k = " << k << " w = " << w
        //   << " delta = " << delta << "\n";
        // }
      }

      // founder symmetry
      for (int i = 0; i < number_of_unique_masks; i++){
        int w = v ^ unique_masks[i];
        bool goes_to_ibd = ibd_state_by_v[w] == stay_in_ibd;
        int delta = (goes_to_ibd) ? 1 : 0;

        F_next(w, k + delta) += F_v_k_times_pr * unique_masks_count[i];

        // if (k ==0){
        //   Rcpp::Rcout << "v = " << v << " k = " << k << " w = " << w
        //               << " delta = " << delta << "\n";
        // }
      }
    }
  }
}

// [[Rcpp::export]]
List get_unique_masks_and_count(IntegerVector masks){

  std::vector<int> unique_masks;
  std::vector<int> unique_masks_count;
  for (int i_mask = 0; i_mask < masks.size(); i_mask++){
    int mask = masks[i_mask];
    bool found = false;

    for (int j = 0; j < unique_masks.size(); j++) {
      if (unique_masks[j] == mask) {
        unique_masks_count[j]++;
        found = true;
        break;
      }
    }

    if (!found) {
      unique_masks.push_back(mask);
      unique_masks_count.push_back(1);
    }
  }

  return(Rcpp::List::create(Rcpp::wrap(unique_masks),
                            Rcpp::wrap(unique_masks_count)));
}

// [[Rcpp::export]]
NumericMatrix pr_number_of_intervals_in_state_by_n(int ibd_state,
                                                   IntegerVector ibd_state_by_v,
                                                   int n_max,
                                                   int number_of_transmissions,
                                                   IntegerVector masks){
  //compute the probability distribution of spending k = 0, 1, 2, ... n_max+1
  //intervals in the ibd_state if there are n = 0, 1, 2, ..., n_max recombinations
  NumericMatrix V(n_max + 2, n_max + 1);

  int number_of_states = ibd_state_by_v.size();
  NumericMatrix F_current(number_of_states, n_max + 2);
  NumericMatrix F_next(number_of_states, n_max + 2);

  // each v is equally probable at start
  double pr_v = 1.0 / ibd_state_by_v.size();
  for (int v = 0; v < ibd_state_by_v.size(); v++){
    int delta = (ibd_state_by_v[v] == ibd_state) ? 1 : 0;
    F_current(v, delta) += pr_v;
    V(delta, 0) += pr_v;
  }

  // in many pedigrees the masks are not uniqe so we can optimise a bit
  // by only considering the unique ones in the computations
  List unique_masks_and_count = get_unique_masks_and_count(masks);
  IntegerVector unique_masks = unique_masks_and_count[0];
  IntegerVector unique_masks_count = unique_masks_and_count[1];

  // recursively compute F: joint pr. of being in state v (row) and
  // having spent m (col) intervals in the ibd state
  for (int n = 1; n <= n_max; n++){
    one_F_step(F_current, F_next, ibd_state, ibd_state_by_v,
               number_of_transmissions, masks.size(),
               unique_masks, unique_masks_count);

    // marginalise
    V(_, n) = Rcpp::colSums(F_next);

    std::swap(F_current, F_next);
  }

  return V;
}

// [[Rcpp::export]]
NumericVector pr_stay_and_leave(int stay_in_ibd,
                                int max_number_of_steps,
                                IntegerVector ibd_state_by_v,
                                int number_of_transmissions,
                                IntegerVector fixed_transmission_masks){

  NumericVector pr_stay_and_leave(max_number_of_steps +1);

  NumericVector a_prev(ibd_state_by_v.size());
  NumericVector a_next(ibd_state_by_v.size());

  // assign uniform prior on transmission vectors that match the ibd state
  int num_matching = 0;
  for (int v = 0; v < ibd_state_by_v.size(); v++){
    if (ibd_state_by_v[v] == stay_in_ibd) num_matching++;
  }
  for (int v = 0; v < ibd_state_by_v.size(); v++){
    if (ibd_state_by_v[v] == stay_in_ibd) {
      a_prev[v] = 1.0 / num_matching;
    }
  }

  double pr = 1.0 / number_of_transmissions;

  // 0 steps
  pr_stay_and_leave[0] = 1 * rate_instant_move(a_prev, stay_in_ibd,
                                               number_of_transmissions,
                                               fixed_transmission_masks,
                                               ibd_state_by_v);

  for (int i_xstep = 1; i_xstep <= max_number_of_steps; i_xstep++){
    one_p_step(a_prev, a_next, stay_in_ibd, ibd_state_by_v,
               number_of_transmissions,
               fixed_transmission_masks);

    for (int i = 0; i < a_next.size(); i++){
      pr_stay_and_leave[i_xstep] += a_next[i];
    }

    pr_stay_and_leave[i_xstep] *= rate_instant_move(a_next, stay_in_ibd,
                                                    number_of_transmissions,
                                                    fixed_transmission_masks,
                                                    ibd_state_by_v);

    std::swap(a_next, a_prev);
  }

  return(pr_stay_and_leave);
}

// [[Rcpp::export]]
double log10_ibd_segment_pr_cpp(NumericVector obs_cM,
                                IntegerVector obs_ibd,
                                IntegerVector ibd_state_by_v,
                                int number_of_transmissions,
                                IntegerVector fixed_transmission_masks){

  for (int i_segment = 0; i_segment < obs_cM.size(); i_segment++){
    if (obs_cM[i_segment] < 0) return std::log10(0);
  }

  int number_of_fixed_transmissions = fixed_transmission_masks.size();

  double log10_pr_total = 0.0;

  // assign uniform prior on transmission vectors
  int number_of_v = ibd_state_by_v.size();
  NumericVector v_prior(number_of_v, 1.0 / number_of_v);

  // reserve memory for a posterior
  NumericVector v_posterior(number_of_v);
  NumericVector v_pr(number_of_v);

  for (int i_segment = 0; i_segment < obs_ibd.size(); i_segment++){
    int obs_ibd_segment = obs_ibd[i_segment];
    double segment_pr = 0;

    // add to the segment pr
    for (int i_v = 0; i_v < number_of_v; i_v++){

      double pr = obs_ibd_segment == ibd_state_by_v[i_v] ? v_prior[i_v] : 0;
      v_pr[i_v] = pr;
      segment_pr += pr;
    }

    log10_pr_total += std::log10(segment_pr);

    // a zero-likelihood knocks out all other locus contributions
    if (segment_pr == 0) return(std::log10(0));

    // normalise v_pr to obtain the posterior pr's of v at the start of the segment
    for (int i_v = 0; i_v < number_of_v; i_v++){
      v_pr[i_v] /= segment_pr;
    }

    // stay in the same observed ibd status for the length of the segment
    double lambda = 0.01 * obs_cM[i_segment] * number_of_transmissions;

    NumericVector v_next = p_step(v_pr, lambda,
                                  ibd_state_by_v, obs_ibd_segment,
                                  fixed_transmission_masks, number_of_transmissions,
                                  number_of_fixed_transmissions);

    // pr of staying in this observed ibd state
    double v_next_sum = Rcpp::sum(v_next);
    log10_pr_total += std::log10(v_next_sum);

    bool is_last_segment = i_segment == obs_ibd.size() - 1;
    if (!is_last_segment){

      if (v_next_sum == 0){
        return std::log10(0);
      }

      // obtain posterior distribution of IBD vectors at the end of the segment
      for (int i_v = 0; i_v < number_of_v; i_v++){
        v_next[i_v] /= v_next_sum;
      }

      v_prior = q_step(v_next, number_of_transmissions,
                       number_of_fixed_transmissions, fixed_transmission_masks);

    }
  }

  return log10_pr_total;
}

// [[Rcpp::export]]
NumericVector log10_ibd_segment_pr_vectorised_cpp(
    IntegerVector sample, IntegerVector chromosome,
    NumericVector obs_cM, IntegerVector obs_ibd,
    IntegerVector ibd_state_by_v, int number_of_transmissions,
    IntegerVector fixed_transmission_masks)
{
  std::vector<double> log10_pr;

  int i_sample_start = 0;
  int i_chromosome_start = 0;

  int n_minus_one = sample.size() - 1;

  double log10_pr_sample = 0;

  for (int i = 0; i <= n_minus_one; i++){

    bool input_stop = i == n_minus_one;

    bool chromosome_stop = input_stop || (chromosome[i+1] != chromosome[i_chromosome_start]);
    bool sample_stop = input_stop || (sample[i+1] != sample[i_sample_start]);

    if (chromosome_stop || sample_stop){
      // calculate prob for this sample/chromosome combo
      NumericVector obs_cM_i = NumericVector(obs_cM.begin() + i_chromosome_start,
                                             obs_cM.begin() + i + 1);
      IntegerVector obs_ibd_i = IntegerVector(obs_ibd.begin() + i_chromosome_start,
                                              obs_ibd.begin() + i + 1);

      double log10_pr_i = log10_ibd_segment_pr_cpp(obs_cM_i, obs_ibd_i,
                                                   ibd_state_by_v, number_of_transmissions,
                                                   fixed_transmission_masks);

      log10_pr_sample += log10_pr_i;
      i_chromosome_start = i + 1;
    }
    if (sample_stop){
      i_sample_start = i + 1;

      log10_pr.push_back(log10_pr_sample);
      log10_pr_sample = 0;
    }
  }

  return wrap(log10_pr);
}

// [[Rcpp::export]]
double pr_never_in_state(NumericMatrix V,
                         NumericVector n_pr,
                         int n_max){
  // compute the pr. of not being in the IBD state at all

  if (n_pr.size() < (n_max + 1)){
    Rcpp::stop("n_pr needs to have a length of at least n_max + 1");
  }
  if (V.ncol() < (n_max + 1)){
    Rcpp::stop("V needs to have at least n_max + 1 columns");
  }
  if (V.nrow() < 1){
    Rcpp::stop("V needs to have at least 1 row");
  }

  double pr = 0;

  for (int n = 0; n <= n_max; n++){
    pr += n_pr[n] * V(0, n);
  }

  return pr;
}

// [[Rcpp::export]]
double pr_always_in_state(NumericMatrix V,
                          NumericVector n_pr,
                          int n_max){
  // compute the pr. of being in the IBD state the whole time

  if (n_pr.size() < (n_max + 1)){
    Rcpp::stop("n_pr needs to have a length of at least n_max + 1");
  }
  if (V.ncol() < (n_max + 1)){
    Rcpp::stop("V needs to have at least n_max + 1 columns");
  }

  if (V.nrow() < (n_max + 2)){
    Rcpp::stop("V needs to have at least n_max + 2 rows");
  }

  double pr = 0;

  for (int n = 0; n <= n_max; n++){
    pr += n_pr[n] * V(n+1, n);
  }

  return pr;
}

// [[Rcpp::export]]
NumericMatrix precompute_V_lbeta(NumericMatrix V){

  NumericMatrix V_lbetas(V.nrow(), V.ncol());


  for (int number_of_recombinations = 1;
       number_of_recombinations < (V.ncol()-1);
       number_of_recombinations++){

    for (int number_of_intervals_in_state = 1;
         number_of_intervals_in_state <= (number_of_recombinations+1) &&
           (number_of_intervals_in_state < V.nrow());
         number_of_intervals_in_state++){

      double pr_v = V(number_of_intervals_in_state,
                      number_of_recombinations);

      if (pr_v == 0) continue;

      double a = number_of_intervals_in_state;
      double b = number_of_recombinations - a + 1;

      double l_beta = R::lbeta(a,b);

      V_lbetas(number_of_intervals_in_state, number_of_recombinations) = l_beta;
    }
  }

  return V_lbetas;
}

// [[Rcpp::export]]
NumericVector d_fraction_ibd_state(NumericVector s,
                                    double L,
                                    NumericVector n_pr,
                                    NumericMatrix V,
                                    NumericMatrix V_lbeta,
                                    double point_mass_0,
                                    double point_mass_1){

  double conditional_fact = 1.0 - point_mass_0 - point_mass_1;

  int s_n = s.size();
  NumericVector f(s_n);
  NumericVector log_s(s_n);
  NumericVector log_1p_ms(s_n);

  // pre-compute some constants
  for(int i = 0; i < s_n; i++){
    double s_i = s[i];

    // smooth out the edges
    if ((s_i < 1e-16) && (s_i >= 0)){
      s[i] = 1e-16;
      s_i = s[i];
    }
    if ((s_i > 1 - 1e-16) && (s_i <= 1)){
      s[i] = 1 - 1e-16;
      s_i = s[i];
    }

    log_s[i] = std::log(s_i);
    log_1p_ms[i] = std::log1p(-s_i);
  }

  for (int number_of_recombinations = 1;
       number_of_recombinations < (n_pr.size() - 1);
       number_of_recombinations++){

    double pr_number_of_recombinations = n_pr[number_of_recombinations];

    for (int number_of_intervals_in_state = 1;
         number_of_intervals_in_state <= (number_of_recombinations+1) &&
           (number_of_intervals_in_state < V.nrow());
         number_of_intervals_in_state++){

      double constant = V(number_of_intervals_in_state,
                          number_of_recombinations);

      if (constant == 0) continue;

      constant *= pr_number_of_recombinations;
      constant /= conditional_fact;

      double a = number_of_intervals_in_state;
      double b = number_of_recombinations - a + 1;

      double l_beta = V_lbeta(number_of_intervals_in_state, number_of_recombinations);

      for (int i_s = 0; i_s < s_n; i_s++){
        double s_i = s[i_s];

        if ((s_i > 0) && (s_i < 1)){

          double conditional_pr_fraction_ibd = std::exp((a-1) * log_s[i_s] +
                                                        (b-1) * log_1p_ms[i_s] -
                                                        l_beta);

          f[i_s] += constant * conditional_pr_fraction_ibd;
        }
      }
    }
  }

  return f;
}
