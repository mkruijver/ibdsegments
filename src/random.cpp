#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List random_ibd(int n,
                NumericVector chromosome_length,
                const IntegerVector ibd_state_by_v,
                int number_of_transmissions,
                IntegerVector fixed_transmission_masks,
                int state_stats = 1,
                bool merge_ibd = true) {

  int number_of_fixed_transmissions = fixed_transmission_masks.size();
  int number_of_non_fixed_transmissions = number_of_transmissions -
    number_of_fixed_transmissions;

  int number_of_inheritance_vectors = 1 << (number_of_non_fixed_transmissions);

  int number_of_chromosomes = chromosome_length.size();

  // data structures for collecting sampling results
  CharacterVector col_names {"sample", "chromosome", "start", "end",
                             "length", "state"};

  std::vector<int> sample_number;
  std::vector<int> chromosome;
  std::vector<double> start;
  std::vector<double> end;
  std::vector<double> length;
  std::vector<double> state;

  // data structures for collecting summary stats by sample
  std::vector<double> stats_total_length_in_state;
  std::vector<int> stats_number_of_segments_in_state;
  stats_total_length_in_state.resize(n);
  stats_number_of_segments_in_state.resize(n);

  // determine transmissions that can be ignored
  // because they are fixed and don't affect other transmissions
  LogicalVector transmission_is_ignored(number_of_transmissions);
  int number_ignored = 0;
  for (int i_fixed = 0; i_fixed < number_of_fixed_transmissions; i_fixed++){
    if (fixed_transmission_masks[i_fixed] == 0){
      // this transmission can be ignored
      transmission_is_ignored[number_of_non_fixed_transmissions + i_fixed] = true;
      number_ignored++;
    }
  }

  if (number_ignored == number_of_transmissions){
    Rcpp::stop("all transmissions are ignored");
  }

  // collect indices of not ignored transmissions
  std::vector<int> idx_transmissions_not_ignored;
  int number_of_transmissions_not_ignored = number_of_transmissions - number_ignored;

  idx_transmissions_not_ignored.reserve(number_of_transmissions_not_ignored);
  for (int i = 0; i < transmission_is_ignored.size(); i++){
    if (!transmission_is_ignored[i]){
      idx_transmissions_not_ignored.push_back(i);
    }
  }

  double rate = 0.01 * (number_of_transmissions_not_ignored);
  double mu = 1.0 / rate;

  for (int i = 0; i < n; i++){

    double total_length_in_state = 0;
    int number_of_segments_in_state = 0;

    for (int i_chromosome = 0; i_chromosome < number_of_chromosomes; i_chromosome++){
      double L = chromosome_length[i_chromosome];

      double segment_start = 0;
      double sub_segment_start = 0;

      int segment_v = sample(number_of_inheritance_vectors, 1,
      true, R_NilValue, false)[0];
      // TODO: below is safe unless number_of_inheritance_vectors is too large
      // int segment_v = std::floor(R::runif(0, number_of_inheritance_vectors));

      bool end_reached = false;

      while (!end_reached){
        double l = R::rexp(mu);
        double segment_end = sub_segment_start + l;

        if (segment_end > L){
          end_reached = true;
          segment_end = L;
        }

        int v_next;

        int segment_ibd = ibd_state_by_v[segment_v];
        bool ibd_status_changed = false;

        if (!end_reached){
          // flip a bit
          int i_transmission_not_ignored = std::floor(
                R::runif(0, number_of_transmissions_not_ignored));
          int i_transmission = idx_transmissions_not_ignored[i_transmission_not_ignored];

          if (i_transmission < number_of_non_fixed_transmissions){
            int flip = 1 << i_transmission;
            v_next = segment_v ^ flip;
          }
          else{
            int i_fixed_transmission = i_transmission -
              number_of_non_fixed_transmissions;
            int flip = fixed_transmission_masks[i_fixed_transmission];

            v_next = segment_v ^ flip;
          }

          int ibd_next = ibd_state_by_v[v_next];
          ibd_status_changed = segment_ibd != ibd_next;
        }

        bool write_segment = ibd_status_changed || (!merge_ibd) ||
          (end_reached);

        if (write_segment){
          double segment_length = segment_end - segment_start;

          sample_number.push_back(i + 1);
          chromosome.push_back(i_chromosome + 1);
          start.push_back(segment_start);
          end.push_back(segment_end);
          length.push_back(segment_length);
          state.push_back(segment_ibd);

          // collect state
          if (segment_ibd == state_stats){
            number_of_segments_in_state++;
            total_length_in_state += segment_length;
          }

          segment_start = segment_end;
          sub_segment_start = segment_end;
        }
        else{
          sub_segment_start += l;
        }

        segment_v = v_next;
      }
    }

    // collect stats
    stats_total_length_in_state[i] = total_length_in_state;
    stats_number_of_segments_in_state[i] = number_of_segments_in_state;
  }

  List samples = List::create(sample_number, chromosome, start, end, length, state);
  samples.attr("names") = col_names;
  samples.attr("row.names") = Rcpp::seq(1, sample_number.size());
  samples.attr("class") = "data.frame";

  List stats = List::create(Named("total_length") = stats_total_length_in_state,
                            Named("segments") = stats_number_of_segments_in_state);
  stats.attr("row.names") = Rcpp::seq(1, n);
  stats.attr("class") = "data.frame";

  return List::create(Named("samples") = samples,
                      Named("stats") = stats);
}
