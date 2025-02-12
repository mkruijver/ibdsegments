#include <Rcpp.h>
using namespace Rcpp;

const int COEFF_KAPPA = 2;
const int COEFF_IDENTITY = 9;
const int COEFF_DETAILED = 15;

// [[Rcpp::export]]
int get_kappa_state(IntegerVector x, int person_idx1, int person_idx2){
  int a = x[2 * person_idx1 - 2];
  int b = x[2 * person_idx1 - 1];

  int c = x[2 * person_idx2 - 2];
  int d = x[2 * person_idx2 - 1];

  if (a==b){
    int ibd = (a == c) + (a == d);

    return ibd;
  }
  else{
    int ibd = (a == c) + (a == d) - (a == c && c == d) +
      (b == c) + (b == d) - (b == c && c==d);

    return ibd;
  }
}

// [[Rcpp::export]]
int get_joint_kappa_state(IntegerVector x, IntegerVector persons_idx){
  if (persons_idx.size() < 2){
    Rcpp::stop("need at least two persons");
  }

  int a = x[2 * persons_idx[0] - 2];
  int b = x[2 * persons_idx[0] - 1];

  // check if a is present across all persons
  bool a_present_everywhere = true;
  bool b_present_everywhere = true;

  bool any_hets_in_other_persons = false;

  for (int i_person = 1; i_person < persons_idx.size(); i_person++){
    int c = x[2 * persons_idx[i_person] - 2];
    int d = x[2 * persons_idx[i_person] - 1];

    a_present_everywhere &= (c == a || d == a);
    b_present_everywhere &= (c == b || d == b);

    any_hets_in_other_persons |= (c != d);
  }

  if (a == b && a_present_everywhere && b_present_everywhere) {
    if (any_hets_in_other_persons){
      return 1;
    }else{
      return 2;
    }
  }

  if (a!=b){
    if (a_present_everywhere && b_present_everywhere) return 2;
    if (a_present_everywhere || b_present_everywhere) return 1;
  }

  return 0;
}

int get_comparison_mask(int a, int b, int c, int d){
  // we compare clockwise: a=b, a=d, a=c,
  //                       b=d, b=c, d=c

  int ab = (a == b);
  int ad = (a == d) << 1;
  int ac = (a == c) << 2;
  int bd = (b == d) << 3;
  int bc = (b == c) << 4;
  int cd = (c == d) << 5;

  int comparison_mask = ab + ad + ac + bd + bc + cd;

  return comparison_mask;
}

// [[Rcpp::export]]
int get_Jacquard_state(IntegerVector x, int person_idx1, int person_idx2){
  int a = x[2 * person_idx1 - 2];
  int b = x[2 * person_idx1 - 1];

  int c = x[2 * person_idx2 - 2];
  int d = x[2 * person_idx2 - 1];

  int comparison_mask = get_comparison_mask(a, b, c, d);
  switch(comparison_mask){
  case 0b111111:
    return(1);

  case 0b100001:
    return(2);

  case 0b010101:
  case 0b001011:
    return(3);

  case 0b000001:
    return(4);

  case 0b100110:
  case 0b111000:
    return(5);

  case 0b100000:
    return(6);

  case 0b001100:
  case 0b010010:
    return(7);

  case 0b000100:
  case 0b000010:
  case 0b010000:
  case 0b001000:
    return(8);

  case 0b000000:
    return(9);

  default:
    return -1;
  }
}

// [[Rcpp::export]]
int get_detailed_Jacquard_state(IntegerVector x, int person_idx1, int person_idx2){
  int a = x[2 * person_idx1 - 2];
  int b = x[2 * person_idx1 - 1];

  int c = x[2 * person_idx2 - 2];
  int d = x[2 * person_idx2 - 1];

  int comparison_mask = get_comparison_mask(a, b, c, d);
  switch(comparison_mask){

  case 0b111111:
    return(1); // 1
  case 0b010101:
    return(2);
  case 0b001011:
    return(3); // 3 -> 2 or 3
  case 0b100110:
    return(4);
  case 0b111000:
    return(5); // 5 -> 4 or 5
  case 0b100001:
    return(6); // 2 -> 6
  case 0b000001:
    return(7); // 4 -> 7
  case 0b100000:
    return(8); // 6 -> 8
  case 0b001100:
    return(9);
  case 0b010010:
    return(12); // 7 -> 9 or 12
  case 0b000100:
    return(10);
  case 0b001000:
    return(11);
  case 0b000010:
    return(13);
  case 0b010000:
    return(14); // 8 -> 10, 11, 13, 14
  case 0b000000:
    return(15);  // 9 -> 15
  default:
    return -1;
  }
}

// [[Rcpp::export]]
int get_ibd_state(IntegerVector x, int coeff, IntegerVector persons_idx){

  switch(coeff){
  case COEFF_KAPPA:
    if (persons_idx.size() == 2){
      return get_kappa_state(x, persons_idx[0], persons_idx[1]);
    }
    else{
      return get_joint_kappa_state(x, persons_idx);
    }
  case COEFF_IDENTITY:
    return get_Jacquard_state(x, persons_idx[0], persons_idx[1]);
  case COEFF_DETAILED:
    return get_detailed_Jacquard_state(x, persons_idx[0], persons_idx[1]);
  default:
    Rcpp::stop("Unknown ibd coeff");
  }
}

IntegerVector assign_founder_alleles(int number_of_persons,
                                     IntegerVector ped_row_is_founder_idx){

  IntegerVector x(2 * number_of_persons);

  for (int i_founder = 0; i_founder < ped_row_is_founder_idx.size(); i_founder++){
    int idx_1based = ped_row_is_founder_idx[i_founder];

    x[2 * idx_1based - 2] = 2 * i_founder + 1;
    x[2 * idx_1based - 1] = 2 * i_founder + 2;
  }

  return x;
}

void drop_founder_alleles(IntegerVector x,
                          int v,
                          IntegerVector from_allele_idx,
                          IntegerVector to_allele_idx,
                          IntegerVector top_to_bottom_order){

  for (int i_order = 0; i_order < top_to_bottom_order.size(); i_order++){

    int i_transmission = top_to_bottom_order[i_order] - 1;

    int from_idx = from_allele_idx[i_transmission] + ((v >> i_transmission) & 1);
    int to_idx = to_allele_idx[i_transmission];

    x[to_idx - 1] = x[from_idx -1];
  }
}

// [[Rcpp::export]]
CharacterVector get_alleles_for_v(int v,
                                int number_of_ped_members,
                                IntegerVector ped_row_is_founder_idx,
                                IntegerVector from_allele_idx,
                                IntegerVector to_allele_idx,
                                IntegerVector persons_idx,
                                int number_of_fixed_transmissions,
                                IntegerVector top_to_bottom_order){

  // reserve return value
  int number_of_transmissions = from_allele_idx.size();
  int number_of_non_fixed_transmissions = number_of_transmissions - number_of_fixed_transmissions;

  int number_of_canonical_inheritance_vectors = 1 << number_of_non_fixed_transmissions;

  if (v < 0 || v >= number_of_canonical_inheritance_vectors){
    Rcpp::stop("invalid v: ", v);
  }

  // assign allele vector
  IntegerVector x = assign_founder_alleles(number_of_ped_members, ped_row_is_founder_idx);
  drop_founder_alleles(x, v, from_allele_idx, to_allele_idx, top_to_bottom_order);

  // convert to CharacterVector for display
  CharacterVector x1(number_of_ped_members);
  for (int i_ped_member = 0; i_ped_member < number_of_ped_members; i_ped_member++){
    int a = x[2 * i_ped_member];
    int b = x[2 * i_ped_member + 1];

    x1[i_ped_member] =  std::to_string(a) + " " + std::to_string(b);
  }

  return x1;
}

// [[Rcpp::export]]
IntegerVector get_ibd_states_by_v(int number_of_ped_members,
                                  IntegerVector ped_row_is_founder_idx,
                                  IntegerVector from_allele_idx,
                                  IntegerVector to_allele_idx,
                                  IntegerVector persons_idx,
                                  int number_of_fixed_transmissions,
                                  IntegerVector top_to_bottom_order,
                                  int coeff){

  // reserve return value
  int number_of_transmissions = from_allele_idx.size();
  int number_of_non_fixed_transmissions = number_of_transmissions - number_of_fixed_transmissions;

  int number_of_canonical_inheritance_vectors = 1 << number_of_non_fixed_transmissions;;
  IntegerVector ibd_states(number_of_canonical_inheritance_vectors);

  // assign allele vector
  IntegerVector x = assign_founder_alleles(number_of_ped_members, ped_row_is_founder_idx);

  // drop alleles for each choice of inheritance vector
  for (int v = 0; v < number_of_canonical_inheritance_vectors; v++){

    // drop alleles down the pedigree for this inheritance vector
    drop_founder_alleles(x, v, from_allele_idx, to_allele_idx, top_to_bottom_order);

    // obtain ibd state for this inheritance vector
    ibd_states[v] = get_ibd_state(x, coeff, persons_idx);
  }

  return ibd_states;
}


// [[Rcpp::export]]
IntegerMatrix get_multi_ibd_patterns_by_v(int number_of_ped_members,
                                          IntegerVector ped_row_is_founder_idx,
                                          IntegerVector from_allele_idx,
                                          IntegerVector to_allele_idx,
                                          IntegerVector person_ids,
                                          int number_of_fixed_transmissions,
                                          IntegerVector top_to_bottom_order,
                                          bool minimal = true){
  // reserve return value
  int number_of_transmissions = from_allele_idx.size();
  int number_of_non_fixed_transmissions = number_of_transmissions - number_of_fixed_transmissions;

  int number_of_canonical_inheritance_vectors = 1 << number_of_non_fixed_transmissions;;
  IntegerVector ibd_states(number_of_canonical_inheritance_vectors);

  // assign allele vector
  IntegerVector x(2 * number_of_ped_members);

  // assign founder alleles
  for (int i_founder = 0; i_founder < ped_row_is_founder_idx.size(); i_founder++){
    int idx_1based = ped_row_is_founder_idx[i_founder];

    x[2 * idx_1based - 2] = 2 * i_founder + 1;
    x[2 * idx_1based - 1] = 2 * i_founder + 2;
  }

  IntegerMatrix ibd_patterns(2 * person_ids.size(), number_of_canonical_inheritance_vectors);

  // drop alleles for each choicce of inheritance vector
  for (int v = 0; v < number_of_canonical_inheritance_vectors; v++){

    // drop alleles down the pedigree for this inheritance vector
    for (int i_order = 0; i_order < top_to_bottom_order.size(); i_order++){

      int i_transmission = top_to_bottom_order[i_order] - 1;

      int from_idx = from_allele_idx[i_transmission] + ((v >> i_transmission) & 1);
      int to_idx = to_allele_idx[i_transmission];

      x[to_idx - 1] = x[from_idx -1];
    }

    // determine multi person ibd patterns
    if (minimal){
      std::unordered_map<int, int> distinct_allele_num_by_placeholder;
      int distinct_allele_num = 1;

      for (int i_person = 0; i_person < person_ids.size(); i_person++){

        // grab allele placeholders for this person
        int idx_person = person_ids[i_person];

        int a = x[2 * (idx_person - 1)];
        int b = x[2 * (idx_person - 1) + 1];

        // do we need to swap labels to obtain the minimal pattern?
        auto it_a = distinct_allele_num_by_placeholder.find(a);
        auto it_b = distinct_allele_num_by_placeholder.find(b);

        bool a_is_new = it_a == distinct_allele_num_by_placeholder.end();
        bool b_is_new = it_b == distinct_allele_num_by_placeholder.end();

        if ((a_is_new && b_is_new) && (a != b)){
          // both new so it may be necessary to swap their labels
          for (int j_person = i_person + 1; j_person < person_ids.size(); j_person++){
            int idx_person2 = person_ids[j_person];

            int c =  x[2 * (idx_person2 - 1)];
            int d =  x[2 * (idx_person2 - 1) + 1];

            // if this is the first genotype with exactly one of a,b
            // then swap if it's b
            bool one_of_a = (a == c) | (a == d);
            bool one_of_b = (b == c) | (b == d);

            if (one_of_a ^ one_of_b){

              if (one_of_b){
                std::swap(a, b);
              }
              break;
            }
          }
        }

        // count distinct alleles
        for (int i_person_allele = 0; i_person_allele < 2; i_person_allele++){
          int allele = i_person_allele == 0 ? a : b;

          // do we have a number for this allele yet?
          auto it = distinct_allele_num_by_placeholder.find(allele);

          if (it == distinct_allele_num_by_placeholder.end()){
            // new distinct allele
            distinct_allele_num_by_placeholder[allele] = distinct_allele_num;
            ibd_patterns(2 * i_person + i_person_allele, v) = distinct_allele_num;

            distinct_allele_num++;
          }
          else{
            // known allele
            ibd_patterns(2 * i_person + i_person_allele, v) = it->second;
          }
        }

        // keep minimal pattern ordered
        if (ibd_patterns(2 * i_person, v) > ibd_patterns(2 * i_person + 1, v)){
          std::swap(ibd_patterns(2 * i_person, v), ibd_patterns(2 * i_person + 1, v));
        }
      }
    }
    else{
      // for debugging, do not attempt to find the minimal pattern
      for (int i_person = 0; i_person < person_ids.size(); i_person++){

        // grab allele placeholders for this person
        int idx_person = person_ids[i_person];

        int a = x[2 * (idx_person - 1)];
        int b = x[2 * (idx_person - 1) + 1];

        ibd_patterns(2 * i_person, v) = a;
        ibd_patterns(2 * i_person + 1, v) = b;
      }
    }


  }

  return ibd_patterns;
}
