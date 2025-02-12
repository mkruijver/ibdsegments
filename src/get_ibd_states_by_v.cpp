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
std::vector<int> minimal_pattern(IntegerVector x, IntegerVector person_idx_1_based){

  std::vector<int> ibd_pattern(person_idx_1_based.size()*2);

  std::unordered_map<int, int> distinct_allele_num_by_placeholder;
  int distinct_allele_num = 1;

  for (int i_person = 0; i_person < person_idx_1_based.size(); i_person++){

    // grab allele placeholders for this person
    int idx_person = person_idx_1_based[i_person];

    int a = x[2 * (idx_person - 1)];
    int b = x[2 * (idx_person - 1) + 1];

    // do we need to swap labels to obtain the minimal pattern?
    auto it_a = distinct_allele_num_by_placeholder.find(a);
    auto it_b = distinct_allele_num_by_placeholder.find(b);

    bool a_is_new = it_a == distinct_allele_num_by_placeholder.end();
    bool b_is_new = it_b == distinct_allele_num_by_placeholder.end();

    if ((a_is_new && b_is_new) && (a != b)){
      // both new so it may be necessary to swap their labels
      for (int j_person = i_person + 1; j_person < person_idx_1_based.size(); j_person++){
        int idx_person2 = person_idx_1_based[j_person];

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
        ibd_pattern[2 * i_person + i_person_allele] = distinct_allele_num;

        distinct_allele_num++;
      }
      else{
        // known allele
        ibd_pattern[2 * i_person + i_person_allele] = it->second;
      }
    }

    // keep minimal pattern ordered
    if (ibd_pattern[2 * i_person] > ibd_pattern[2 * i_person + 1]){
      std::swap(ibd_pattern[2 * i_person], ibd_pattern[2 * i_person + 1]);
    }
  }

  return ibd_pattern;
}


// [[Rcpp::export]]
List get_multi_ibd_patterns_by_v(int number_of_ped_members,
                                  IntegerVector ped_row_is_founder_idx,
                                  IntegerVector from_allele_idx,
                                  IntegerVector to_allele_idx,
                                  IntegerVector person_ids,
                                  int number_of_fixed_transmissions,
                                  IntegerVector top_to_bottom_order,
                                  bool minimal = true){

  int pattern_size = person_ids.size() * 2;

  int number_of_transmissions = from_allele_idx.size();
  int number_of_non_fixed_transmissions = number_of_transmissions - number_of_fixed_transmissions;
  int number_of_canonical_inheritance_vectors = 1 << number_of_non_fixed_transmissions;

  std::set<std::vector<int>> unique_patterns_set;
  std::vector<const std::vector<int>*> pattern_references_by_v(number_of_canonical_inheritance_vectors);

  // assign allele vector
  IntegerVector x = assign_founder_alleles(number_of_ped_members, ped_row_is_founder_idx);

  // drop alleles for each choice of inheritance vector
  for (int v = 0; v < number_of_canonical_inheritance_vectors; v++){

    // drop alleles down the pedigree for this inheritance vector
    drop_founder_alleles(x, v, from_allele_idx, to_allele_idx, top_to_bottom_order);

    if (minimal){
      // determine multi person ibd patterns
      std::vector<int> m = minimal_pattern(x, person_ids);

      auto insertion_result = unique_patterns_set.insert(m);
      pattern_references_by_v[v] = &(*insertion_result.first);
    }
    else{
      std::vector<int> m(pattern_size);

      for (int i_person = 0; i_person < person_ids.size(); i_person++){

        // grab allele placeholders for this person
        int idx_person = person_ids[i_person];

        int a = x[2 * (idx_person - 1)];
        int b = x[2 * (idx_person - 1) + 1];

        m[2 * i_person] = a;
        m[2 * i_person + 1] = b;
      }

      auto insertion_result = unique_patterns_set.insert(m);
      pattern_references_by_v[v] = &(*insertion_result.first);
    }
  }

  // list out all unique patterns
  int number_of_unique_patterns = unique_patterns_set.size();
  IntegerMatrix unique_patterns(pattern_size, number_of_unique_patterns);

  auto dest = unique_patterns.begin();
  for (const auto& element : unique_patterns_set) {
    std::copy(element.begin(), element.end(), dest);
    dest += pattern_size;
  }

  // list out index among unique patterns by v
  IntegerVector pattern_idx_by_v(number_of_canonical_inheritance_vectors);

  for (int v = 0; v < number_of_canonical_inheritance_vectors; v++) {

    const std::vector<int>* ref = pattern_references_by_v[v];
    auto it = unique_patterns_set.find(*ref);

    int index = std::distance(unique_patterns_set.begin(), it);
    pattern_idx_by_v[v] = index + 1;
  }

  return List::create( _["unique_patterns"] = unique_patterns,
                       _["pattern_idx_by_v"] = pattern_idx_by_v);
}

// [[Rcpp::export]]
DataFrame multi_ibd_patterns_by_v_df(IntegerMatrix unique_patterns,
                                     IntegerVector pattern_idx_by_v,
                                     CharacterVector ids,
                                     double pr_v){

  int number_of_ids = ids.size();
  int number_of_patterns = unique_patterns.ncol();

  if (unique_patterns.nrow() != 2 * number_of_ids){
    Rcpp::stop("unique_patterns and ids have incompatible dimensions");
  }

  // start building list to be converted to df
  Rcpp::List result(1 + number_of_ids);
  Rcpp::CharacterVector col_names(result.size());

  // probability column
  col_names[0] = "Prob";
  std::vector<long> count_by_pattern(number_of_patterns);
  for (const auto pattern_idx: pattern_idx_by_v){
    count_by_pattern[pattern_idx-1]++;
  }

  Rcpp::NumericVector pr(number_of_patterns);
  for (int i = 0; i < number_of_patterns; i++){
    pr[i] = count_by_pattern[i] * pr_v;
  }
  result[0] = pr;

  // id columns
  for (int i = 0; i < number_of_ids; i++){
    CharacterVector id_column(number_of_patterns);

    for (int i_pattern = 0; i_pattern < number_of_patterns; i_pattern++){
      int a = unique_patterns(2*i, i_pattern);
      int b = unique_patterns(2*i + 1, i_pattern);

      id_column[i_pattern] = std::to_string(a) + " " + std::to_string(b);
    }

    result[1+i] = id_column;
    col_names[1+i] = ids[i];
  }

  DataFrame df = Rcpp::DataFrame(result);
  df.attr("names") = col_names;

  return df;
}

std::vector<int> which(LogicalVector x){
  std::vector<int> these;
  these.reserve(x.size());

  for (int i = 0; i < x.size(); i++){
    if (x[i]){
      these.push_back(i);
    }
  }

  return these;
}

// [[Rcpp::export]]
DataFrame multi_ibd_patterns_df(NumericVector prob,
                                IntegerMatrix multi_locus_m_idx,
                                IntegerMatrix unique_patterns,
                                CharacterVector ids){

  int number_of_ids = ids.size();
  int number_of_patterns = unique_patterns.ncol();

  if (unique_patterns.nrow() != 2 * number_of_ids){
    Rcpp::stop("unique_patterns and ids have incompatible dimensions");
  }

  // only return rows with positive probability
  LogicalVector prob_is_positive = prob > 0;
  int number_of_rows = Rcpp::sum(prob_is_positive);
  std::vector<int> idx_positive = which(prob_is_positive);

  // start building list to be converted to df
  int number_of_loci = multi_locus_m_idx.ncol();
  Rcpp::List result(1 + number_of_ids * number_of_loci);
  Rcpp::CharacterVector col_names(result.size());

  // probability column
  col_names[0] = "Prob";
  result[0] = prob[prob_is_positive];

  // id columns by locus
  int i_df_col = 0;
  for (int i_locus = 0; i_locus < number_of_loci; i_locus++){
    for (int i_person = 0; i_person < number_of_ids; i_person++){
      CharacterVector id_column(number_of_rows);

      for (int i_row = 0; i_row < number_of_rows; i_row++){
        int idx = idx_positive[i_row];

        int i_pattern = multi_locus_m_idx(idx, i_locus) - 1;

        int a = unique_patterns(2 * i_person, i_pattern);
        int b = unique_patterns(2 * i_person + 1, i_pattern);

        id_column[i_row] = std::to_string(a) + " " + std::to_string(b);
      }

      result[1 + i_df_col] = id_column;

      std::string id = Rcpp::as<std::string>(ids[i_person]);
      col_names[1 + i_df_col] =  id + "@"+ std::to_string(i_locus+1) +"";

      i_df_col++;
    }
  }

  DataFrame df = Rcpp::DataFrame(result);
  df.attr("names") = col_names;

  return df;
}
