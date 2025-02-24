
v_prior_with_f <- function(pedigree, i){

  if (is.null(i$relevant_masks_from_id_idx)){
    stop("inheritance_space was not constructed with exploit_symmetries=TRUE")
  }

  founder_id_idx_with_multiple_offspring <- i$relevant_masks_from_id_idx
  founder_f <- pedtools::founderInbreeding(pedigree, ids=founder_id_idx_with_multiple_offspring)

  v_prior_with_f_cpp(founder_masks = i$relevant_masks_from_id_idx,
                     founder_f = founder_f,
                     number_of_transmissions = i$number_of_relevant_transmissions,
                     number_of_fixed_transmissions = length(i$relevant_masks))
}
