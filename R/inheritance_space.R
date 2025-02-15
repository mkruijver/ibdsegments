#' Inheritance space for pedigree
#'
#' The `inheritance_space` function determines the IBD vectors for a pedigree.
#'
#' @export
inheritance_space <- function(pedigree, persons, coefficients = "kappa",
                              exploit_symmetries = TRUE){

  .validate_pedigree(pedigree)

  transmissions <- list_pedigree_transmissions(pedigree,
                                               exploit_symmetries = exploit_symmetries,
                                               symmetries_last = TRUE)

  # determine for each fixed transmission in which other transmissions this founder is involved
  transmissions_fixed_idx <- which(transmissions$is_fixed)

  fixed_transmission_masks <- sapply(transmissions_fixed_idx, function(i)
    to_mask(which(transmissions$from_person_idx == transmissions$from_person_idx[i] &
                    !transmissions$is_fixed) - 1)
  )
  if (length(transmissions_fixed_idx) == 0){
    fixed_transmission_masks <- integer()
  }
  transmissions$masks <- -1
  transmissions$masks[transmissions_fixed_idx] <- fixed_transmission_masks

  # if a transmission is fixed and the mask is 0, then it can be ignored
  transmissions$is_ignorable <- transmissions$masks == 0

  if (any(fixed_transmission_masks == 0)){
    transmissions <- transmissions[c(which(!transmissions$is_ignorable),
                                     which(transmissions$is_ignorable)),]
  }

  i <- list(pedigree = pedigree,
            transmissions = transmissions,
            fixed_transmission_masks = fixed_transmission_masks)

  # determine the IBD status for the persons of interest for each transmission vector
  if (!missing(persons)){
    persons_idx <- match(persons, pedigree$ID)
    if (anyNA(persons_idx)) {
      stop("Person(s) ",
           paste(c( persons[is.na(persons_idx)]), collapse = ", "),
           , " not found in pedigree")
    }

    i$persons <- persons
    i$persons_idx <- persons_idx

    if (length(persons) > 0){
      i$coefficients <- coefficients
      coeff <- .validate_coefficients(coefficients)

      ibd_state_by_v <- get_ibd_states_by_v(number_of_ped_members = length(pedigree$ID),
                                            ped_row_is_founder_idx = which(pedigree$FIDX == 0),
                                            from_allele_idx = transmissions$from_allele_idx,
                                            to_allele_idx = transmissions$to_allele_idx,
                                            persons_idx = persons_idx,
                                            number_of_fixed_transmissions = sum(transmissions$is_fixed),
                                            top_to_bottom_order = transmissions$top_to_bottom_order,
                                            coeff = coeff)

      i$ibd_state_by_v <- ibd_state_by_v
    }
  }

  if (get_option_ignore_irrelevant_transmissions()){
    i$number_of_relevant_transmissions <- sum(!i$transmissions$is_ignorable)
    i$relevant_masks <- i$fixed_transmission_masks[i$fixed_transmission_masks > 0]
  }
  else{
    i$number_of_relevant_transmissions <- nrow(i$transmissions)
    i$relevant_masks <- i$fixed_transmission_masks
  }

  class(i) <- "inheritance_space"
  i
}
