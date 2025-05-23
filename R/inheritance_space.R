#' Inheritance space for pedigree
#'
#' The `inheritance_space` function determines the space of IBD vectors for a pedigree.
#' This is mostly for internal use but may be interesting by itself.
#'
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param ids Ids for which IBD is observed. Default is `pedtools::leaves(pedigree)`.
#' @param states One of `"ibd"` (default), `"kappa"`, `"identity"` or `"detailed"`.
#' @param exploit_symmetries Should symmetries be used to reduce to state space?
#'                           This can be set to `FALSE` for debugging purposes.
#'
#' @returns Object of class `inheritance_space`.
#' @examples
#' # set up inheritance space for half sib pedigree
#' i <- inheritance_space(pedigree = pedtools::halfSibPed())
#'
#' # since there are 2 non-founders, there are 2^4 IBD vectors
#' # but only 2 distinct states are considered because of symmetries
#' i
#'
#' # pry into the internals to see individual pedigree transmissions
#' i$transmissions
#' @export
inheritance_space <- function(pedigree, ids, states = "ibd",
                              exploit_symmetries = TRUE){

  .validate_pedigree(pedigree)

  transmissions <- list_pedigree_transmissions(pedigree,
                                               exploit_symmetries = exploit_symmetries,
                                               symmetries_last = TRUE)

  # determine for each fixed transmission in which other transmissions this founder is involved
  transmissions_fixed_idx <- which(transmissions$is_fixed)

  fixed_transmission_masks <- sapply(transmissions_fixed_idx, function(i)
    to_mask(which(transmissions$from_id_idx == transmissions$from_id_idx[i] &
                    !transmissions$is_fixed) - 1)
  )
  if (length(transmissions_fixed_idx) == 0){
    fixed_transmission_masks <- integer()
  }
  fixed_transmission_masks_from_id_idx <- transmissions$from_id_idx[transmissions_fixed_idx]

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
            symmetries = exploit_symmetries,
            fixed_transmission_masks = fixed_transmission_masks)

  # determine the IBD status for the ids of interest for each transmission vector
  if (!missing(ids)){
    ids_idx <- match(ids, pedigree$ID)
    if (anyNA(ids_idx)) {
      stop("Id(s) ",
           paste(c( ids[is.na(ids_idx)]), collapse = ", "),
           , " not found in pedigree")
    }

    i$ids <- ids
    i$ids_idx <- ids_idx

    if (length(ids) > 0){
      i$states <- states
      states_idx <- .validate_states(states)

      ibd_state_by_v <- get_ibd_states_by_v(number_of_ped_members = length(pedigree$ID),
                                            ped_row_is_founder_idx = which(pedigree$FIDX == 0),
                                            from_allele_idx = transmissions$from_allele_idx,
                                            to_allele_idx = transmissions$to_allele_idx,
                                            ids_idx = ids_idx,
                                            number_of_fixed_transmissions = sum(transmissions$is_fixed),
                                            top_to_bottom_order = transmissions$top_to_bottom_order,
                                            states_idx = states_idx)

      i$ibd_state_by_v <- ibd_state_by_v
    }
  }

  if (get_option_ignore_irrelevant_transmissions()){
    i$number_of_relevant_transmissions <- sum(!i$transmissions$is_ignorable)
    i$relevant_masks <- i$fixed_transmission_masks[i$fixed_transmission_masks > 0]
    i$relevant_masks_from_id_idx <- fixed_transmission_masks_from_id_idx[i$fixed_transmission_masks > 0]

  }
  else{
    i$number_of_relevant_transmissions <- nrow(i$transmissions)
    i$relevant_masks <- i$fixed_transmission_masks
    i$relevant_masks_from_id_idx <- fixed_transmission_masks_from_id_idx
  }

  class(i) <- "inheritance_space"
  i
}


#' @export
print.inheritance_space <- function(x, ...){
  number_of_pedigree_members <- length(x$pedigree$ID)
  number_of_pedigree_founders <- sum(x$pedigree$FIDX==0)
  number_of_pedigree_non_founders <- number_of_pedigree_members - number_of_pedigree_founders

  cat("Inheritance space for pedigree of size", number_of_pedigree_members,
      "with",number_of_pedigree_founders, "founder(s) and",
      number_of_pedigree_non_founders, "non-founder(s)\n")

  number_of_ibd_vectors <- 4^number_of_pedigree_non_founders

  cat("# IBD vectors including symmetries:", number_of_ibd_vectors, "\n")

  number_of_transmissions <- nrow(x$transmissions)
  number_of_fixed_transmissions <- length(x$fixed_transmission_masks)
  number_of_active_transmissions <- number_of_transmissions - number_of_fixed_transmissions
  number_of_ibd_vectors_considered <- 2^number_of_active_transmissions

  cat("# canonical IBD vectors:", number_of_ibd_vectors_considered, "\n")
}


