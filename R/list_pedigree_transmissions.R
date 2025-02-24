list_pedigree_transmissions <- function(pedigree,
                                        exploit_symmetries = TRUE,
                                        symmetries_last = FALSE){

  ped_row_is_founder <- pedigree$FIDX == 0 & pedigree$MIDX == 0
  ped_row_is_non_founder <- !ped_row_is_founder

  number_of_non_founders <- sum(ped_row_is_non_founder)

  ped_row_is_founder_idx <- which(ped_row_is_founder)

  ped_non_founder_row_idx <- which(ped_row_is_non_founder)
  to_id_idx <-  rep(ped_non_founder_row_idx, each = 2)

  ped_non_founder_fidx <- pedigree$FIDX[ped_non_founder_row_idx]
  ped_non_founder_midx <- pedigree$MIDX[ped_non_founder_row_idx]

  allele <- rep(c(1L, 2L), number_of_non_founders)

  father_idx <- rep(ped_non_founder_fidx, each = 2)
  mother_idx <- rep(ped_non_founder_midx, each = 2)

  from_id_idx <- ifelse(allele == 1L, father_idx, mother_idx)

  from_allele_idx <- 2 * from_id_idx - 1

  to_allele_idx <- 2 * to_id_idx - 2 + allele

  # founder symmetry: fix first transmissions for founders
  is_fixed <- rep(FALSE, length(father_idx))

  if (exploit_symmetries){
    is_fixed[ match(ped_row_is_founder_idx, table = from_id_idx)] <- TRUE
  }

  transmissions <- data.frame(to_id_idx = to_id_idx,
                              from_id_idx = from_id_idx,
                              father_idx = father_idx,
                              mother_idx = mother_idx,
                              allele = allele,
                              from_allele_idx = from_allele_idx,
                              to_allele_idx = to_allele_idx,
                              is_fixed = is_fixed,
                              is_from_male = allele == 1L)

  if (symmetries_last){
    # put the fixed ones at the end
    transmissions_sorted <- rbind(transmissions[!transmissions$is_fixed,],
                                  transmissions[transmissions$is_fixed,])

    transmissions_sorted$top_to_bottom_order <- order(as.integer(rownames(transmissions_sorted)))

    return(transmissions_sorted)
  }
  else{
    return(transmissions)
  }
}
