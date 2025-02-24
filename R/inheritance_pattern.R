#' Inheritance pattern for inheritance vectors
#'
#' The `inheritance_pattern` function determines the inheritance pattern for all
#' pedigree members by dropping the founder allele labels down the pedigree
#' according to the IBD vector `v`.
#'
#' @examples
#' ped_fs <- pedtools::nuclearPed(nch = 2)
#' i <- inheritance_space(ped_fs, ids = 3:4)
#'
#' # show the inheritance pattern and IBD state for all canonical IBD vectors
#' inheritance_pattern(i, v = 0:3)
#'
#' # wihtout exploiting founder symmetry
#' i2 <- inheritance_space(ped_fs, ids = 3:4, exploit_symmetries = FALSE)
#' inheritance_pattern(i2, v = 0:15)
#' @export
inheritance_pattern <- function(inheritance_space, v){

  if (length(v) > 1){
    return(do.call(rbind, lapply(v, function(v0) inheritance_pattern(inheritance_space, v0))))
  }

  founder_labels <- get_founder_labels_for_v(v = v,
                                             number_of_ped_members = length(inheritance_space$pedigree$ID),
                                             ped_row_is_founder_idx = which(inheritance_space$pedigree$FIDX == 0),
                                             from_allele_idx = inheritance_space$transmissions$from_allele_idx,
                                             to_allele_idx = inheritance_space$transmissions$to_allele_idx,
                                             number_of_fixed_transmissions = sum(inheritance_space$transmissions$is_fixed),
                                             top_to_bottom_order = inheritance_space$transmissions$top_to_bottom_order)

  df <- data.frame(v = v,
                   t(setNames(founder_labels,
                              nm = inheritance_space$pedigree$ID)),
                   check.names = FALSE, row.names = NULL)

  if (!is.null(inheritance_space$ibd_state_by_v)){
    df$state <- inheritance_space$ibd_state_by_v[v + 1]
  }

  df
}
