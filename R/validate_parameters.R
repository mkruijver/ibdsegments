.validate_coefficients <- function(coefficients){
  coefficient_type <- match.arg(coefficients,
                                choices = names(COEFF_BY_NAME))

  if (coefficient_type %in% names(COEFF_BY_NAME)){
    return(COEFF_BY_NAME[coefficient_type])
  } else{
    stop("Unknown coefficient type:", coefficient_type)
  }
}

.check_ids_compatible_with_coeff <- function(ids, coeff){

  if (coeff == COEFF_IBD){
    if (length(ids) < 1){
      stop("ids needs to have at least length one for IBD status coefficients")
    }
  }
  else{
    if (length(ids) != 2){
      stop("ids need to have length 2 for the chosen coefficients. ",
           "Actual length: ", length(ids))
    }
  }
}

.assert_no_founder_inbreeding <- function(pedigree, reason){
  if (pedtools::hasInbredFounders(pedigree)){
    stop(reason)
  }
}

.assert_ids_are_not_inbred_founders <- function(pedigree, ids){

  ids_idx <- match(ids, pedigree$ID)
  ids_idx_is_founder <- pedigree$FIDX[ids_idx] == 0L

  f <- pedtools::founderInbreeding(pedigree, ids[ids_idx_is_founder],)
  if (any(f > 0)){
    stop("ids should not include inbred founders")
  }
}

.validate_obs_compatible_with_coeff <- function(obs, argument_name, coeff){
  nm <- NAME_BY_COEFF_NAME[as.character(coeff)]
  if (is.null(nm)) stop("Unknown coefficient value: ", coeff)

  valid_obs <- VALID_OBS_BY_COEFF_NAME[[nm]]
  if (is.null(valid_obs)) stop("Could not obtain valid values for coefficient: ", coeff)

  is_unchecked <- FALSE
  if (length(valid_obs) == 1){
    if (isTRUE(valid_obs[1] == "unchecked")){
      is_unchecked <- TRUE
    }
  }

  if (!is_unchecked){
    obs_is_invalid <- is.na(match(obs, table = valid_obs))

    if (any(obs_is_invalid)){
      stop(argument_name, " should only take values ",
           paste(valid_obs, collapse = ", "), " for ", nm, " coefficients. ",
           "Invalid value(s): ", paste0(head(obs[obs_is_invalid]), collapse = ", "))
    }
  }
}

.validate_pedigree <- function(pedigree, continuous_genome = FALSE){
  if (!inherits(pedigree, "ped")){
    stop("pedigree should be of class ped")
  }

  if (!pedtools::hasParentsBeforeChildren(pedigree)){
    stop("pedigree should list parents before children")
  }

  if (continuous_genome){
    .assert_no_founder_inbreeding(pedigree, "founder inbreeding is not supported for continuous ibd methods")
  }
}

.validate_pedigree_ids <- function(ids, pedigree){

  ids_idx <- match(ids, table = pedigree$ID)

  if (anyNA(ids_idx)) {
    stop("Id(s) ",
         paste(c( ids[is.na(ids_idx)]), collapse = ", "),
         , " not found in pedigree")
  }

  ids_idx
}

.validate_recombination_rates_compatible_with_obs <- function(obs, argument_name, recombination_rates){
  expected_length_of_recombination_rates <- length(obs) - 1L

  if (length(recombination_rates) != expected_length_of_recombination_rates){
    stop(argument_name, " has length ", length(recombination_rates),
         "; ", expected_length_of_recombination_rates, " was expected")
  }
}

.validate_logical <- function(x, argument_name){
  if (length(x) != 1L || !is.logical(x)){
    stop(argument_name, " should be TRUE or FALSE")
  }
}
