.validate_coefficients <- function(coefficients){
  coefficient_type <- match.arg(coefficients,
                                choices = names(COEFF_BY_NAME))

  if (coefficient_type %in% names(COEFF_BY_NAME)){
    return(COEFF_BY_NAME[coefficient_type])
  } else{
    stop("Unknown coefficient type:", coefficient_type)
  }
}

.check_persons_compatible_with_coeff <- function(persons, coeff){

  if (coeff == COEFF_IBD){
    if (length(persons) < 2){
      stop("Persons needs to have at least length two for IBD status coefficients")
    }
  }
  else{
    if (length(persons) != 2){
      stop("Persons need to have length 2 for the chosen coefficients. ",
           "Actual length: ", length(persons))
    }
  }
}

.validate_obs_compatible_with_coeff <- function(obs, argument_name, coeff){
  nm <- NAME_BY_COEFF_NAME[as.character(coeff)]
  if (is.null(nm)) stop("Unknown coefficient value: ", coeff)

  valid_obs <- VALID_OBS_BY_COEFF_NAME[[nm]]
  if (is.null(valid_obs)) stop("Could not obtain valid values for coefficient: ", coeff)

  obs_is_invalid <- is.na(match(obs, table = valid_obs))

  if (any(obs_is_invalid)){
    stop(argument_name, " should only take values ",
         paste(valid_obs, collapse = ", "), " for ", nm, " coefficients. ",
         "Invalid value(s): ", paste0(head(obs[obs_is_invalid]), collapse = ", "))
  }
}

.validate_pedigree <- function(pedigree){
  if (!inherits(pedigree, "ped")){
    stop("pedigree should be of class ped")
  }

  if (!pedtools::hasParentsBeforeChildren(pedigree)){
    stop("pedigree should list parents before children")
  }
}

.validate_pedigree_ids <- function(ids, pedigree){

  person_idx <- match(ids, table = pedigree$ID)

  if (anyNA(person_idx)) {
    stop("Id(s) ",
         paste(c( ids[is.na(person_idx)]), collapse = ", "),
         , " not found in pedigree")
  }

  person_idx
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
