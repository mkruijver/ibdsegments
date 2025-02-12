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

  if (coeff == COEFF_KAPPA){
    if (length(persons) < 2){
      stop("Persons needs to have at least length two for kappa coefficients")
    }
  }
  else{
    stop("Persons need to have length two for the chosen coefficients")
  }
}

.validate_obs_compatible_with_coeff <- function(obs, argument_name, coeff){
  nm <- NAME_BY_COEFF[coeff]
  if (is.null(nm)) stop("Unknown coefficient value: ", coeff)

  valid_obs <- VALID_OBS_BY_COEFF_NAME[[nm]]
  if (is.null(valid_obs)) stop("Could not obtain valid values for coefficient: ", coeff)

  if (anyNA(match(obs, table = valid_obs))){
    stop(argument_nm, " should only take values ",
         paste(valid_obs, collapse = ", "), " for ", nm, " coefficients")
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
