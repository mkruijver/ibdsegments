.validate_states <- function(states){
  states_value <- match.arg(states,
                                choices = names(STATES_BY_NAME))

  if (states_value %in% names(STATES_BY_NAME)){
    return(STATES_BY_NAME[states_value])
  } else{
    stop("Unknown states value:", states_value)
  }
}

.check_ids_compatible_with_states_idx <- function(ids, states_idx){

  if (states_idx == STATES_IBD){
    if (length(ids) < 1){
      stop("ids needs to have at least length one for IBD states")
    }
  }
  else{
    if (length(ids) != 2){
      stop("ids need to have length 2 for the chosen states ",
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

.validate_obs_compatible_with_states_idx <- function(obs, argument_name, states_idx){
  nm <- NAME_BY_STATES_NAME[as.character(states_idx)]
  if (is.null(nm)) stop("Unknown states value: ", states_idx)

  valid_obs <- VALID_OBS_BY_STATES_NAME[[nm]]
  if (is.null(valid_obs)) stop("Could not obtain valid values for states value: ", states_idx)

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
           paste(valid_obs, collapse = ", "), " for ", nm, " states. ",
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

.validate_ibs_dist <- function(x){
  if (!inherits(x, "ibd_dist")){
    stop("distribution should be of class ibd_dist")
  }
}

.validate_not_ibd_fraction_dist <- function(x){
  if (x$fraction){
    stop("IBD fraction distributions are not supported in convolution")
  }
}

.validate_integer <- function(x, argument_name){
  if (!is.integer(x)){
    if (is.numeric(x)){
      x_as_integer <- as.integer(x)

      if (isTRUE(all.equal(x, x_as_integer))){
        return()
      }
    }
  }
  else{
    return()
  }

  stop(argument_name, " should be an integer vector")
}

.validate_chromosome_length <- function(chromosome_length){
  if (!is.numeric(chromosome_length)){
    stop("chromosome_length needs to be numeric")
  }

  if (any(chromosome_length <= 0)){
    stop("chromosome_length needs to be strictly positive")
  }
}

.validate_numeric <- function(x, argument_name){
  if (!is.numeric(x)){
    stop(argument_name, " should be a numeric vector")
  }
}
.validate_rcibd_result <- function(x){

  if (is.null(x$samples)){
    stop("x$samples is NULL")
  }

  .validate_integer(x$samples$sample, "x$samples$sample")
  .validate_integer(x$samples$chromosome, "x$samples$chromosome")
  .validate_numeric(x$samples$length, "x$samples$length")
  .validate_integer(x$samples$state, "x$samples$state")
}
