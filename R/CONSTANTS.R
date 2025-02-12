COEFF_BY_NAME <- setNames(c(2L, 9L, 15L), nm = c("kappa", "identity", "detailed"))

NAME_BY_COEFF <- setNames(names(COEFF_BY_NAME), COEFF_BY_NAME)

COEFF_KAPPA <- COEFF_BY_NAME["kappa"]
COEFF_IDENTITY <- COEFF_BY_NAME["identity"]
COEFF_DETAILED <- COEFF_BY_NAME["detailed"]

VALID_OBS_BY_COEFF_NAME <- list(kappa = 0:2,
                                identity = 1:9,
                                detailed = 1:15)
