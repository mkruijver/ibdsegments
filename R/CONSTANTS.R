COEFF_BY_NAME <- setNames(c(1L, 2L, 9L, 15L, 99L),
                          nm = c("ibd", "kappa", "identity", "detailed", "v"))

NAME_BY_COEFF_NAME <- setNames(names(COEFF_BY_NAME), COEFF_BY_NAME)

COEFF_IBD <- COEFF_BY_NAME["ibd"]
COEFF_KAPPA <- COEFF_BY_NAME["kappa"]
COEFF_IDENTITY <- COEFF_BY_NAME["identity"]
COEFF_DETAILED <- COEFF_BY_NAME["detailed"]
COEFF_V <- COEFF_BY_NAME["v"]

VALID_OBS_BY_COEFF_NAME <- list(ibd = 0:2,
                                kappa = 0:2,
                                identity = 1:9,
                                detailed = 1:15,
                                v = "unchecked")
