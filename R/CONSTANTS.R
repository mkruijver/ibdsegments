STATES_BY_NAME <- stats::setNames(c(1L, 2L, 9L, 15L, 99L),
                          nm = c("ibd", "kappa", "identity", "detailed", "v"))

NAME_BY_STATES_NAME <- stats::setNames(names(STATES_BY_NAME), STATES_BY_NAME)

STATES_IBD <- STATES_BY_NAME["ibd"]
STATES_KAPPA <- STATES_BY_NAME["kappa"]
STATES_IDENTITY <- STATES_BY_NAME["identity"]
STATES_DETAILED <- STATES_BY_NAME["detailed"]
STATES_V <- STATES_BY_NAME["v"]

VALID_OBS_BY_STATES_NAME <- list(ibd = 0:2,
                                kappa = 0:2,
                                identity = 1:9,
                                detailed = 1:15,
                                v = "unchecked")
