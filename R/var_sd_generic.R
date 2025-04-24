#' @export
var <- function(x, ...){
  UseMethod("var")
}

#' @export
var.default <- function(x, ...){
  stats::var(x, ...)
}


#' @export
sd <- function(x, ...){
  UseMethod("sd")
}

#' @export
sd.default <- function(x, ...){
  stats::sd(x, ...)
}
