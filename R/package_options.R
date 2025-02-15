.packageOptions <- new.env(parent = emptyenv())
.packageOptions$ignore_irrelevant_transmissions <- TRUE

set_option_ignore_irrelevant_transmissions <- function(value) {
  .packageOptions$ignore_irrelevant_transmissions <- value
}

get_option_ignore_irrelevant_transmissions <- function() {
  .packageOptions$ignore_irrelevant_transmissions
}
