to_mask <- function(set_bits) {
  sum(bitwShiftL(1, set_bits))
}
