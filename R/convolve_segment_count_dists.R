convolve_segment_count_dists <- function(...){
  x <- list(...)
  if (length(x) == 1){
    if (is.list(x[[1]])){
      x <- x[[1]]
    }
  }

  # start with the smallest and work your way up
  scales <- sapply(x, function(d0) length(d0$x))
  o <- order(scales)

  partial_sum <- x[[o[1]]]

  for(i in o[-1]){
    partial_sum <- convolve_two_segment_count_dists(partial_sum,
                                          x[[i]])
  }

  partial_sum
}

convolve_two_segment_count_dists <- function(d1, d2){
  dist_sum <- pmf_of_sum(x1 = d1$x,
                         p1 = d1$px,
                         x2 = d2$x,
                         p2 = d2$px,
                         eps = .Machine$double.xmin)


  class(dist_sum) <- "segment_count_dist"

  dist_sum
}
