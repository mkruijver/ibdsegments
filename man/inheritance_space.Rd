% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inheritance_space.R
\name{inheritance_space}
\alias{inheritance_space}
\title{Inheritance space for pedigree}
\usage{
inheritance_space(pedigree, ids, states = "ibd", exploit_symmetries = TRUE)
}
\arguments{
\item{pedigree}{Pedigree in \code{\link[pedtools:ped]{pedtools::ped}} form.}

\item{ids}{Ids for which IBD is observed. Default is \code{pedtools::leaves(pedigree)}.}

\item{states}{One of \code{"ibd"} (default), \code{"kappa"}, \code{"identity"} or \code{"detailed"}.}

\item{exploit_symmetries}{Should symmetries be used to reduce to state space?
This can be set to \code{FALSE} for debugging purposes.}
}
\value{
Object of class \code{inheritance_space}.
}
\description{
The \code{inheritance_space} function determines the space of IBD vectors for a pedigree.
This is mostly for internal use but may be interesting by itself.
}
\examples{
# set up inheritance space for half sib pedigree
i <- inheritance_space(pedigree = pedtools::halfSibPed())

# since there are 2 non-founders, there are 2^4 IBD vectors
# but only 2 distinct states are considered because of symmetries
i

# pry into the internals to see individual pedigree transmissions
i$transmissions
}
