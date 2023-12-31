% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.r
\name{.allocate}
\alias{.allocate}
\title{Allocate samples to soil classes}
\usage{
.allocate(classes, n = 15, method = "random", weights = NULL)
}
\arguments{
\item{classes}{The classes to allocate to.}

\item{n}{The number of samples to allocate.}

\item{method}{The method of allocation. Valid values are \code{"weighted"},
for weighted-random allocation using the weights in \code{weights}, and
\code{"random"} for random allocation (the default).}
}
\description{
\code{.allocate} randomly or weighted-randomly draws a sample of size
\code{n} from \code{classes}.
}
