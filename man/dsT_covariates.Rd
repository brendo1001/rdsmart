\name{dsT_covariates}
\docType{data}
\alias{dsT_covariates}
\title{Environmental covariate suite of rasters for testing DSMART}
\description{A library of 20 environmental covariate rasters (GeoTIFF format). These data are described in Odgers et al. (2014). The assembled covariates are surrogates for soil forming factors and are commonly referred to as scorpan factors, where scorpan is a mnemonic of the soil forming factors (s) soil, (c) climate, (o) organism, (r) relief, (p) parent material, (a) age, and (n) spatial position (McBratney et al. 2003). The raster represent a small subset from the Dalrymple Shire in Central Queensland, Australia. The subset area is a 7km by 7km square, and the rasters have a 30m pixel resolution. The coordinate reference system of these maps is GDA 1994 Australia Albers.}

\usage{terra::rast(list.files(path = system.file("extdata/gtiff_covs/",package="rdsmart"),pattern = ".tif",full.names = TRUE))}
\format{
\code{dsT_covariates} are 20 layer small raster objects.}
\details{This is an example data set in order to familiarize users to the working and inputs required for operationalising DSMART which is implemented in the companion functions \code{disaggregate} and \code{summarise}. }

\references{
\itemize{
\item McBratney, A.B., Mendonca Santos, M.L., Minasny, B., (2003) \href{http://dx.doi.org/10.1016/S0016-7061(03)00223-4}{On digital soil mapping}. Geoderma, 117: 3-52.
\item Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., (2014) \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{Disaggregating and harmonising soil map units through resampled classification trees}. Geoderma, 214-215: 91-100.
}
}

\examples{
library(rdsmart)
library(terra)
terra::rast(list.files(path = system.file("extdata/gtiff_covs/",package="rdsmart"),pattern = ".tif",full.names = TRUE))


}
\keyword{datasets}
