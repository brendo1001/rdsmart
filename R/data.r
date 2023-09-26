#' Composition of soil map units across a small part of the Land Resources of 
#' the Dalrymple Shire map.
#' 
#' A dataset containing the composition of soil map units across a small part of
#' the Land Resources of the Dalrymple Shire map (Rogers \emph{et al}., 1999a, 
#' b) from Queensland, Australia.
#' 
#' @format A data frame with 43 rows and 4 variables: \describe{ 
#'   \item{poly}{polygon id, integer} \item{mapunit}{code of the soil map unit
#'   that the polygon belongs to, character string} \item{soil_class}{code of
#'   the soil class, character string} \item{proportion}{proportion of the soil
#'   map unit that the soil class is assumed to occupy, percent} }
#'   
#' @references Rogers, L.G., Cannon, M.G., Barry, E.V., 1999a. Land Resources of
#'   the Dalrymple Shire, Volume 1. Land Resources Bulletin DNRQ980090. 
#'   Queensland Department of Natural Resources, Brisbane, Queensland. 120 pp.
#'   
#'   Rogers, L.G., Cannon, M.G., Barry, E.V., 1999b. Land Resources of the 
#'   Dalrymple Shire, Volume 2 (Appendices). Land Resources Bulletin DNRQ980090.
#'   Department of Natural Resources, Brisbane, Queensland. 156 pp.
#'   
#'   
"dalrymple_composition"

#' \emph{scorpan} environmental covariates for a small part of the former
#' Dalrymple Shire, Queensland, Australia.
#'
"dalrymple_covariates"

#' Lookup table that relates soil class codes in the soil map unit composition
#' to integer codes in the soil class realisation rasters.
#'
"dalrymple_lookup"

#' Location and soil class of some soil profiles observed across a small part of
#' the former Dalrymple Shire, Queensland, Australia.
#'
"dalrymple_observations"

#' A small part of the Land Resources of the Dalrymple Shire soil map from
#' Queensland, Australia.
#'
"dalrymple_polygons"

#' Several realisations of the potential soil class distribution across a small
#' part of the former Dalrymple Shire, Queensland, Australia.
#'
"dalrymple_realisations"

