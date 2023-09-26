library(raster)

setwd("/home/brendo/myWork/dsmart/rPackage/dsmart/pkg/data")
list.files()

#dsmart outputs
load("dsmartOutMaps.rda")
dsmartOutMaps

#dsmart covariates
load("dsT_covariates.rda")
dsT_covariates

#dsmart composition
load("dsT_composition.rda")
dsT_composition

#dsmart polygons
load("dsT_polygons.rda")
dsT_polygons

#dsmart polygons
load("dsT_lookup.rda")
dsT_lookup
setwd("/home/brendo/myWork/")


#Run Functions
library(rgdal); library(raster); library(sp); library(gtools); library(C50)
setwd("/home/brendo/myWork/dsmart/data")
#dsmart
dsmart(covariates = dsT_covariates, polygons = dsT_polygons, composition = dsT_composition, n=15, reals = 5, cpus=2)
#dsmartR
dsmartR(rLocs= dsmartOutMaps, nprob = 2, sepP=TRUE, lookup=dsT_lookup , cpus=2)

