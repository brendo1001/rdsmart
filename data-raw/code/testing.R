library(raster)
library(rgdal)
library(sp)

# Load covariates to raster stack
covariates <- stack()
for(filename in list.files(path = "C:/Users/OdgersN/Dropbox/datasets/burdekin/test25km/covariates/",
                           pattern = ".img$", full.names = TRUE))
{
  r <- raster(filename)
  
  # Load raster to stack
  covariates <- stack(covariates, r)
}
plot(covariates[[1]])

shp <- readOGR(dsn="C:/Users/OdgersN/Dropbox/datasets/burdekin/test25km/polygons/test25km_polys.shp",
               layer="test25km_polys")
plot(shp)

composition = read.table("C:/Users/OdgersN/Dropbox/datasets/burdekin/test25km/polygons/formatted_attribute_table.txt",
                         sep=",", header=FALSE)
colnames(composition) = c("poly","mapunit","soil_class","proportion")


# start <- Sys.time()
# disagg(covariates, shp, composition, n = 15, reals = 10, cpus = 6)
# finish <- Sys.time()
# print(finish - start)

# On the 25-km Burdekin test area:
# 10 realisations took 11.86173 mins
# 20 realisations took 18.01235 mins
# 40 realisations took 30.06257 mins
# 80 realisations took 57.39468 mins
# 100 realisations took 70.51134 mins

# # Load realisations to stack
# realisations <- stack()
# for(filename in list.files(path = "D:/R/dsmart/rPackage/dsmart/pkg/output/realisations/",
#                            pattern = ".tif$", full.names = TRUE))
# {
#   r <- raster(filename)
#   
#   # Load raster to stack
#   realisations <- stack(realisations, r)
# }
# 
# # Read lookup table
# lookup <- read.table("D:/R/dsmart/rPackage/dsmart/pkg/classLookupTable.txt", header = TRUE, sep = ",")
# 
# start <- Sys.time()
# dsmartR(realisations, lookup, nprob = 5, cpus = 6)
# finish <- Sys.time()
# print(finish - start)
# 
# 
# 
# 
# observations <- read.table("C:/Users/OdgersN/Dropbox/datasets/burdekin/test25km/obs/burdekin_25km_obs_3col.txt", 
#                            header = TRUE, sep = ",")
# 
#   
# s.test <- getSamples(covariates, shp, composition, n.realisations = 5, n.samples = 15, method = "weighted")

#####
start <- Sys.time()
dsmart(covariates, shp, composition, rate = 15, method.sample = "by_polygon", reals = 5, cpus = 6, 
       observations = NULL, nprob = 3, outputdir = "D:/burdekin/test25km",
       stub = "test_25km")
finish <- Sys.time()
print(finish - start)
#####


writeOGR(burdekin_polygons, dsn = "D:/burdekin/7kmpolys.shp", "7kmpolys", driver = "ESRI Shapefile")

burdekin_observations <- read.table("D:/burdekin/7kmobservations.csv", header = TRUE, sep = ",")
burdekin_observations <- burdekin_observations[, c(1, 2, 9)]
colnames(burdekin_observations) <- c("x", "y", "soil_class")
devtools::use_data(burdekin_observations)

data(burdekin_composition)
dalrymple_composition <- burdekin_composition
devtools::use_data(dalrymple_composition)

data(burdekin_covariates)
dalrymple_covariates <- burdekin_covariates
devtools::use_data(dalrymple_covariates)

data(burdekin_polygons)
dalrymple_polygons <- burdekin_polygons
devtools::use_data(dalrymple_polygons)

dalrymple_observations <- burdekin_observations
devtools::use_data(dalrymple_observations)

dalrymple_lookup <- lookup
devtools::use_data(dalrymple_lookup)

dalrymple_realisations <- realisations
devtools::use_data(dalrymple_realisations)
