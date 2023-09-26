library(raster)
library(rgdal)

polygons <- readOGR(dsn = "D:/LCR/projects/337001-0005 SNutrient/fourpeaks/GIS/soil/nzlri_hurunui_41a_single.shp",
                    layer = "nzlri_hurunui_41a_single")

composition <- read.table("D:/LCR/projects/337001-0005 SNutrient/fourpeaks/hurunui_composition_dsmart_by_geomorphon_updated.csv",
                          header = T, sep = ",")

strata <- raster("D:/LCR/projects/337001-0005 SNutrient/fourpeaks/GIS/geomorphons/hurunui_25m_geomorphons_or15.tif")

# Load covariates
covr <- stack(list.files(path = "D:/LCR/projects/337001-0005 SNutrient/fourpeaks/GIS/dsmart_covariates/",
                         pattern = ".tif$", full.names = TRUE))

outputdir <- "D:/LCR/projects/337001-0005 SNutrient/fourpeaks/"

# samps <- stratifiedVirtualSamples(covariates, polygons, composition, strata,
#                                   n.realisations = 10, rate = 1, 
#                                   method.sample = "by_area", 
#                                   method.allocate = "weighted")

start <- Sys.time()
dsmart(covr, polygons, composition, rate = 1, reals = 5, 
       observations = NULL, method.sample = "by_area", 
       method.allocate = "weighted", strata = strata, nprob = 3, 
       outputdir = outputdir, stub = NULL, cpus = 6)
finish <- Sys.time()
print(finish - start)