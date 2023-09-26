#data
#covariates
load(file = "C:/rdev/dsmart/rPackage/dsmart/pkg/data/dalrymple_covariates.rda")

#polygons
load(file = "C:/rdev/dsmart/rPackage/dsmart/pkg/data/dalrymple_polygons.rda")

#compostions
load(file = "C:/rdev/dsmart/rPackage/dsmart/pkg/data/dalrymple_composition.rda")


setwd("C:/temp")

library(C50)
library(ranger)
library(raster)
library(caret)

myargs<-list(num.trees = 100, trControl = trainControl(method = "oob"))




## Disagg with no factors
disag_out<-disaggregate(covariates = dalrymple_covariates, polygons = dalrymple_polygons, composition = dalrymple_composition, rate = 15,
                        reals = 5,
                        method.sample = "by_area", 
                        method.allocate = "weighted",
                        method.model = NULL,
                        strata = NULL,
                        outputdir = getwd(), stub = "C5", cpus = 4,
                        factors = NULL)


# Do a raster classification (7 classes)
load(file = "C:/rdev/dsmart/rPackage/dsmart/pkg/data/dalrymple_covariates.rda")
dalrymple_covariates[[1]]
m <- c(300, 325, 1,  325, 350, 2,  350, 375, 3, 375,400,4,400,425,5,425,450,6,450,475,7)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- as.factor(reclassify(dalrymple_covariates[[1]], rclmat))
dalrymple_covariates[[1]]<- rc

# Do a raster classification (3 classes)
load(file = "C:/rdev/dsmart/rPackage/dsmart/pkg/data/dalrymple_covariates.rda")
dalrymple_covariates[[1]]
m <- c(300, 350, 1,  350, 400, 2,  400, 500, 3)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- as.factor(reclassify(dalrymple_covariates[[1]], rclmat))
dalrymple_covariates[[1]]<- rc
names(dalrymple_covariates[[1]])<- "classes1"

dalrymple_covariates[[8]]
m <- c(0, 0.3, 1,  0.3, 0.6, 2,  0.6, 1, 3)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- as.factor(reclassify(dalrymple_covariates[[8]], rclmat))
dalrymple_covariates[[8]]<- rc
names(dalrymple_covariates[[8]])<- "classes2"





disag_out<-disaggregate(covariates = dalrymple_covariates, polygons = dalrymple_polygons, composition = dalrymple_composition, rate = 15,
                        reals = 3,
                        method.sample = "by_area", 
                        method.allocate = "weighted",
                        strata = NULL,
                        outputdir = getwd(), stub = "factor", cpus = 3,
                        factors = c("classes1", "classes2"))


