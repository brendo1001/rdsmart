root<-"/home/brendo1001/mywork"

datdir<-paste0(root,"/dsmart/rPackage/dsmart/pkg/data/")

setwd(datdir)

file_names1<-as.list(dir(pattern="rda"))

lapply(file_names1,load,.GlobalEnv)

sourcedir<-paste0(root,"/dsmart/rPackage/dsmart/pkg/R/")

setwd(sourcedir)

file_names2<-as.list(dir())

lapply(file_names2,source)

rm(file_names1,file_names2)

setwd(paste0(root,"/muddles/"))

library(C50)
library(ranger)
library(foreach)
library(caret)
library(magrittr)
library(raster)

# 1: Test the bugfix relating to caret::train predictions with raw class predictions (rate set 
# to 1 to ensure that levels are dropped from soil_class)

set.seed(1)

test1<-dsmart(  covariates = dalrymple_covariates
              , polygons = dalrymple_polygons
              , composition = dalrymple_composition
              , rate = 1
              , reals = 3
              , method.model = "ranger"
              , args.model = list(  num.trees = 100
                                  , trControl = trainControl(method = "oob")
                                  )
              )

# Seems to work :)

# 2: Test probabilistic predictions with C5.0.

set.seed(1)

test2<-dsmart(  covariates = dalrymple_covariates
              , polygons = dalrymple_polygons
              , composition = dalrymple_composition
              , rate = 1
              , reals = 3
              , type = "prob"
              )

br<-brick(paste0(root,"/muddles/output/realisations/realisation_1.tif"))
br
plot(br[[1]])
plot(br[[2]])

# Works as well :)

# 3: Test probabilistic predictions with caret.

set.seed(1)

test3<-dsmart(  covariates = dalrymple_covariates
                , polygons = dalrymple_polygons
                , composition = dalrymple_composition
                , rate = 15 # Needs to be fairly high, otherwise ranger fails because of rare classes.
                , reals = 3
                , method.model = "ranger"
                , args.model = list(  num.trees = 100
                                    , trControl = trainControl(  method = "LGOCV"
                                                                 , number = 10
                                                                 , p = 0.9
                                                                 , classProbs = TRUE
                                                                 , savePredictions = "final"
                                                               ))
                , type = "prob"
                )

br<-brick(paste0(root,"/muddles/output/realisations/realisation_1.tif"))
br
plot(br[[1]])
plot(br[[2]])

# Booyakasha.

# 4: Test Test probabilistic predictions with caret; only one realisation.

set.seed(1)

test4<-dsmart(  covariates = dalrymple_covariates
                , polygons = dalrymple_polygons
                , composition = dalrymple_composition
                , rate = 15
                , reals = 1
                , method.model = "ranger"
                , args.model = list(  num.trees = 100
                                      , trControl = trainControl(  method = "LGOCV"
                                                                   , number = 10
                                                                   , p = 0.9
                                                                   , classProbs = TRUE
                                                                   , savePredictions = "final"
                                      ))
                , type = "prob"
                )


# Works just fine :)
