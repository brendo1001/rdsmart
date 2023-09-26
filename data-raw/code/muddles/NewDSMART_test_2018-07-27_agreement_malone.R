## test out puts of dsmart

library(raster);library(sp);library(rgdal)


#load most probable rasters

r1<- raster("/home/brendo1001/mywork/muddles/test1/mostprobable/mostprob_01_class.tif")
r2<- raster("/home/brendo1001/mywork/muddles/test2/mostprobable/mostprob_01_class.tif")
r3<- raster("/home/brendo1001/mywork/muddles/test3/mostprobable/mostprob_01_class.tif")
r4<- raster("/home/brendo1001/mywork/muddles/test4/mostprobable/mostprob_01_class.tif")

s1<- stack(r1,r2,r3,r4)
s1

# convert to data.frame
# Covert rasters to table
tempD <- data.frame(cellNos = seq(1:ncell(s1)))
vals <- as.data.frame(getValues(s1))
tempD<- cbind(tempD, vals)
tempD <- tempD[complete.cases(tempD), ]
cellNos <- c(tempD$cellNos)
gXY <- data.frame(xyFromCell(s1, cellNos, spatial = FALSE))
tempD<- cbind(gXY, tempD)
str(tempD)



# check agreement
# 100% test1,2,3,4
xx<- c()
for(i in 1:nrow(tempD)){
  xx[i]<- sum(tempD[i,4]== tempD[i,4:7])/4}
sum(xx==1)/nrow(tempD)

# 100% test 1 and 3
xx<- c()
for(i in 1:nrow(tempD)){
  xx[i]<- sum(tempD[i,4]== tempD[i,6])}
sum(xx==1)/nrow(tempD)







