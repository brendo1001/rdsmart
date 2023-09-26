library(rgdal); library(raster); library(sp); library(gtools); library(C50)


#Files location
setwd("/home/brendo/myWork/dsmart/data")
# Identify map unit polygon shapefile (attribute table structured the same way as for Python DSMART)
shapefile = "dlr_polys_alb_7km.shp"
# Map unit composition file
compositionFile = "formatted_attribute_table2.txt"


# Load data
# Load shapefile
shp = readOGR(dsn=paste0("./polygons/", shapefile), layer=substr(shapefile, 1, nchar(shapefile) - 4))
polygons<- shp


#########################################################################################################3
# Load covariates as raster stack
covariates=stack()
for(filename in dir(path=paste0(getwd(),"/covariates/"), pattern=".img$"))
{
  r = raster(paste0(getwd(),"/covariates/",filename))
  # Load raster to stack
  covariates=stack(covariates, raster(paste0(getwd(),"/covariates/",filename)))
}

setwd("/home/brendo/myWork/dsmart/data/covariates")
fileN<- list.files(getwd(),  pattern="img$", full.names=FALSE)
files<- list.files(getwd(), pattern='img$',full.names=T)

#Load rasters into memory
s12<- stack()
for (i in 1:length(files)){
d1<- readGDAL(files[i])
r1<- raster(d1)
names(r1)<- str_sub(fileN[i],1,-5)
s12<- stack(s12,r1)}

##################################################################################################################
# Load map unit composition text file
composition = read.table(paste0(getwd(), "/polygons/", compositionFile), sep=",", header=TRUE)
colnames(composition) = c("poly","mapunit","soil_class","proportion")


####### DSMART tree creation #######
dsmart(covariates = covariates, polygons = shp, composition = composition, n=15, reals = 5, cpus=4)

####### DSMART probability rasters #######
setwd("/home/brendo/myWork/dsmart/data/dsmartOuts/rasters")
list.files(getwd(),  pattern="tif$", full.names=FALSE)
files<- list.files(getwd(), pattern='tif$',full.names=T)

s8<- stack()
for (i in 1:length(files)){
  r1<- raster(files[i])
  s8<- stack(s8, r1)}
lookup <- read.table("classLookupTable.txt",sep=",", header=TRUE)



test1<- dsmartR(rLocs= s8, nprob = 3, sepP=TRUE, lookup = lookup,cpus=2, param=nrow(lookup), param2 = nlayers(s8))


####Prepare data for inclusion in package
#composition file
dsT_composition <- composition
save(dsT_composition, file= "/home/brendo/myWork/dsmart/rPackage/dsmart/pkg/data/dsT_composition.rda")

#Polygon File
dsT_polygons <- polygons
save(dsT_polygons, file= "/home/brendo/myWork/dsmart/rPackage/dsmart/pkg/data/dsT_polygons.rda")

#covariates File
dsT_covariates <- s12
save(dsT_covariates, file= "/home/brendo/myWork/dsmart/rPackage/dsmart/pkg/data/dsT_covariates.rda")

#Lookup table
dsT_lookup <- lookup
save(dsT_lookup, file= "/home/brendo/myWork/dsmart/rPackage/dsmart/pkg/data/dsT_lookup.rda")


#Load dsmart rasters into memory
fileN<- list.files(getwd(),  pattern="img$", full.names=FALSE)
files<- list.files(getwd(), pattern='img$',full.names=T)

#Load rasters into memory
setwd("/home/brendo/myWork/dsmart/data/dsmartOuts/rasters")
fileN<-list.files(getwd(),  pattern="tif$", full.names=FALSE)
files<- list.files(getwd(), pattern='tif$',full.names=T)
files
s13<- stack()
for (i in 1:length(files)){
  d1<- readGDAL(files[i])
  r1<- raster(d1)
  names(r1)<- str_sub(fileN[i],1,-5)
  s13<- stack(s13,r1)}

#covariates File
dsmartOutMaps <- s13
save(dsmartOutMaps, file= "/home/brendo/myWork/dsmart/rPackage/dsmart/pkg/data/dsmartOutMaps.rda")


#examples
  library(dsmart)
  library(raster)
  
  #load the relevent data
  
  #Covariates
  data(dsT_covariates)
  
  #Polygons
  data(dsT_polygons)
  
  #Map unit compositions
  data(dsT_composition)
  
  #Run dsmart (with 15 samples per polygon, 5 C5 model realisations, using 1 compute node)
  testRun<- dsmart(covariates = dsT_covariates, polygons = dsT_polygons, composition = dsT_composition, obsdat=obsLocs,  n=15, reals = 5, cpus=1)
  testRun1<- dsmart(covariates = dsT_covariates, polygons = dsT_polygons, composition = dsT_composition, obsdat=NULL,  n=15, reals = 5, cpus=1)




###Making a small dataset of observed data to include into the DSMART Algorithm (added 290716)
training <- sample(nrow(locs), 0.7 * nrow(locs))
obsLocs<- locs[training,]
#covariates File
dsmartOutMaps <- s13
save(obsLocs, file= "C:/rdev/dsmart/rPackage/dsmart/pkg/data/obsLocs.rda")
