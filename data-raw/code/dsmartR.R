#dsmartR

###Based on the C5 tree modelling the function determines:
#1. The pixel based probability to each class
#2. The n most probable classes at each pixel
##All outputs are save to file in raster format
# Class probability maps
# n most probable maps
#R pbjects are saved:
#1. The probability rasters (rasterStack) is requested
#2. The n most probable classes (raster Stack)
#FUnction Requires:
# rLocs: The locations of the rasters where the dsmart algorithm deposited its outputs
# nprob: the n most probable class maps to be produced
# sepP: logical of whether class probability maps should be produced
# lookup: lookup table produced from dsmart that numerically links soil class codes to a number

dsmartR<- function(rLocs= NULL, nprob = NULL, sepP=FALSE, lookup= NULL){
  
  setwd(rLocs)
  dir.create("counts/",showWarnings = F)
  dir.create("probabilities/",showWarnings = F)
  dir.create("nProbable/",showWarnings = F)
  strc<- paste(getwd(),"/counts/",sep="")
  strp<- paste(getwd(),"/probabilities/",sep="")
  strn<- paste(getwd(),"/nProbable/",sep="")
  lookup = read.table("classLookupTable.txt", sep=",", header=TRUE)
  
  
  #List the rasters
  #list.files(strg,  pattern="tif$", full.names=FALSE)
  files<- list.files(rLocs, pattern='tif$',full.names=T)
  
  #Make a stack of the rasters
  s1<- stack()
  for (kk in 1: length(files)){
    r1<- raster(files[kk])
    s1<- stack(s1,r1)}
  
  
  #counts
  nme1<- paste(strc,"countOuts.tif" ,sep="") 
  counts = calc(s1, function(s1) tabulate(s1, nbins=nrow(lookup)), filename=nme1,format="GTiff",progress="text",overwrite=T)
  
  
  #probabilities
  nme2<- paste(strp,"countOutsPropbs.tif" ,sep="") 
  probs= calc(counts, function(counts) (counts/nlayers(s1)),filename=nme2,format="GTiff",progress="text",overwrite=T )
  
  if (sepP==TRUE) {s3<- stack()
                   for (np in 1:nlayers(probs)){
                     nme5<- paste(paste(strp, as.character(lookup[np,1]), sep=""), "_probs.tif", sep="")
                     names(probs[[np]])<- as.character(lookup[np,1])
                     s3<- stack(s3,probs[[np]])
                     writeRaster(probs[[np]],filename=nme5,format="GTiff",progress="text",overwrite=T)}}
  
  #Most probable
  nme3<- paste(strn,"nProbable.tif",sep="")
  ordered.indices = calc(counts, fun=function(counts) order(counts, decreasing=TRUE, na.last=TRUE),filename=nme3,format="GTiff",progress="text",overwrite=T)
  s4<- stack()
  for (zz in 1:nprob){
    nme4<- paste(strn,paste(zz,"_probable.tif",sep=""),sep="")
    s4<- stack(s4,ordered.indices[[zz]])
    writeRaster(ordered.indices[[zz]],filename=nme4,format="GTiff",progress="text",overwrite=T)}
  
  if (sepP==TRUE){retval<-list(s3,s4)} else{ retval<-list(s4)}
  return(retval)}

#END