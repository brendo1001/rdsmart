#dsmart

# Randomly samples spatial polygons and builds a See5 classification tree
# with the sampling points.
#
# Args:
#   covariates: A RasterStack of covariate layers.
#   polygons: A SpatialPolygonsDataFrame containing the polygons to be
#      randomly sampled.
#   composition: A data frame containing the map unit composition for each
#      map unit polygon. First column is the polygon number (corresponding
#      to the first field of the SpatialPolygonsDataFrame attribute table);
#      second column is the map unit code; third column is the soil class
#      code; fourth column is the areal proportion of the map unit the soil
#      class is assumed to occupy.
#   n: The number of samples to draw from each polygon. With area-proportional
#      sampling n indicates the number of samples to draw per square kilometer
#      (for meter-based and latitude/longitude coordinate systems).
#   reals: Number of model resamplings to execute
#   cpus: Number of compute nodes to use
#   factors: Character vector with the names of the covariates that the
#      model should treat as factors.
#   sampling: Name of the sampling scheme. Per-polygon sampling ("PP")
#      draws a set number of samples for each polygon. Area-proportional
#      sampling ("AP") draws samples in proportion to the areas of the
#      polygons.
#   minrate: Minimum number of samples to draw from a polygon with area-
#      proportional sampling.

# Returns:
#   A number of items saved to file:
#          1. R date object that holds to model parameters from each fitted C5 model
#          2. Text files of the model summary output from each fitted C5 model
#          3. Rasters of each C5 model realisation. 
#          4. class and unique number lookup table
#   


dsmart<-function(covariates = NULL, polygons = NULL, composition = NULL, n=NULL, reals = NULL, cpus=1,factors = NULL,sampling = "PP",minrate = 0){
  beginCluster(cpus)
  # Generate lookup table
  lookup = as.data.frame(sort(unique(composition$soil_class)))
  lookup$code = seq(from=1, to=nrow(lookup), by=1)
  colnames(lookup) = c("name", "code")
  
  #create output repositories
  model_lists<- vector("list", reals) #empty list 
  dir.create("dsmartOuts/",showWarnings = F)
  dir.create("dsmartOuts/rasters",showWarnings = F)
  dir.create("dsmartOuts/models",showWarnings = F)
  dir.create("dsmartOuts/summaries",showWarnings = F)
  strg<- paste(getwd(),"/dsmartOuts/rasters/",sep="")
  strm<- paste(getwd(),"/dsmartOuts/models/",sep="")
  strs<- paste(getwd(),"/dsmartOuts/summaries/",sep="")
  write.table(lookup, paste(strg,"classLookupTable.txt",sep=""),sep=",", col.names=T,row.names=F) 
  
  # For area-proportional sampling, calculate polygon areas in kilometers squared
  if(sampling == "AP"){
  areas <- raster::area(polygons)/1e6
  sample.rate <- n
  number <- 0
  }
  
  for (j in 1:reals){
    # Empty data frame to store samples
    coordF<- matrix(NA, nrow=1000, ncol=3)
    coordF<- data.frame(coordF)
    names(coordF)<- c("x", "y", "class")
    cf<- 1
    for(poly.id in polygons@data[,1])
    {
      #print(poly.id)
      
      # For area-proportional sampling, calculate number of samples to draw
      if(sampling == "AP"){
      number <- number + 1
      n <- ceiling(areas[number] / sample_rate))
      n <- max(minrate, n)
      }
      
      # Subset a single polygon
      poly = subset(polygons, polygons@data[,1]==poly.id)
      coordF[cf:(cf+(n-1)),1:2] = as.data.frame(spsample(poly, n , type="random", iter=10))
      
      # Allocate soil classes from within map unit
      poly.comp=subset(composition, composition$poly==poly.id)
      # Draw from Dirichlet distribution
      s=rdirichlet(1, poly.comp$proportion)
      
      # Weighted-random sample
      coordF$class[cf:(cf+(n-1))]=sample(poly.comp$soil_class, size=n, replace=TRUE, prob=s[1,])
      cf<- cf+n}
      #spatial object
      locs<- as.data.frame(coordF[complete.cases(coordF),])
      coordinates(locs)<- ~ x + y  
      # Extract covariate values for the sampling locations
      values=extract(covariates,locs)
    
      #sample frame
      samples = cbind(as.data.frame(values),as.data.frame(locs)[,3])  
      names(samples)[ncol(samples)]<- "soil_class"
      samples$soil_class<- as.factor(samples$soil_class)
      
      #Convert designated covariates to factors
      if(is.character(factors)){
      fcols <- c(1:ncol(samples))[colnames(samples) %in% factors]
      frasters <- c(1:length(names(covariates)))[names(covariates) %in% factors]
      for(i in 1:length(fcols)){
        samples[,fcols[i]]<-factor(samples[,fcols[i]],levels = as.character(levels(var.stack[[frasters[i]]])[[1]]$Value))
      }}
      
      #Fit model####
      res = C5.0(samples[,-ncol(samples)], y=samples$soil_class)
      model_lists[[j]]<- res
      #Capture output
      out<-capture.output(summary(res))
      f2<- paste(strm,paste("C5_model_", j,".txt",sep=""),sep="" )
      cat(out,file=f2,sep="\n",append=TRUE)
    
    nme<- paste(paste(paste(strg,"map",sep=""),"_",j,sep=""), ".tif", sep="")
    r1 <- clusterR(covariates, predict, args=list(res),filename=nme,format="GTiff",overwrite=T, datatype="INT2S")}
  
  #Save models to file
  save(model_lists, file = paste(paste(getwd(),"/dsmartOuts/",sep=""),"dsmartModels.RData", sep="") )
  endCluster()
  message(paste(paste("DSMART outputs can be located at:",getwd(), sep=" "), "/dsmartOuts/",sep="") )}

#END