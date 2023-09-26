#' Disaggregating and harmonising soil map units through resampled
#' classification trees
#'
#' This function, together with companion function \code{summarise} implements the
#' DSMART (Disaggregating and harmonising soil map units through resampled
#' classification trees) algorithm as described in Odgers et al. (2014). This is
#' the workhorse function that involves multiple resampling, C5 decision tree
#' model fitting, and subsequent mapping in order to realize potential candidate
#' soil classes within aggregated soil mapping units. There is also added
#' facility to incorporate point observed data into the algorithm too.
#'
#' There is no prescription for the layers that compose \code{covariates} other
#' than that they should represent the \emph{scorpan} factors (McBratney
#' \emph{et al.}, 2003) as effectively as possible. The order of the layers in
#' the RasterStack is not important. See \code{data(dsT_covariates)} for an
#' example.
#'
#' DSMART assumes, but does not currently check, that \code{covariates},
#' \code{polygons} and any \code{observations} are projected to the same
#' coordinate reference system.
#'
#' DSMART assumes that the soil classes in \code{composition} and those in
#' \code{observations} belong to the same taxonomic level of the same soil
#' classification system. For example if the soil classes in \code{composition}
#' are all soil series, those in \code{observations} should also be series
#' rather than, for example, orders or great groups.
#'
#' \code{dsmart} produces a number of outputs which are saved to subdirectories
#' in the current working directory. The base folder for the outputs is
#' \code{output}. Inside this folder, rasters of the realisations are saved into
#' the \code{realisations} subfolder and the classification models are saved
#' into the \code{models} subfolder.
#'
#' @param covariates A \code{RasterStack} of \emph{scorpan} environmental
#'   covariates to calibrate the \code{C50} classification trees against. See
#'   \emph{Details} for more information.
#' @param polygons A \code{SpatialPolygonsDataFrame} containing the soil map
#'   unit polygons that will be disaggregated. The first field of the data frame
#'   must be an integer that identifies each polygon.
#' @param composition A \code{data.frame} that contains information on the
#'   soil-class composition of each polygon in \code{polygons}. Each row
#'   contains information about one soil class component of one polygon, which
#'   belongs to one soil map unit. First field contains the integer that
#'   identifies the polygon. Second field contains a code that identifies the
#'   soil map unit that the polygon belongs to.
#'
#'   If \code{strata = NULL} (the default), third column contains a code that
#'   identifies the soil class and fourth column contains a number in the range
#'   \code{(0, 100)} that identifies the proportion of the \strong{map unit}
#'   that the soil class corresponds to. See the example data
#'   \code{data(dalrymple_composition)}.
#'
#'   If \code{strata} is a \code{RasterLayer}, third column contains an integer
#'   that identifies the stratum in \code{strata}, fourth column contains a code
#'   that identifies the soil class and fifth column contains a number in the
#'   range \code{(0, 100)} that identifies the proportion of the
#'   \strong{stratum} that the soil class corresponds to.
#' @param rate An integer that identifies the number of virtual samples to draw
#'   from each polygon in each realisation. If \code{method.sample =
#'   "by_polygon"}, the number of samples to draw from each polygon in
#'   \code{polygons}. If \code{method.sample = "by_area"}, the sampling density
#'   in number of samples per square kilometre.
#' @param reals An integer that identifies the number of realisations of the
#'   soil class distribution that DSMART should compute.
#' @param observations \emph{Optional} A \code{data.frame} that contains actual
#'   observations of the soil class at locations across the soil map area. These
#'   data augment the virtual samples and are used in each realisation. Each row
#'   contains information about one soil class observation. First and second
#'   fields contain the \emph{x-} and \emph{y-}components of the observation's
#'   spatial location. Third field is the soil class code. See \emph{Details}.
#' @param method.sample Identifies the sampling method. Valid values are
#'   \code{"by_polygon"} (the default), in which case the same number of samples
#'   are taken from each polygon; or \code{"by_area"}, in which case the number
#'   of samples per polygon depends on the area of the polygon.
#' @param method.allocate Method of allocation of virtual samples to soil
#'   classes. Valid values are \code{"weighted"}, for weighted-random allocation
#'   to a soil class from within the virtual sample's map unit;
#'   \code{"random_mapunit"}, for completely random allocation to a soil class
#'   from within the virtual sample's map unit; and \code{"random_all"}, for
#'   completely random allocation to a soil class from within the entire map
#'   area.
#' @param method.model Method to be used for the classification model. If no
#'   value is passed, a C5.0 decision tree is built. Otherwise, the value must
#'   match a valid 'method' argument in the caret::train function.
#' @param method.args A list of arguments to be passed to the caret::train
#'   function. The list can include a trainControl object, arguments to be
#'   passed directly to the train function and arguments to be passed to the
#'   predictive model.
#' @param strata \emph{optional} An integer-valued \code{RasterLayer} that will
#'   be used to stratify the allocation of virtual samples to soil classes.
#'   Integer values could represent classes of slope position (e.g. crest,
#'   backslope, footslope, etc.) or land use (e.g. cropland, native vegetation,
#'   etc.) or some other variable deemed to be an important discriminator of the
#'   occurrence of soil classes within a map unit.
#' @param outputdir A character string that identifies the location of the main
#'   output directory. The folder \code{output} and its subfolders will be
#'   placed here. Default is the current working directory, \code{getwd()}.
#' @param stub \emph{optional} A character string that identifies a short name
#'   that will be prepended to all output.
#' @param cpus An integer that identifies the number of CPU processors to use
#'   for parallel processing.
#' @param factors A character vector with the names of the covariates that
#'   should be treated as factors.
#' @param prob A character vector to specify the type of the predictions. By
#'   default, raw class predictions are used. If set to "prob", a rasterbrick
#'   with class probabilities will be produced for each realisation.
#'
#' @return A list that contains metadata about the current run of
#'   \code{disaggregate}.
#'
#' @examples
#' # Load datasets
#' data(dalrymple_composition)
#' data(dalrymple_covariates)
#' data(dalrymple_observations)
#' data(dalrymple_polygons)
#'
#' # Run disaggregate without adding observations
#' disaggregate(dalrymple_covariates, dalrymple_polygons, dalrymple_composition,
#'  rate = 15, reals = 10, cpus = 6)
#'
#' # Run disaggregate with extra observations
#' disaggregate(dalrymple_covariates, dalrymple_polygons, dalrymple_composition,
#'  observations = dalrymple_observations, rate = 15, reals = 10)
#'
#' @references McBratney, A.B., Mendonca Santos, M. de L., Minasny, B., 2003. On
#'   digital soil mapping. Geoderma 117, 3--52. doi:
#'   \href{https://doi.org/10.1016/S0016-7061(03)00223-4}{10.1016/S0016-7061(03)00223-4}
#'   
#'   Odgers, N.P., McBratney, A.B., Minasny, B., Sun, W., Clifford, D., 2014.
#'   DSMART: An algorithm to spatially disaggregate soil map units, \emph{in:}
#'   Arrouays, D., McKenzie, N.J., Hempel, J.W., Richer de Forges, A.,
#'   McBratney, A.B. (Eds.), GlobalSoilMap: Basis of the Global Spatial Soil
#'   Information System. Taylor & Francis, London, pp. 261--266.
#'
#'   Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., 2014.
#'   Disaggregating and harmonising soil map units through resampled
#'   classification trees. Geoderma 214, 91--100. doi:
#'   \href{https://doi.org/10.1016/j.geoderma.2013.09.024}{10.1016/j.geoderma.2013.09.024}
#'
#' @export

disaggregate <- function(covariates, polygons, composition, rate = 15,
                         reals = 100, observations = NULL,
                         method.sample = "by_polygon", 
                         method.allocate = "weighted",
                         method.model = NULL, args.model = NULL,
                         strata = NULL,
                         outputdir = getwd(), stub = NULL, cpus = 1,
                         factors = NULL, type = "raw", predict = TRUE)
{
  # Create list to store output
  output <- base::list()
  
  # Save start time
  output$timing <- base::list(start = base::date())
  
  # Check arguments before proceeding
  messages <- c("Attention is required with the following arguments:\n")
  if(!(class(covariates) == "RasterStack"))
  {
    messages <- append(messages, "'covariates': Not a valid RasterStack.\n")
  }
  if(!(class(polygons) == "SpatialPolygonsDataFrame"))
  {
    messages <- append(messages, 
                       "'polygons': Not a valid SpatialPolygonsDataFrame.\n")
  }
  if(!(class(composition) == "data.frame"))
  {
    messages <- append(messages, "'composition': Not a valid data.frame.\n")
  }
  if(rate <= 0)
  {
    messages <- append(messages, "'n': Value must be greater than 0.\n")
  }
  if(reals <= 0)
  {
    messages <- append(messages, "'reals': Value be greater than 0.\n")
  }
  if(cpus <= 0)
  {
    messages <- append(messages, "cpus must be greater than 0.\n")
  }
  if(!(is.null(observations)))
  {
    if(!(class(observations) == "data.frame"))
    {
      messages <- append(messages, "'observations': Not a valid data.frame.\n")
    }
  }
  if(!(file.exists(outputdir)))
  {
    messages <- append(messages, "'outputdir': Output directory does not exist.\n")
  }
  if(!(is.null(strata)))
  {
    if(!(class(strata) == "RasterLayer"))
    {
      messages <- append(messages, "'strata': Not a valid RasterLayer.\n")
    }
  }
  if(!(is.null(method.model)))
  {
    if(!(is.character(method.model) & length(method.model) == 1))
    {
      messages <- append(messages, "'method.model' must be NULL or a single character value.\n")
    }
  }
  if(length(messages) > 1)
  {
    stop(messages)
  }
  
  # If method.model is a character value, load the caret package
  if(is.character(method.model) & length(method.model) == 1){
  require(caret)
  }
  
  # Set stub to "" if NULL
  if(is.null(stub))
  {
    stub <- ""
    
  } else if (stub == "") {
    
    stub <- ""
    
  } else if(!(substr(stub, nchar(stub), nchar(stub)) == "_")) {
    
    stub <- paste0(stub, "_")
  }
  
  # Save function call
  output$call <- base::match.call()
  
  # Save parameters
  output$parameters <- base::list(rate = rate, reals = reals,
                                  method.sample = method.sample,
                                  method.allocate = method.allocate,
                                  method.model = method.model, 
                                  args.model = args.model, stub = stub,
                                  cpus = cpus, factors = factors, type = type)
  
  # Create subdirectories to store results in
  outputdir <- file.path(outputdir)
  dir.create(file.path(outputdir, "output"), showWarnings = FALSE)
  dir.create(file.path(outputdir, "output", "realisations"), showWarnings = FALSE)
  dir.create(file.path(outputdir, "output", "models"), showWarnings = FALSE)
  
  # Save output locations
  output$locations <- list(root = file.path(outputdir, "output"),
                           realisations = file.path(outputdir, "output", "realisations"),
                           models = file.path(outputdir, "output", "models"))
  
  # Rename composition column names
  if(!(is.null(strata))) {
    names(composition) <- c("poly", "mapunit", "stratum", "soil_class", "proportion")
  } else {
    names(composition) <- c("poly", "mapunit", "soil_class", "proportion")
  }
  
  # Since we only really need the first column of the polygons SpatialPolygonsDataFrame,
  # drop all its other attributes and rename the first column to "poly"
  polygons <-
    polygons %>%                                        # With polygons:
    sf::st_as_sf() %>%                                  # Cast to an sf data frame
    dplyr::select(1, geometry) %>%                      # Select just the id column and the geometry
    magrittr::set_colnames(c("poly", "geometry")) %>%   # Rename columns
    as("Spatial")                                       # Cast back to a SpatialPolygonsDataFrame

  # Make sure that there are no missing values in map unit composition
  polys_to_remove <-
    composition %>%                 # With the map unit composition:
    stats::complete.cases() %>%     # Identify complete cases
    magrittr::equals(FALSE) %>%     # Invert the complete cases identification
    base::which() %>%               # Find the rows of composition that are incomplete
    composition$poly[.] %>%         # Find the ids of the polygons that contain incomplete data
    base::unique()                  # Get the unique polygon ids

  # Remove the polygons that contain missing information from both the polygons
  # SpatialPolygonsDataFrame and the composition data frame
  if(length(polys_to_remove) > 0) {
    
    polygons <-
      polygons %>%                                      # With polygons:
      sf::st_as_sf() %>%                                # Cast to an sf data frame
      dplyr::filter(!(poly %in% polys_to_remove)) %>%   # Filter out the polygons whose ids are in polys_to_remove
      sf::st_sf() %>%                                   # Shouldn't have to do this step
      as("Spatial")                                     # Cast back to a SpatialPolygonsDataFrame

    composition <-
      composition %>%                                   # With composition:
      dplyr::filter(!(poly %in% polys_to_remove))       # Filter out the polygons whose ids are in polys_to_remove

    # Let the user know about what we've done
    warning(base::paste0("The following polygons were removed from further analysis because they have incomplete or undefined map unit compositions: ",
                 base::paste(base::as.character(polys_to_remove), collapse = ", ")))
  }
    
  
  # Write covariate names to file
  write.table(names(covariates), file.path(outputdir, "output",
                                           paste0(stub, "covariate_names.txt")),
              quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
  
  # Get samples for all realisations
  message(Sys.time(), " Generating samples for ", reals, " realisations")
  samples <- data.frame()
  
  if(is.null(strata))
  {
    # Get samples without allocation stratification
    samples <- .getVirtualSamples(covariates, polygons, composition,
                                  n.realisations = reals, rate = rate,
                                  method.sample = method.sample,
                                  method.allocate = method.allocate,
                                  cpus = cpus)
  } else {
    samples <- .getStratifiedVirtualSamples(covariates, polygons, composition,
                                            strata, n.realisations = reals,
                                            rate = rate, 
                                            method.sample = method.sample, 
                                            method.allocate = method.allocate)
  }
  
  # Make sure that there are no missing values in samples
  samples <-
    samples %>%                       # With the samples data frame:
    complete.cases() %>%              # Determine which records are complete
    dplyr::filter(samples, .)         # Filter for only the complete records
  
  # If there are observations, get their covariates
  if(!(is.null(observations)))
  { 
    names(observations) <- c("x", "y", "class")
    observations <- .observations(observations, covariates)
    write.table(observations, file.path(outputdir, "output",
                                        paste0(stub, "observations_with_covariates.txt")),
                sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  # Write samples to text file
  write.table(samples, file.path(outputdir, "output",
                                 paste0(stub, "virtual_samples.txt")),
              sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  # We submit the target classes to C5.0 as a factor. To do that, we need to 
  # make sure that the factor has the full set of levels, since due to the 
  # randomness of the allocation procedure we can't assume that they are all
  # represented in the samples drawn for a particular realisation. If the levels
  # are not specified correctly and they are not all represented in all sets of
  # samples, it is possible that the integer values that the class predictions
  # are coded to by raster::predict will not refer to the same soil class from
  # realisation to realisation.
  levs <- as.character(unique(composition$soil_class))
  
  if(!(is.null(observations)))
  {
    # Union map unit levels with observation levels
    levs <- base::union(levs, as.character(unique(observations$soil_class)))
  }
  
  # Sort levels
  levs <- sort(levs)
  
  # Generate lookup table
  lookup = as.data.frame(levs)
  lookup$code = seq(from=1, to=nrow(lookup), by=1)
  colnames(lookup) = c("name", "code")
  
  # Write lookup table to file
  write.table(lookup, file.path(outputdir, "output", paste0(stub, "lookup.txt")),
              sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
  output$locations$lookup <- file.path(outputdir, "output", paste0(stub, "lookup.txt"))
  
  # Compute tuning parameter for raster::clusterR
  tuning <- .blocks_per_node(raster::nrow(covariates),
                             raster::ncol(covariates),
                             cpus = cpus)
  
  # Process realisations
  for (j in 1:reals)
  {
    message(Sys.time(), " Realisation ", j)
    
    # Column number of samples' first covariate
    startcol <- numeric(0)
    if(is.null(strata))
    {
      startcol = 8
    } else {
      startcol = 9
    }
    
    # Extract sample covariates for the current realisation
    s <- samples[which(samples$realisation == j), startcol:ncol(samples)]
    soil_class <- as.character(samples$soil_class[which(samples$realisation == j)])
    
    # Concatenate observations if they exist
    if(!(is.null(observations)))
    {
      s <- rbind(s, observations[, 8:ncol(observations)])
      soil_class <- append(soil_class, as.character(observations$soil_class))
    }
    
    # Sort levels and convert soil_class back to factor
    soil_class <- factor(soil_class, levels = levs)

    # Convert designated covariates to factors.
    # if(is.character(factors)){
    #   fcols <- c(1:ncol(s))[colnames(s) %in% factors]
    #   for(i in 1:length(fcols)){
    #     s[,fcols[i]]<-as.factor(s[,fcols[i]])
    #   }}
    if(is.list(factors)) {
      for(f in names(factors)){
        s %<>% dplyr::mutate_at(f, factor, levels = factors[[f]])
      }
    }
    
    # Test if there are levels in soil_class with 0 cases
    # (this is used for a workaround, see below).
    zeroes <- sum(table(soil_class) == 0) > 0
    
    # If there are levels with 0 cases, make a reclassification table for the prediction raster.
    rclt<-cbind(c(1:sum(table(soil_class) != 0))
                , c(1:length(table(soil_class)))[table(soil_class) != 0]
    )

    # Fit model
    if(is.null(method.model)){
      model = C50::C5.0(s, y = soil_class)
    }else{
      # The train function does not accept factors in which there are levels which have 0 cases.
      # These levels must therefore be dropped.
      soil_class <- base::droplevels(soil_class)
      model = base::do.call(caret::train,
                            c(list(x = s,
                                   y = soil_class,
                                   method = method.model),
                              args.model))
    }
    
    # Save model to text file
    if(is.null(method.model)){
      out <- utils::capture.output(summary(model))
    }else{
      out <- utils::capture.output(model$finalModel)
    }
    cat(out, file = paste0(outputdir, "/output/models/", stub, "model_",
                           formatC(j, width = nchar(reals), format = "d",
                                   flag = "0"), ".txt"),
        sep = "\n", append = TRUE)
    
    # Save model to rdata file
    save(model, 
         file = file.path(outputdir, "output", "models", 
                          paste0(stub, "model_",
                                 formatC(j, width = nchar(reals), format = "d", flag = "0"),
                                 ".RData")))
    
    
    
    if(predict){
      if(type != "prob")
      {
        # If raw class predictions are specified (default), predict realisation and save it to 
        # raster.
        raster::beginCluster(cpus)
        r1 <- raster::clusterR(covariates, 
                               raster::predict,
                               args = list(model = model, factors = factors,
                                           inf.rm = TRUE),
                               m = tuning,
                               filename = file.path(outputdir, "output", "realisations",
                                                    paste0(stub, "realisation_", 
                                                           formatC(j, width = nchar(reals), format = "d", flag = "0"),
                                                           ".tif")),
                               format = "GTiff", overwrite = TRUE, datatype = "INT2S")
        
        # If levels were dropped from soil_class in order to use the train function,
        # the prediction raster must be reclassified in order to ensure that the integer
        # values represent the same soil types across the realizations.
        if(zeroes == TRUE & is.null(method.model) == FALSE)
        {
          r1 <- raster::clusterR(r1, reclassify, args = list(rcl = rclt),
                                 m = tuning,
                                 filename = file.path(outputdir, "output", "realisations",
                                                      paste0(stub, "realisation_", 
                                                             formatC(j, width = nchar(reals), format = "d", flag = "0"),
                                                             ".tif")),
                                 format = "GTiff", overwrite = TRUE, datatype = "INT2S")
        }
        raster::endCluster()
        
      } else {
        # If probabilistic predictions are specified, produce a raster brick with the
        # probabilities of each soil class.
        if(zeroes == FALSE | is.null(method.model)){
          # If no levels were dropped from soil_class in order to use the train function,
          # the rasterbrick can be used as it is.
          raster::beginCluster(cpus)
          r1 <- raster::clusterR(covariates, predict, 
                                 args = list(model, type = "prob",
                                             index = 1:nrow(lookup),
                                             na.rm = TRUE,
                                             factors = factors),
                                 m = tuning,
                                 filename = file.path(outputdir, "output", "realisations",
                                                      paste0(stub, "realisation_",
                                                             formatC(j, width = nchar(reals), format = "d", flag = "0"), ".tif")),
                                 format = "GTiff", overwrite = TRUE, datatype = "FLT4S")
          raster::endCluster()
          
        }else{
          # If levels were dropped from soil_class in order to use the train function, rasters
          # with probabilities for the missing levels must be inserted in the rasterbrick.
          # First, the rasterbrick is predicted as above.
          raster::beginCluster(cpus)
          tmp1 <- raster::clusterR(covariates,
                                   predict,
                                   factors = factors,
                                   args = list(model,
                                               type = "prob",
                                               index = 1:base::length(base::levels(soil_class)),
                                               na.rm = TRUE),
                                   m = tuning)
          # Then, a raster with 0's is calculated (the missing levels have 0 probability for the
          # realisation in question)
          tmp2 <- raster::clusterR(tmp1[[1]], calc, args = list(function(x) x*0), 
                                   m = tuning)
          raster::endCluster()
          # The rasters are then all arranged in a rasterstack, which is converted into the
          # final rasterbrick.
          rlist<-list()
          for(i in 1:nrow(lookup))
          {
            if(i %in% rclt[,2])
            {
              rlist[[i]]<-tmp1[[which(rclt[,2] == i)]]
            }else{
              rlist[[i]]<-tmp2
            }
          }
          rlist<-stack(rlist)
          r1<-raster::brick(rlist, filename = file.path(outputdir, "output", "realisations",
                                                        paste0(stub, "realisation_", formatC(j, width = nchar(reals), format = "d", flag = "0"),
                                                               ".tif")),
                            format = "GTiff", overwrite = TRUE, datatype = "FLT4S")
          rm(tmp1,tmp2)
        }
      }
    }
  }
  
  # Save finish time
  output$timing$finish <- base::date()
  
  # Return output
  return(output)
}

#END