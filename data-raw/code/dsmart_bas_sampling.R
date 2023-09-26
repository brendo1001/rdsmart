library(dplyr)
library(data.table)
library(doParallel)
library(magrittr)
library(microbenchmark)
library(raster)
library(sf)
library(sp)

# Load shapefile
lowlands_smap <- 
  st_read("D:/LCR/projects/682208-0136 Hill Country DSM/data/waikato/smap/lowlands_smap.shp") %>% 
  mutate(polygon = seq(1, nrow(.))) %>% 
  dplyr::select(polygon, everything())


realisations <- 50
rate <- 15




cl <- makeCluster(6)
registerDoParallel(cl)

x <- foreach(i = seq(1, realisations),
             .packages = c('data.table', 'dplyr', 'sf', 'SDraw')) %dopar%
{
  # List to hold intermediate output
  output <- list()
  
  # Loop through polygons
  for(j in unique(lowlands_smap$OBJECTID)) {
    s <-
      lowlands_smap %>% 
      dplyr::filter(OBJECTID == j) %>% 
      as("Spatial") %>% 
      SDraw::bas.polygon(15) %>% 
      st_as_sf()
    
    output <- append(output, list(s))
  }
  
  # Merge polygon samples into one sf data frame
  output <- data.table::rbindlist(output) %>%
    st_sf() %>% 
    mutate(realisation = i)
  
  return(output)
}

stopCluster(cl)

x <- 
  rbindlist(x) %>% 
  #set_colnames(c("polygon", "geometry")) %>% 
  st_sf()

# Load covariates
covariates <- 
  list.files("D:/LCR/projects/682208-0136 Hill Country DSM/data/waikato/covariates/", pattern = ".tif$", full.names = TRUE) %>% 
  stack()

# Get covariates and remove rows where there is NA in the covariates
x_covariates <- 
  raster::extract(covariates, as(x, "Spatial")) %>% 
  cbind(x, .) %>% 
  na.omit()

# Sample 
x_samp <- 
  x_covariates %>%
  group_by(polygon) %>% 
  sample_n(size = (rate * realisations), replace = TRUE) %>% 
  mutate(i = rep(seq(1, realisations), rate)) %>% 
  arrange(polygon, i) %>% 
  dplyr::select(polygon, i, everything())



bas_x <- bas.polygon(as(lowlands_smap, "Spatial"), 15)
nrow(bas_x)

lowlands_smap %>% 
  filter(OBJECTID == 11950) %>% 
  as("Spatial") %>% 
  bas.polygon(15) %>% 
  st_as_sf()





################################################################################

microbenchmark({
# List to hold intermediate output
output <- list()

# Loop through polygons
for(j in unique(lowlands_smap$OBJECTID)) {
  s <-
    lowlands_smap %>% 
    dplyr::filter(OBJECTID == j) %>% 
    as("Spatial") %>% 
    SDraw::bas.polygon(15) %>% 
    st_as_sf()
  
  output <- append(output, list(s))
}

# Merge polygon samples into one sf data frame
output <- 
  data.table::rbindlist(output) %>%
  st_sf()

}, times = 1)

# 32 seconds
# 34

output <- rbind(output, output2)

microbenchmark({
# Get covariates and remove rows where there is NA in the covariates
output_covariates <- 
  raster::extract(covariates, as(output, "Spatial")) %>% 
  cbind(dplyr::select(output, polygon), .)
}, times = 1)
# 89 seconds
# 90 seconds

# Remove samples that have NA in the covariates
output_covariates <-
  output_covariates %>% 
  st_set_geometry(NULL) %>% 
  complete.cases() %>% 
  filter(output_covariates, .)

# Remove map units that do not have the right number of samples
n <- output_covariates %>% 
  group_by(polygon) %>% 
  summarise(n = n()) %>% 
  filter(n == 15) %>% 
  st_set_geometry(NULL) %>% 
  left_join(output_covariates, by = "polygon")






################################################################################



# Load covariates
covariates <- 
  list.files("D:/LCR/projects/682208-0136 Hill Country DSM/data/waikato/covariates/", pattern = ".tif$", full.names = TRUE) %>% 
  stack()


cl <- makeCluster(5)
registerDoParallel(cl)

microbenchmark({
x <- foreach(i = seq(1, 100),
             .packages = c('raster')) %dopar%
             {
               return(covariates[[1]][5, i])
             }
}, times = 10)

stopCluster(cl)
