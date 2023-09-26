.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Writing of raster files may fail under GDAL >= 3 and PROJ >= 6! This will be fixed in a future release of the rdsmart package.")
}