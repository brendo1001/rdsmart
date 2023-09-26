# rdsmart 2.0.28

* Created `NEWS.md`.
* Updated doi links in documentation for `disaggregate()`, `dsmart()` and `summarise()` to https.
* Implemented `shannon_entropy()` to compute Shannon's entropy (rather than compute it inside `summarise()`) and fixed a bug in its previous implementation.
* Improved documentation for `confusion_index()`.
* The function `summarise()` now uses `confusion_index()` and `shannon_entropy()` to compute the confusion index and Shannon's entropy rather than computing them itself.
* Implemented `.onAttach()` to display a message when the package is loaded to explain that writing rasters may fail under GDAL >= 3 and PROJ >= 6.