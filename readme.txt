# This is a place to put working notes for r package development of dsmart

##12.5.15
dsmart has complete parallelisation and can be controlled by user
dsmartR was originally in sequential but now can be run completely in parallel.
NEED TO DO:
1. Write up the manuals for the companion functions
2. Write description and namespace files

##14.05.15
Writing up the manuals presently for the functions and data.
Did first build of package. Worked ok on unix. Didnt work for windows.
Built package on windows computer ok
Some manual typos attended to.

##17.05.15
Package seems to be building ok.
Need to:
	1. Calculate probabilties of n-probables
	2. Calculate confusion index.


##18.05.15
Add the functionality of the previous days tasks (both 1 and 2)
Add progress bars

##27.05.15
Tried function on data set from NZ. Seems a bit slower than python version
Improvements to make:
	1. paralellise the extraction process on dsmart
	2. make more general the composition table input.
