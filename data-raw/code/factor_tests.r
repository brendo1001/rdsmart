library(C50)
library(raster)
library(rgdal)

##### Load point observations #####
obs <- readOGR(dsn = "D:/LCR/projects/337001-0005 SNutrient/fourpeaks/GIS/observations/site_obs_from_jim.shp",
               layer = "site_obs_from_jim")


##### Load covariates #####

covr <- stack()
for(filename in list.files(path = "D:/LCR/projects/337001-0005 SNutrient/fourpeaks/GIS/dsmart_covariates/",
                           pattern = ".tif$", full.names = TRUE))
{
  r <- raster(filename)
  
  # Load raster to stack
  covr <- stack(covr, r)
}

# Extract covariates for all points
obs.covr <- data.frame(obs[, c(1:6, 35)], extract(covr, obs, df = TRUE)) # 199 observations

# Fit tree to all observations
tree <- C5.0(obs.covr[, 12:ncol(obs.covr)], obs.covr$NZSC)
summary(tree)

# (a)   (b)   (c)   (d)   (e)   (f)   (g)    <-classified as
# ----  ----  ----  ----  ----  ----  ----
#         1                                  (a): class BFT
#        83           1     1     1     2    (b): class BOT
#         2     6           1                (c): class PIT
#         1           3                      (d): class RFT
#                          34           3    (e): class ROK
#               1                 2          (f): class ROT
#        28     1     1     2          25    (g): class RXT

tree$levels
#[1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"

levels(obs.covr$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"


# What if I remove all RXT from the calibration dataset?

which(obs.covr$NZSC == "RXT")
obs.covr.noRXT <- obs.covr[-which(obs.covr$NZSC == "RXT"),] # 142 observations

# Fit tree
tree.noRXT <- C5.0(obs.covr.noRXT[, 12:ncol(obs.covr.noRXT)], obs.covr.noRXT$NZSC)
summary(tree.noRXT)

# (a)   (b)   (c)   (d)   (e)   (f)   (g)    <-classified as
# ----  ----  ----  ----  ----  ----  ----
#          1                                  (a): class BFT
#         87                       1          (b): class BOT
#          2     7                            (c): class PIT
#                      4                      (d): class RFT
#          3                34                (e): class ROK
#                1                 2          (f): class ROT
#                                             (g): class RXT

tree.noRXT$levels
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"

levels(obs.covr.noRXT$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"

# CONCLUSION: Because the target classes given to C5.0 are a factor, RXT is
# included in the contingency table even though there are no observations of it,
# as long as it is still a factor level. To confirm, let's also remove RXT as a
# factor level.

levels(obs.covr.noRXT$NZSC) <- c("BFT", "BOT", "PIT", "RFT", "ROK", "ROT")
# That command doesn't work:
# Error in `levels<-.factor`(`*tmp*`, value = c("BFT", "BOT", "PIT", "RFT",  : 
#   number of levels differs

# To confirm:
levels(obs.covr.noRXT$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"

# If I give it the same number of levels but replace RXT with a duplicate RFT:
levels(obs.covr.noRXT$NZSC) <- c("BFT", "BOT", "PIT", "RFT", "ROK", "ROT", "RFT")
# This command works, and removes RXT from the levels
levels(obs.covr.noRXT$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT"

# Let's try building another tree now, and confirm that RXT is removed from the
# contingency table
tree.noRXT.2 <- C5.0(obs.covr.noRXT[, 12:ncol(obs.covr.noRXT)], obs.covr.noRXT$NZSC)
summary(tree.noRXT.2)

# (a)   (b)   (c)   (d)   (e)   (f)    <-classified as
# ----  ----  ----  ----  ----  ----
#          1                            (a): class BFT
#         87                       1    (b): class BOT
#          2     7                      (c): class PIT
#                      4                (d): class RFT
#          3                34          (e): class ROK
#                1                 2    (f): class ROT

tree.noRXT.2$levels
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT"

levels(obs.covr.noRXT$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT"










# Let's subset the observations into two groups: one containing all the RXT and 
# ROK, and the other containing all the other observations, and make sure that
# each subset's NZSC factor contains only the levels that are represented in it.

obs.covr.sub1 <- obs.covr[-which((obs.covr$NZSC == "RXT") | (obs.covr$NZSC == "ROK")), ] # 105 observations
obs.covr.sub2 <- obs.covr[which((obs.covr$NZSC == "RXT") | (obs.covr$NZSC == "ROK")), ] # 94 observations

levels(obs.covr.sub1$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"
# Still contains the full range of levels

levels(obs.covr.sub2$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"
# Still contains the full range of levels

# As an aside, let's make a tree based on obs.covr.sub2 and take a look at its
# contingency table, considering all but two of the classes are missing
tree.sub2.alllevels <- C5.0(obs.covr.sub2[, 12:ncol(obs.covr.sub2)], obs.covr.sub2$NZSC)
summary(tree.sub2.alllevels)

# (a)   (b)   (c)   (d)   (e)   (f)   (g)    <-classified as
# ----  ----  ----  ----  ----  ----  ----
#                                            (a): class BFT
#                                            (b): class BOT
#                                            (c): class PIT
#                                            (d): class RFT
#                          31           6    (e): class ROK
#                                            (f): class ROT
#                           2          55    (g): class RXT

# Just making sure:
tree.sub2.alllevels$levels
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"

levels(obs.covr.sub2$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"

# Also, just making sure:
levels(obs.covr.sub2$NZSC) <- c("ROK", "RXT")
# Error in `levels<-.factor`(`*tmp*`, value = c("ROK", "RXT")) : 
#   number of levels differs

# What happens to the factor levels if we convert obs.covr.sub2$NZSC to a 
# character and back to a factor?
obs.covr.sub2$NZSC <- as.character(obs.covr.sub2$NZSC)
class(obs.covr.sub2$NZSC)
# [1] "character"

levels(obs.covr.sub2$NZSC)
# NULL

obs.covr.sub2$NZSC <- as.factor(obs.covr.sub2$NZSC)
class(obs.covr.sub2$NZSC)
# [1] "factor"

# Aha!
levels(obs.covr.sub2$NZSC)
# [1] "ROK" "RXT"

# CONCLUSION: So, the way to rationalise a factor's levels is to turn it into
# a character vector and back to a factor again.

# What happens to the factor levels when we:
# 1. rbind two tables, each with rationalised factor levels?
# 2. rbind two tables, one with a full list of factors, but missing
#    representatives of some levels, and the other with the missing
#    representatives but a rationalised factor levels list?
# 3. rbind two tables, one with a full list of factors, but missing
#    representatives of some levels, and the other with the missing
#    representatives AND an extended factor levels list that contains levels
#    with no representatives in either table?
#
# Is there a difference between rbinding two tables with a factor NZSC column
# and a character NZSC column?

##### 1.

# First, rationalise the levels of obs.covr.sub1
obs.covr.sub1$NZSC <- as.character(obs.covr.sub1$NZSC)
obs.covr.sub1$NZSC <- as.factor(obs.covr.sub1$NZSC)

# To confirm:
levels(obs.covr.sub1$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROT"

levels(obs.covr.sub2$NZSC)
# [1] "ROK" "RXT"

# Good. Now rbind them.

obs.covr.sub12.1 <- rbind(obs.covr.sub1, obs.covr.sub2) # as expected, 199 observations
class(obs.covr.sub12.1$NZSC)
# [1] "factor"

levels(obs.covr.sub12.1$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROT" "ROK" "RXT"

# Good. The factor levels of the rbinded table are the total of the levels of
# the individual tables.



##### 2.

# First, reset obs.covr.sub1 so that it has the full list of levels in NZSC.
obs.covr.sub1 <- obs.covr[-which((obs.covr$NZSC == "RXT") | (obs.covr$NZSC == "ROK")), ] # 105 observations

levels(obs.covr.sub1$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"

levels(obs.covr.sub2$NZSC)
# [1] "ROK" "RXT"

# rbind the tables
obs.covr.sub12.2 <- rbind(obs.covr.sub1, obs.covr.sub2) # as expected, 199 observations
class(obs.covr.sub12.2$NZSC)
# [1] "factor"

# As we expect:
levels(obs.covr.sub12.2$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"



##### 3.

levels(obs.covr.sub2$NZSC)
# [1] "ROK" "RXT"

# Add more levels that are not represented in either table
levels(obs.covr.sub2$NZSC) <- c("ROK", "RXT", "XXA", "XXB", "XXC", "XXD")

# Good!
levels(obs.covr.sub2$NZSC)
# [1] "ROK" "RXT" "XXA" "XXB" "XXC" "XXD"

# Now let's rbind the tables.
obs.covr.sub12.3 <- rbind(obs.covr.sub1, obs.covr.sub2) # as expected, 199 observations

# That's interesting. So it actually merges the levels of both factors, instead
# of rationalising the levels to those that are actually represented in both
# original tables.
levels(obs.covr.sub12.3$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT" "XXA" "XXB" "XXC" "XXD"

# Let's make a tree with this table and examine the contingency table.
tree.sub12.3 <- C5.0(obs.covr.sub12.3[, 12:ncol(obs.covr.sub12.3)], obs.covr.sub12.3$NZSC)
summary(tree.sub12.3)

# Good!

# (a)   (b)   (c)   (d)   (e)   (f)   (g)   (h)   (i)   (j)   (k)    <-classified as
# ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
#         1                                                          (a): class BFT
#        83           1     1     1     2                            (b): class BOT
#         2     6           1                                        (c): class PIT
#         1           3                                              (d): class RFT
#                          34           3                            (e): class ROK
#               1                 2                                  (f): class ROT
#        28     1     1     2          25                            (g): class RXT
#                                                                    (h): class XXA
#                                                                    (i): class XXB
#                                                                    (j): class XXC
#                                                                    (k): class XXD








soilclass.1 <- obs.covr$NZSC

class(soilclass.1)
# [1] "factor"

levels(soilclass.1)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"

# Can we convert to a character direct from the data frame?
soilclass.2 <- as.character(obs.covr$NZSC)
class(soilclass.2)
#[1] "character"
# Yes we can.


# Extract soil classes from two tables as characters, concatenate them, turn
# them into factors and update levels

# Let's reset the subsets
obs.covr.sub1 <- obs.covr[-which((obs.covr$NZSC == "RXT") | (obs.covr$NZSC == "ROK")), ] # 105 observations
obs.covr.sub2 <- obs.covr[which((obs.covr$NZSC == "RXT") | (obs.covr$NZSC == "ROK")), ] # 94 observations

# Rationalise their factor levels

obs.covr.sub1$NZSC <- as.character(obs.covr.sub1$NZSC)
obs.covr.sub1$NZSC <- as.factor(obs.covr.sub1$NZSC)

obs.covr.sub2$NZSC <- as.character(obs.covr.sub2$NZSC)
obs.covr.sub2$NZSC <- as.factor(obs.covr.sub2$NZSC)

levels(obs.covr.sub1$NZSC)
# [1] "BFT" "BOT" "PIT" "RFT" "ROT"

levels(obs.covr.sub2$NZSC)
# [1] "ROK" "RXT"

soilclass.sub1 <- as.character(obs.covr.sub1$NZSC)
soilclass.sub2 <- as.character(obs.covr.sub2$NZSC)

soilclass.sub12 <- c(soilclass.sub1, soilclass.sub2)
soilclass.sub12.f <- as.factor(soilclass.sub12)

levels(soilclass.sub12.f)
# [1] "BFT" "BOT" "PIT" "RFT" "ROK" "ROT" "RXT"

# Pretend that sub1 and sub2 came from data that had more levels than are
# represented
lev1 <- c("BFT", "BOT", "PIT", "RFT", "ROT", "XXA", "XXB", "XXC")
lev2 <- c("ROK", "RXT", "YYA", "YYB", "YYC", "YYD", "YYE")
lev12 <- c(lev1, lev2)

# So it preserves the levels in the order that they were specified by levels(),
# not alphabetical order. That's interesting.
levels(soilclass.sub12.f) <- lev12
# [1] "BFT" "BOT" "PIT" "RFT" "ROT" "XXA" "XXB" "XXC" "ROK" "RXT" "YYA" "YYB" "YYC" "YYD" "YYE"

# Let's rbind the covariates then build a tree and examine its contingency table
covr.sub12 <- rbind(obs.covr.sub1[, 12:ncol(obs.covr.sub1)], obs.covr.sub2[, 12:ncol(obs.covr.sub2)])

tree.sub12.sf <- C5.0(covr.sub12, soilclass.sub12.f)
summary(tree.sub12.sf)

# (a)   (b)   (c)   (d)   (e)   (f)   (g)   (h)   (i)   (j)   (k)   (l)   (m)   (n)   (o)    <-classified as
# ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
#         1                                                                                  (a): class BFT
#        83           1     1     1     2                                                    (b): class BOT
#         2     6           1                                                                (c): class PIT
#         1           3                                                                      (d): class RFT
#                          34           3                                                    (e): class ROT
#               1                 2                                                          (f): class XXA
#        28     1     1     2          25                                                    (g): class XXB
#                                                                                            (h): class XXC
#                                                                                            (i): class ROK
#                                                                                            (j): class RXT
#                                                                                            (k): class YYA
#                                                                                            (l): class YYB
#                                                                                            (m): class YYC
#                                                                                            (n): class YYD
#                                                                                            (o): class YYE

# That's interesting. There are positive values in the contingency table for XXA
# and XXB even though there aren't actually any of these classes in the
# calibration dataset. Looks like this occurred because the concatenated factor
# levels in lev 12 were not sorted into alphabetical order:
lev12
# [1] "BFT" "BOT" "PIT" "RFT" "ROT" "XXA" "XXB" "XXC" "ROK" "RXT" "YYA" "YYB" "YYC" "YYD" "YYE"

# When they are unsorted and I amend the factor levels, it somehow alters some
# of the factor values:
levels(soilclass.sub12.f) <- lev12
soilclass.sub12.f

# [1] BOT BOT BOT BOT BOT RFT BOT BOT BOT BOT BOT BOT BOT XXA PIT PIT PIT XXA XXA RFT PIT PIT BOT PIT PIT PIT BOT BOT BOT BOT BOT BOT BOT PIT BOT BOT BOT
# [38] BOT BFT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT RFT RFT BOT BOT BOT BOT
# [75] BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT XXB XXB XXB XXB XXB XXB
# [112] XXB XXB ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT XXB XXB XXB XXB XXB XXB
# [149] XXB XXB XXB XXB XXB XXB XXB ROT XXB ROT XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB XXB
# [186] XXB XXB XXB XXB XXB XXB XXB XXB ROT ROT ROT ROT ROT ROT
# Levels: BFT BOT PIT RFT ROT XXA XXB XXC ROK RXT YYA YYB YYC YYD YYE

table(soilclass.sub12.f)
# soilclass.sub12.f
# BFT BOT PIT RFT ROT XXA XXB XXC ROK RXT YYA YYB YYC YYD YYE 
#   1  88   9   4  37   3  57   0   0   0   0   0   0   0   0 

# Notice that the tabulated values are still the same as for the original data.
# Compare the soilclass.sub12, which is a character vector of the original data.
table(soilclass.sub12)
# soilclass.sub12
# BFT BOT PIT RFT ROK ROT RXT 
#   1  88   9   4  37   3  57

# The problem appears to be solved if I sort the levels alphabetically before
# amending soilclass.sub12.f:
lev12 <- sort(lev12)
levels(soilclass.sub12.f) <- lev12
soilclass.sub12.f

# [1] BOT BOT BOT BOT BOT RFT BOT BOT BOT BOT BOT BOT BOT ROT PIT PIT PIT ROT ROT RFT PIT PIT BOT PIT PIT PIT BOT BOT BOT BOT BOT BOT BOT PIT BOT BOT BOT
# [38] BOT BFT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT RFT RFT BOT BOT BOT BOT
# [75] BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT RXT RXT RXT RXT RXT RXT
# [112] RXT RXT ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK RXT RXT RXT RXT RXT RXT
# [149] RXT RXT RXT RXT RXT RXT RXT ROK RXT ROK RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT
# [186] RXT RXT RXT RXT RXT RXT RXT RXT ROK ROK ROK ROK ROK ROK
# Levels: BFT BOT PIT RFT ROK ROT RXT XXA XXB XXC YYA YYB YYC YYD YYE

table(soilclass.sub12.f)
# soilclass.sub12.f
# BFT BOT PIT RFT ROK ROT RXT XXA XXB XXC YYA YYB YYC YYD YYE 
#   1  88   9   4  37   3  57   0   0   0   0   0   0   0   0

# But what if I add another class to the levels that doesn't sit at the end of
# the list?
lev12 <- c(lev12, "BBX")
lev12 <- sort(lev12)

rm(soilclass.sub12.f)
soilclass.sub12.f <- as.factor(soilclass.sub12)
levels(soilclass.sub12.f) <- lev12
soilclass.sub12.f

# Now some values are incorrectly substituted by BBX!

# [1] BFT BFT BFT BFT BFT PIT BFT BFT BFT BFT BFT BFT BFT ROK BOT BOT BOT ROK ROK PIT BOT BOT BFT BOT BOT BOT BFT BFT BFT BFT BFT BFT BFT BOT BFT BFT BFT
# [38] BFT BBX BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT PIT PIT BFT BFT BFT BFT
# [75] BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT BFT ROT ROT ROT ROT ROT ROT
# [112] ROT ROT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT RFT ROT ROT ROT ROT ROT ROT
# [149] ROT ROT ROT ROT ROT ROT ROT RFT ROT RFT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT ROT
# [186] ROT ROT ROT ROT ROT ROT ROT ROT RFT RFT RFT RFT RFT RFT
# Levels: BBX BFT BOT PIT RFT ROK ROT RXT XXA XXB XXC YYA YYB YYC YYD YYE

table(soilclass.sub12.f)
# soilclass.sub12.f
# BBX BFT BOT PIT RFT ROK ROT RXT XXA XXB XXC YYA YYB YYC YYD YYE 
#   1  88   9   4  37   3  57   0   0   0   0   0   0   0   0   0


# What if I create the factor differently?
soilclass.sub12.f.2 <- factor(soilclass.sub12, levels = lev12)
soilclass.sub12.f.2

# [1] BOT BOT BOT BOT BOT RFT BOT BOT BOT BOT BOT BOT BOT ROT PIT PIT PIT ROT ROT RFT PIT PIT BOT PIT PIT PIT BOT BOT BOT BOT BOT BOT BOT PIT BOT BOT BOT
# [38] BOT BFT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT RFT RFT BOT BOT BOT BOT
# [75] BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT BOT RXT RXT RXT RXT RXT RXT
# [112] RXT RXT ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK ROK RXT RXT RXT RXT RXT RXT
# [149] RXT RXT RXT RXT RXT RXT RXT ROK RXT ROK RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT RXT
# [186] RXT RXT RXT RXT RXT RXT RXT RXT ROK ROK ROK ROK ROK ROK
# Levels: BBX BFT BOT PIT RFT ROK ROT RXT XXA XXB XXC YYA YYB YYC YYD YYE

table(soilclass.sub12.f.2)
# soilclass.sub12.f.2
# BBX BFT BOT PIT RFT ROK ROT RXT XXA XXB XXC YYA YYB YYC YYD YYE 
#   0   1  88   9   4  37   3  57   0   0   0   0   0   0   0   0

# AHA! Creating the factor using factor() instead of as.factor() works properly, 
# because the levels are specified at the outset.