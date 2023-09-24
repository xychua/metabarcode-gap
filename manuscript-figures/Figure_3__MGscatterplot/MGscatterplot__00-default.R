## ****************************************************************************
## Setup and calling script for each of the other 5 primer pair
##
## __________
##   Author: Xin-Yi Chua
##
## ****************************************************************************



## ****************************************************************************
## check working directory

if (basename(getwd()) != "Figure_3__MGscatterplot") {
  stop("Please set the working directory to 'Figure_3__MGscatterplot'")
}



## ****************************************************************************
## load data ----
##
##    * load common plotting parameters e.g. colours, shapes, etc
##    * load data
##    * load helper functions
  
source('../utils/MGscatterplot-setup.R')
source('../utils/load-MGboundary.R')
source("../../data-processing/utils/plot_single_MGscatterplot.R")


## set which species are of interest, for the manuscript we were primarily
## interested in Fishes and Sharks
absLimits[, interest := query.class %in% c('Actinopteri', 
                                           'Chondrichthyes', 
                                           'Cladistia')]

source('MGscatterplot__01-Fish16S.R')
source('MGscatterplot__02-teleo.R')
source('MGscatterplot__03-MiFish.R')
source('MGscatterplot__04-COI-E.R')
source('MGscatterplot__05-COI-L.R')
