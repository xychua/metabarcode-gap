## ****************************************************************************
## Plot MG scatter plot for single primer pair
##
## Author: Xin-Yi Chua (x.chua@connect.qut.edu.au)
##
## This script generates the metabarcoding gap (MG) scatterplot for one primer
## pair data set.
##
## Other example figures can be found in 'manuscript-figures' folder.
##
## !! PLEASE NOTE !!
## This script is to be run from R and not the command line as file paths 
## are hard coded in the script and needs to change accordingly.
##
## ****************************************************************************




## ****************************************************************************
## setup ----

library(data.table)
library(futile.logger)

library(ggrepel)
library(cowplot)
library(ggplot2)
theme_set(theme_light(base_size=12) +
            theme(axis.title.x = element_text(margin=margin(t=5, unit='mm')),
                  axis.title.y = element_text(margin=margin(l=0, r=5, unit='mm')),
                  panel.grid.minor = element_blank()))

## helper functions for plotting
source("utils/plot_single_MGscatterplot.R")


## Plotting parameters
taxaShapes <- c(superkingdom = 16,
                     kingdom = 4,
                      phylum = 0,
                       class = 6,
                       order = 5,
                      family = 2,
                       genus = 1)

MGcolours <- c(overlap = "orange2",
                    on = "purple",
                   gap = "deepskyblue3")


## Saving figure parameters
curr.primer <- 'teleo'
single_fig.height <- 8
single_fig.width <- 9.75
facet_fig.height <- 10
facet_fig.width <- 15.5



## ****************************************************************************
## load data ----
##
##    * users can then define which species are of interest for further
##      highlighting in the plots


abs.limits <- readRDS("MGboundary.rds")
abs.limits[, interest:=query.class %in% c('Actinopteri', 'Chondrichthyes', 'Cladistia')]



## ****************************************************************************
## base plot ----

flog.info("Plotting single ...")
if (hasName(abs.limits, 'show.lbl')) {
  abs.limits$show.lbl <- NULL
}
abs.limits[interest == T & intra.minP < 70, show.lbl:=query.species]
abs.limits[interest == T & intra.minP < 100 & inter.maxP < 80, show.lbl:=query.species]

pbase <- plotMGscatter(abs.limits) + labs(title = curr.primer)

## add labels
plt <- add_labels(pbase, abs.limits, wrap=F)


figFile <- sprintf("%s_MG-scatterplot.png", curr.primer)
flog.info("Plotting: %s", figFile)
ggsave(plt, 
       file = figFile, 
       bg = "white", 
       height = single_fig.height, 
       width = single_fig.width)



## ****************************************************************************
## split by rank ----
##
##    * split the scatterplot by taxonomic ranks
##    * readjust labels


flog.info("Plotting by ranks")
if (hasName(abs.limits, 'show.lbl')) {
  abs.limits$show.lbl <- NULL
}
abs.limits[interest == T & intra.minP < 70, show.lbl:=query.species]
abs.limits[interest == T & intra.minP < 100 & inter.maxP < 80, show.lbl:=query.species]

plt <- add_labels(pbase, abs.limits) + 
  facet_wrap(inter.rank ~ ., nrow=2) +
  labs(title = curr.primer) +
  theme(strip.text = element_text(size=14, face='bold', colour='grey30'),
        strip.background = element_rect(fill='grey85'),
        panel.background = element_rect(fill=NA, colour='grey30', linewidth=.3))


figFile <- sprintf("%s_MG-scatterplot_rank.png", curr.primer)
flot.info("Plotting: %s", figFile)
ggsave(plt, 
       file = figFile,
       bg = "white",
       height = facet_fig.height, 
       width = facet_fig.width)
