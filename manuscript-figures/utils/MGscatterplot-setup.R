## ****************************************************************************
## Common setup parameters for generating single MG scatter plot for each 
## primer pair.
##
## __________
##   Author: Xin-Yi Chua
##
## ****************************************************************************


## setup ----
library(data.table)
library(futile.logger)
invisible(flog.threshold(DEBUG))

library(ggrepel)
library(cowplot)
library(ggplot2)
theme_set(theme_light(base_size=12) +
            theme(axis.title.x = element_text(margin=margin(t=5, unit='mm')),
                  axis.title.y = element_text(margin=margin(l=0, r=5, unit='mm')),
                  panel.grid.minor = element_blank()))


TAXA_RANKS <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')


taxaShapes <- c(superkingdom = 16,
                     kingdom = 4,
                      phylum = 0,
                       class = 6,
                       order = 5,
                      family = 2,
                       genus = 1)

MGcolours <- c(overlap = 'orange2',
                    on = 'purple',
                   gap = 'deepskyblue3')


PRIMER_LABELS <- c(Fish16S = 'Fish16S',
                     teleo = 'Teleo',
                    MiFish = 'MiFish',
                     miniE = 'COI_E',
                     miniL = 'COI_L')

taxaColours <- c('black', RColorBrewer::brewer.pal(6,name='Dark2'))
names(taxaColours) <- TAXA_RANKS
