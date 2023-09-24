## ****************************************************************************
## Case study: taxa below the diagonal
##
## Histograms for example taxa below the diagonal
##
## Due to the large number of pairwise alignments (and size), only the relevant
## data for example species is provided in this section.
##
## Format of tables are explained in the data-processing documentation
##
## __________
##   Author: Xin-Yi Chua
##
## ****************************************************************************


## ****************************************************************************
## check working directory

if (basename(getwd()) != "Figure_6__case-study__Below-diagonal") {
  stop("Please set the working directory to 'Figure_6__case-study__Below-diagonal'")
}



## ****************************************************************************
## setup ----

library(data.table)
library(futile.logger)

library(cowplot)
library(ggplot2)
theme_set(theme_light(base_size = 14) +
            theme(axis.title.x = element_text(size=12, margin=margin(t=5, unit='mm')),
                  axis.title.y = element_text(size=12, margin=margin(l=0, r=5, unit='mm')),
                  panel.border = element_rect(colour='black'),
                  panel.grid.minor = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.spacing.y = unit(5, 'mm'),
                  plot.title = element_text(size=12, hjust=-.1),
                  strip.background = element_rect(fill='grey85'),
                  strip.text = element_text(size=12, face='bold', colour='grey10')))

## common plotting parameters, e.g., colours, shapes etc
source("../utils/MGscatterplot-setup.R")
LOGLABELS <- function(x) ifelse(x >= 5e5, 
                                paste0(round(x/1e6,2),'M'),
                                ifelse(x >= 1e3,paste0(round(x/1e3),'K'),x))


## define case studies and data folder where annotated pairwise alignments are
CASE.STUDIES <- c("Herklotsichthys quadrimaculatus",
                  "Pantodon buchholzi",
                  "Schindleria sp. 2 TK-2007")

ANNOT_DIR <- "../data/Fish16S_annot"



## ****************************************************************************
## prepare data -----
##
##    * see the data-processing folder for documentation
##    * usually we will only load in the relevant pairwise alignments for
##      the selected species of interest
          
pwFiles <- list.files('data', pattern = ".tsv.gz", full.names = T)
pwAlign <- rbindlist(lapply(pwFiles, function(fn) {
  flog.info("Loading: %s", fn)
  dat <- fread(fn)
  dat[query.species %in% CASE.STUDIES]
}), use.names = T, fill = T)

flog.info("Loaded: %d alignments belonging to %d species",
          pwAlign[,.N],
          pwAlign[, uniqueN(query.species)])

## reverse taxonomy rank for plotting
pwAlign[,same.rank:=factor(same.rank, levels = rev(TAXA_RANKS))]


## ****************************************************************************
## plot histogram -----
##
##    * PNG files are saved in the figures/ sub-folder
##    * a ggplot2 object is returned for the combined plot below


all.plots <- lapply(CASE.STUDIES, function(qry.species) {
  ss.align <- pwAlign[!is.na(same.rank) & query.species == qry.species]
  
  curr.primer <- ss.align[,unique(primer)]
  if (length(curr.primer) > 1) {
    ss.align[,primer:=factor(primer,
                             levels=names(PRIMER_LABELS),
                             labels=PRIMER_LABELS)]
    fig.width <- 12
    PREFIX <- 'all'
  } else {
    fig.width <- 4.5
    PREFIX <- curr.primer
  }
  
  plt <- ggplot(ss.align, aes(x=pident, fill=same.rank)) +
    geom_histogram(alpha=.7, binwidth=1, linewidth=.2, colour='grey90') +
    scale_fill_manual(values=taxaColours) +
    scale_x_continuous(breaks=seq(0,100,10), name='% identity') +
    scale_y_continuous(labels=LOGLABELS, 
                       name = 'Num alignments',
                       expand=expansion(mult=c(0,.1))) +
    labs(subtitle = qry.species) + 
    facet_grid(same.rank ~ primer, scales='free_y', switch='y', drop=F) +
    theme(legend.position='none',
          axis.line = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.border = element_rect(colour='grey30', linewidth=.3),
          panel.spacing.y=unit(0, 'mm'))
  
  outFile <- sprintf('figures/%s-%s.png', PREFIX, qry.species)
  flog.info("Plotting: %s", outFile)
  ggsave(plt, bg='white', height=8, width=fig.width, filename= outFile)
  
  plt
})



## ****************************************************************************
## combine histograms and MG scatter-plot ----
##
##    * load in the ggplot2 scatter plot object from previous example in
##      Figure_3__MGscatterplot folder
##    * combine the scatterplot with the histograms generated from this script

tweak.plots <- lapply(all.plots, function(plt) {
  plt + theme(panel.spacing.y = unit(0, 'mm'),
              panel.border = element_rect(linewidth=.2),
              plot.subtitle = element_blank(),
              plot.margin = margin(10,20,1,1, unit='mm'))
})

## re-position the legends
inFile <- '../Figure_3__MGscatterplot/figures/Fish16S.rds'
if (!file.exists(inFile)) {
  flog.error("File not found: %s", inFile)
  stop("Run the sripts in Figure__3__MGscatterplot first")
}
mg_scatter <- readRDS('../Figure_3__MGscatterplot/figures/Fish16S.rds')
mg_tweaked <- mg_scatter +
  guides(shape = guide_legend(ncol = 2, order=2),
         size = guide_legend(direction = 'horizontal', 
                             title.position = 'top',
                             label.position = 'bottom',
                             override.aes = list(shape=1), 
                             order=3)) +
  theme(plot.margin = margin(t=2, r=2, l=2, b=10, unit='mm'))

combined <- plot_grid(mg_tweaked,
                      plot_grid(plotlist = tweak.plots,
                                nrow=1, 
                                align='h',
                                axis='tb',
                                labels = CASE.STUDIES,
                                label_x = 0.2,
                                label_size = 14,
                                label_fontface = 'plain',
                                hjust = 0),
                      nrow=2)

figFile <- 'figures/Fish16S-combined.png'
flog.info("Plotting: %s", figFile)
ggsave(combined, 
       filename = figFile,
       bg='white', 
       height=12, 
       width=11)
