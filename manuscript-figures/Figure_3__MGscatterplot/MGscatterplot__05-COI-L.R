## ****************************************************************************
## Do not call this script directly, use the MGscatterplot__00-default.R
##
## Plot MG scatter plot for COI_L primer pair
##
## __________
##   Author: Xin-Yi Chua
##
## ****************************************************************************


if (any(!sapply(c('absLimits', 'PRIMER_LABELS'), exists))) {
  stop("Don't call this script directly, use MGscatterplot__00-default.R")
}


## ****************************************************************************
## setup ----

curr.primer <- 'miniL'
single_fig.height <- 8
single_fig.width <- 11.5
facet_fig.height <- 13
facet_fig.width <- 16


flog.info('Subsetting data ...')
abs.limits <- absLimits[primer == curr.primer]



## ****************************************************************************
## base plot ----

flog.info("Plotting single ...")
if (hasName(abs.limits, 'show.lbl')) {
  abs.limits$show.lbl <- NULL
}
abs.limits[interest == T & intra.minP < 80, show.lbl:=query.species]
abs.limits[interest == T & intra.minP < 83 & inter.maxP < 99, show.lbl:=query.species]
abs.limits[interest == T & intra.minP < 91 & inter.maxP < 90, show.lbl:=query.species]

# case.studies <- c('Butis butis',
#                   'Perccottus glenii')
# qry <- abs.limits[query.species %in% qry]
# qry[,show.lbl:=query.species]


pbase <- plotMGscatter(abs.limits, not.interest.col = 'grey75') +
  labs(title = PRIMER_LABELS[[curr.primer]])

## add labels
plt <- add_labels(pbase, abs.limits, wrap=F)

figFile <- sprintf('figures/%s.png', PRIMER_LABELS[[curr.primer]])
flog.info("\tsaving plot: %s", figFile)
ggsave(plt, file=figFile, bg='white', height=single_fig.height, width=single_fig.width)




## ****************************************************************************
## split by rank ----
##    readjust labels

flog.info("Plotting by ranks")
if (hasName(abs.limits, 'show.lbl')) {
  abs.limits$show.lbl <- NULL
}
abs.limits[interest == T & inter.rank == 'class' & intra.minP < 85, show.lbl:=query.species]
abs.limits[interest == T & inter.rank == 'order' & intra.minP < 85, show.lbl:=query.species]
abs.limits[interest == T & inter.rank == 'family' & 
             intra.minP < 90 & inter.maxP < 99, 
           show.lbl:=query.species]

p <- add_labels(pbase, abs.limits) + 
  facet_wrap(inter.rank ~ ., nrow=3) +
  labs(title = PRIMER_LABELS[[curr.primer]]) +
  theme(strip.text = element_text(size=14, face='bold', colour='grey30'),
        strip.background = element_rect(fill='grey85'),
        panel.background = element_rect(fill=NA, colour='grey30', linewidth=.3))

figFile <- sprintf('figures/%s-by-ranks.png', PRIMER_LABELS[[curr.primer]])
flog.info("\tsaving plot: %s", figFile)
ggsave(p, file=figFile, bg='white', height=facet_fig.height, width=facet_fig.width)
