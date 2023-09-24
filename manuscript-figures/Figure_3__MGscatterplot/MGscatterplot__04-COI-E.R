## ****************************************************************************
## Do not call this script directly, use the MGscatterplot__00-default.R
##
## Plot MG scatter plot for COI_E primer pair
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

curr.primer <- 'miniE'
single_fig.height <- 8
single_fig.width <- 8
facet_fig.height <- 10
facet_fig.width <- 13


flog.info('Subsetting data ...')
abs.limits <- absLimits[primer == curr.primer]

## species in other case studies
CASE.STUDIES <- c("Isurus oxyrinchus")



## ****************************************************************************
## base plot ----

flog.info("Plotting single ...")
if (hasName(abs.limits, 'show.lbl')) {
  abs.limits$show.lbl <- NULL
}
abs.limits[interest == T & intra.minP < 78, show.lbl:=query.species]
abs.limits[query.species %like% 'Johnius', show.lbl:=query.species]
abs.limits[query.species %in% CASE.STUDIES, show.lbl:=query.species]

## simple version
pbase <- plotMGscatter(abs.limits) + labs(title = PRIMER_LABELS[[curr.primer]])

## add labels
tmp <- abs.limits[!is.na(show.lbl), .(show.lbl, interest,
                                      inter.rank, nVariants,
                                      intra.minP, inter.maxP, 
                                      y=c(78, 102, 98, 95))]
plt <-  pbase + 
  geom_text_repel(data = tmp, 
                  aes(label = show.lbl),
                  colour = 'grey10',
                  size = 4.5,
                  box.padding = .1,
                  hjust = 1,
                  nudge_y = tmp$y - tmp$inter.maxP) +
  theme(axis.title.x = element_text(hjust=0, size=14),
        axis.title.y = element_text(hjust=0, size=14))

figFile <- sprintf('figures/%s-SIMPLE.png', PRIMER_LABELS[[curr.primer]])
flog.info("\tsaving plot: %s", figFile)
ggsave(plt, file=figFile, bg='white', height=7, width=7.25)


cropped_plt <- plt + scale_y_continuous(limits=c(90, 105), name='Next closest species') +
  theme(legend.position = 'top',
        plot.margin = margin(t=150, r=3, b=3, l=3),
        plot.title = element_text(hjust=0, margin=margin(b=-150)))
ggsave(cropped_plt, file=sprintf('figures/%s-CROPPED.png', PRIMER_LABELS[[curr.primer]]),
       bg='white', height=7, width=7.25)


## more labels
abs.limits[interest == T & MG.type == 'overlap' & inter.maxP < 90, show.lbl:=query.species]
abs.limits[interest == T & MG.type == 'gap' & 
             intra.minP < 90 & 
             inter.maxP < 90, show.lbl:=query.species]

pbase <- plotMGscatter(abs.limits) +
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
abs.limits[interest == T & intra.minP < 80, show.lbl:=query.species]
abs.limits[interest == T & intra.minP < 95 & inter.maxP < 90, show.lbl:=query.species]
abs.limits[interest == T & inter.rank == 'kingdom' & MG.type == 'gap', show.lbl:=query.species]
abs.limits[interest == T & inter.rank == 'phylum' & inter.maxP < 100 & intra.minP < 99, show.lbl:=query.species]
abs.limits[interest == T & inter.rank == 'class' & inter.maxP == 100 & intra.minP < 83, show.lbl:=query.species]
abs.limits[interest == T & inter.rank == 'order' & inter.maxP == 100 & intra.minP < 83, show.lbl:=query.species]
abs.limits[interest == T & inter.rank == 'family' & inter.maxP == 100 & intra.minP < 85, show.lbl:=query.species]

p <- add_labels(pbase, abs.limits) + 
  facet_wrap(inter.rank ~ ., nrow=2) +
  labs(title = PRIMER_LABELS[[curr.primer]]) +
  theme(strip.text = element_text(size=14, face='bold', colour='grey30'),
        strip.background = element_rect(fill='grey85'),
        panel.background = element_rect(fill=NA, colour='grey30', linewidth=.3))

figFile <- sprintf('figures/%s-by-ranks.png', PRIMER_LABELS[[curr.primer]])
flog.info("\tsaving plot: %s", figFile)
ggsave(p, file=figFile, bg='white', height=facet_fig.height, width=facet_fig.width)
