## ****************************************************************************
## Do not call this script directly, use the MGscatterplot__00-default.R
##
## Plot MG scatter plot for Fish16S primer pair
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

curr.primer <- 'Fish16S'
single_fig.height <- 8
single_fig.width <- 10.75
facet_fig.height <- 10
facet_fig.width <- 18


flog.info('Subsetting data ...')
abs.limits <- absLimits[primer == curr.primer]

## species in other case studies
CASE.STUDIES <- c("Herklotsichthys quadrimaculatus",
                  "Pantodon buchholzi",
                  "Schindleria sp. 2 TK-2007")


## ****************************************************************************
## base plot ----

flog.info("Plotting single ...")
if (hasName(abs.limits, 'show.lbl')) {
  abs.limits$show.lbl <- NULL
}
abs.limits[interest == T & intra.minP %between% c(50, 65) &
             inter.maxP %between% c(85,100), 
           show.lbl:=query.species]

## plot
pbase <- plotMGscatter(abs.limits) + labs(title = PRIMER_LABELS[[curr.primer]])

## for other case studies
case.studies <- abs.limits[query.species %in% CASE.STUDIES]
case.studies[,show.lbl:=query.species]

plt <- add_labels(pbase, abs.limits, wrap=F) +
  geom_text_repel(data = case.studies, 
                  aes(label = show.lbl),
                  colour = 'grey10',
                  size = 4,
                  hjust = 1,
                  vjust = 1, 
                  nudge_x = -1,
                  nudge_y = -2)

figFile <- sprintf('figures/%s.png', PRIMER_LABELS[[curr.primer]])
flog.info("\tsaving plot: %s", figFile)
ggsave(plt, file=figFile, bg='white', height=single_fig.height, width=single_fig.width)

saveRDS(plt, file=gsub('.png', '.rds', figFile))



## ****************************************************************************
## split by rank ----
##    readjust labels


flog.info("Plotting by ranks")
if (hasName(abs.limits, 'show.lbl')) {
  abs.limits$show.lbl <- NULL
}
abs.limits[interest == T & intra.minP < 70, show.lbl:=query.species]
abs.limits[interest == T & intra.minP < 98 & inter.maxP < 78, show.lbl:=query.species]

p <- add_labels(pbase, abs.limits) + 
  facet_wrap(inter.rank ~ ., nrow=2) +
  labs(title = PRIMER_LABELS[[curr.primer]]) +
  theme(strip.text = element_text(size=14, face='bold', colour='grey30'),
        strip.background = element_rect(fill='grey85'),
        panel.background = element_rect(fill=NA, colour='grey30', linewidth=.3))

figFile <- sprintf('figures/%s-by-ranks.png', PRIMER_LABELS[[curr.primer]])
flog.info("\tsaving plot: %s", figFile)
ggsave(p, file=figFile, bg='white', height=facet_fig.height, width=facet_fig.width)
