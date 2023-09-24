## ****************************************************************************
## Case study: Comparing multiple primer pairs using Epinephelus example 
##
## __________
##   Author: Xin-Yi Chua
##
## ****************************************************************************


## ****************************************************************************
## check working directory

if (basename(getwd()) != "Figure_4__case-study__Epinephelus") {
  stop("Please set the working directory to 'Figure_4__case-study__Epinephelus'")
}


## ****************************************************************************
## setup ----

source("../utils/MGscatterplot-setup.R")
source("../utils/plot_MGscatter.R")


## load data
# source("../utils/load-insilico-amplicons.R")
source("../utils/load-MGboundary.R")


QRY <- 'Epinephelus'



## ****************************************************************************
## venn diagram ----

interest_group <- absLimits[query.genus == QRY]
interest_group[, num_primers := uniqueN(primer), query.species]
interest_group[order(num_primers), .(nSpecies=uniqueN(query.species)), num_primers]

tally <- interest_group[, .N, .(primer, query.species)]
groups <- split(tally$query.species, tally$primer)

png(filename='figures/venn.png', res = 300, height=6, width=6, units='in')
venn::venn(groups, box=F, ilcs=2, sncs=2, zcolor=c('green','red','orange','blue','cyan'))
dev.off()



## ****************************************************************************
## MG-limits ----
##
##    * limit to only the Fish16S primer pair

absLimits$interest <- F
absLimits[query.genus == QRY, interest := T]
dcast(absLimits,
      primer ~ interest,
      value.var = 'query.species',
      fun.aggregate = uniqueN)

absLimits[interest == T & primer == 'Fish16S', uniqueN(query.species), MG.type]


## ****************************************************************************
## MG scatter plot ----
##
##    reduced version for main manuscript


## shorten species labels, abbreviate genus name to E.
absLimits$show.lbl <- NULL
if (QRY == 'Epinephelus') {
  absLimits[query.species %in% interest_group[num_primers >= 4, query.species],
            show.lbl:=gsub(QRY, "E.", query.species)]
} else {
  absLimits[query.species %in% interest_group[num_primers >= 4, query.species], show.lbl:=query.species]
}


## figure formatting style
absLimits[!is.na(show.lbl), show.lbl:=paste0("italic(\"", show.lbl, "\")")]
absLimits[!is.na(show.lbl), show.lbl:=gsub(' x ', ' &\n', show.lbl)]

dummy <- data.table(primer = 'NA')
dat <- rbind(dummy, absLimits[intra.minP >= 80 & inter.maxP >= 80], use.names=T, fill=T)

plt <- plotMGscatter(dat, show.lbl = T, lbl.size=3)

figFile <- sprintf('figures/example__%s-1.png', QRY)
flog.info("Plotting: %s", figFile)
ggsave(plt, filename = figFile, bg='white', height=20, width=16)




## ****************************************************************************
## Supplementary figures ----
##
##    * show all labels

absLimits$show.lbl <- NULL
absLimits[query.species %in% interest_group$query.species, 
          show.lbl:=gsub(QRY, "E.", query.species)]
absLimits[!is.na(show.lbl), show.lbl:=paste0("italic('", show.lbl, "')")]

plt <- plotMGscatter(absLimits[primer == 'Fish16S'], show.lbl = T, lbl.size=5, nudge_x = -5, nudge_y=-15)
ggsave(plt, filename=sprintf('figures/example__%s-Fish16S.png', QRY), bg='white', height=13, width=18)


plt <- plotMGscatter(absLimits[primer == 'teleo'], show.lbl = T, lbl.size=5, nudge_x = -5)
ggsave(plt, filename=sprintf('figures/example__%s-teleo.png', QRY), bg='white', height=13, width=16)


plt <- plotMGscatter(absLimits[primer == 'MiFish'], show.lbl = T, lbl.size=5)
ggsave(plt, filename=sprintf('figures/example__%s-MiFish.png', QRY), bg='white', height=13, width=16)


plt <- plotMGscatter(absLimits[primer == 'miniE'], show.lbl = T, lbl.size=5)
ggsave(plt, filename=sprintf('figures/example__%s-COI_E.png', QRY), bg='white', height=13, width=14)


plt <- plotMGscatter(absLimits[primer == 'miniL'], show.lbl = T, lbl.size=5)
ggsave(plt, filename=sprintf('figures/example__%s-COI_L.png', QRY), bg='white', height=11, width=16)



## ****************************************************************************
## Supp: full figures combined ----
##

dat <- rbind(dummy, absLimits[intra.minP >= 80 & inter.maxP >= 80], use.names=T, fill=T)
plt <- plotMGscatter(dat, show.lbl = T, lbl.size=4)
figFile <- sprintf('figures/example__%s-2.png', QRY)
flog.info("Plotting: %s", figFile)
ggsave(plt, filename=figFile, bg='white', height=20, width=16)

