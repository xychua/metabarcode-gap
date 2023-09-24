## ****************************************************************************
## Case study: Hippoglossus hippoglossus detected by the COI_E primer pair
##
## Due to the extremely large size of the full pairwise alignments for COI_E
## (miniE) primer pair, we only extract the subset pairwise alignments relevant
## to the Hippoglossus hippoglossus case study.
##
## Descriptions of the data generation and code for data processing is provided
## in the code repository that describes the whole process.
##
## __________
##   Author: Xin-Yi Chua
##
## ****************************************************************************


## ****************************************************************************
## check working directory

if (basename(getwd()) != "Figure_5__case-study__Hippoglossus") {
  stop("Please set the working directory to 'Figure_5__case-study__Hippoglossus'")
}


## ****************************************************************************
## setup ----
##
##    * load common figure parameters
##    * specify query species
##    * load pairwise alignments

source("../utils/MGscatterplot-setup.R")

QRY <- 'Hippoglossus hippoglossus'

pwAlignments <- fread("data/Actinopteri-07_multi.tsv.gz")


## ****************************************************************************
## get plot data ----
##
##    * prepare data for histogram and heatmap
##    * reverse the taxonomy ranks for plotting

qry_pw <- pwAlignments[query.species == QRY]
qry_pw[, same.rank := factor(same.rank, 
                             levels = rev(TAXA_RANKS),
                             ordered = T)]


## specify which alignments will be used in the heatmap
qry_pw$plot.map <- F
qry_pw[same.rank <= 'genus' | 
         (same.rank >= 'family' & pident >= 95) |
         (same.rank >= 'order' & pident >= 90), 
       plot.map:=T]
qry_pw[query.species == QRY & target %in% qry_pw[plot.map == T, target], plot.map:=T]


## prepare heatmap data
hmap_data <- qry_pw[plot.map == T] 
hmap_data[,query:=as.character(query)]
hmap_data[,target:=as.character(target)]


## add self-matches to heatmap
## self-matches are 100% identical to themselves
self.matches <- unique(hmap_data[same.rank == 'species', 
                                 .(query, query.species,
                                   target=query, 
                                   target.species=query.species,
                                   same.taxa=query.species,
                                   pident=100, 
                                   coverage=100,
                                   query.variant.type,
                                   same.rank)])
hmap_data <- rbind(hmap_data, self.matches, use.names=T, fill=T)
rm(self.matches)

hmap_data[,query.acc:=tstrsplit(query, '_q', keep=1)]
hmap_data[,target.acc:=tstrsplit(target, '_q', keep=1)]
hmap_data[,ylbl:=paste(target.species, target.acc)]

## order rows by decreasing p-ident of the odd one out
odd_one_out <- hmap_data[same.rank == 'species', 
                   .SD[which.min(pident)], query.acc][,.N, target.acc][which.max(N), target.acc]
yorder <- hmap_data[query.acc == odd_one_out][order(pident), unique(ylbl)]
hmap_data[,ylbl:=factor(ylbl, levels=yorder)]

## order the columns
xorder <- hmap_data[query.acc == odd_one_out & same.rank == 'species'][order(-pident), target.acc]
hmap_data[,query.acc:=factor(query.acc, levels=xorder)]


## ****************************************************************************
## plot heatmap ----

p_heat <- ggplot(hmap_data, aes(x=query.acc, y=ylbl, fill=pident)) +
  geom_tile(col='white', linewidth=.1) +
  geom_text(aes(label=pident), size=3) +
  scale_x_discrete(position = 'top', 
                   name = 'Query',
                   expand = expansion(mult=c(0,0))) +
  scale_y_discrete(name = 'Target',
                   expand = expansion(mult=c(0,0))) +
  scale_fill_gradient(low='yellow', high='forestgreen', name='% identity') +
  facet_grid(same.rank ~ ., scales='free_y', space='free') +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size=8),
        axis.text.x = element_text(angle=90),
        legend.position = 'top', 
        legend.direction = 'horizontal',
        legend.key.height = unit(3, 'mm'),
        legend.key.width = unit(7, 'mm'),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        panel.spacing.y = unit(3, 'mm'),
        panel.background = element_rect(fill=NA, colour='grey20', linewidth=1),
        strip.background = element_rect(fill=NA, colour='grey20', linewidth=.3),
        strip.text.y.right = element_text(colour = 'black',
                                          angle = 0, 
                                          hjust = 0, 
                                          vjust = 1, 
                                          face = 'bold'))




## ****************************************************************************
## Plot histogram ----
##
##    * use the alpha channel to show which alignments in the histogram
##      correspond with the heatmap

p_hist <- ggplot(qry_pw, aes(x=pident, fill=same.rank, col=same.rank)) +
  geom_histogram(aes(alpha = plot.map), binwidth=1, linewidth=.1, closed='left', breaks=seq(0,101,1)) +
  stat_bin(geom="text", binwidth=1, breaks=seq(0,101,1), closed='left',
           aes(label = ifelse(after_stat(count) %between% c(1,5), after_stat(count), NA)), 
           size=3, col='black', vjust = - .25) +
  scale_x_continuous(limits = c(58,101)) +
  scale_y_continuous(expand = expansion(mult=c(0, .2))) +
  scale_fill_manual(values = taxaColours,
                    name = 'taxonomy rank',
                    aesthetics = c('fill','col'),
                    drop = F,
                    guide = 'none') +
  scale_alpha_discrete(name = 'shown in heatmap',
                       range = c(.4, .8),
                       drop = F) +
  facet_grid(same.rank ~ ., scales = 'free_y', switch = 'y') +
  theme(legend.position = 'top',
        legend.key.size = unit(4, 'mm'),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        strip.text.y.left = element_text(colour='black',
                                         face = 'bold',
                                         size = 10))



## ****************************************************************************
## combine plots ----

plt <- plot_grid(p_hist, p_heat, 
                 nrow = 1, 
                 align = 'h',
                 axis = 'tb',
                 labels = 'AUTO', 
                 rel_widths = c(1,1.5))
ggsave(plt, 
       filename = 'figures/example__Hippoglossus.png', 
       bg = 'white',
       height = 7, 
       width = 9)

