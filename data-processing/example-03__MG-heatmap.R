## ****************************************************************************
## Plot pairwise histogram
##
## Author: Xin-Yi Chua (x.chua@connect.qut.edu.au)
##
## This script generates a heatmap of key pairwise alignments for a SINGLE
## query species of interest, from data of a single primer pair.
##
## The heatmap generated display all sequences belonging to the species of 
## interest as columns and the target sequences belonging to all other species
## as rows. As pairwise alignments can be very large data sets, generating
## a big heatmap is also not useful. The script will select only the top targets
## per taxonomic rank for the visualisation as the goal is to find out which
## other species are the closest match.
##
## Other example figures can be found in 'manuscript-figures' folder.
##
## !! PLEASE NOTE !! 
## This script is to be run from R and not the command line as file paths 
## are hard coded in the script and needs to change accordingly.
##
##
## ****************************************************************************




## ****************************************************************************
## setup ----

library(data.table)
library(futile.logger)
library(forcats)

library(cowplot)
library(ggplot2)

TAXA_RANKS <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')



## file parameters
INPUT_METADATA <- "03-taxaAnnot/metadata-variant-type.rds"
DIR_ALIGNMENTS <- "03-taxaAnnot"
PRIMER_PAIR <- "teleo"
OUTPUT_FIG <- sprintf("%s_MG-heatmap.png", PRIMER_PAIR)



## plotting parameters
taxaColours <- c('black', RColorBrewer::brewer.pal(6,name='Dark2'))
names(taxaColours) <- TAXA_RANKS

FIG_HEIGHT <- 10
FIG_WIDTH <- 5.5


## ****************************************************************************
## load data ----
##
##    * specify a query species of interest e.g., Bagarius yarrelli
##    * find the group files that contain this species of interest
##    * load in the pairwise alignments related to this species only
##    * concatenate the pairwise alignments


SPECIES_OF_INTEREST <- "Bagarius yarrelli"


flog.info("Loading metadata: %s", INPUT_METADATA)
metadata <- readRDS(INPUT_METADATA)

pwFiles <- metadata[species == SPECIES_OF_INTEREST]
pwFiles[,rdsFile := sprintf("%s/%s_%s.rds",  
                            DIR_ALIGNMENTS, 
                            group,
                            variant.type)]
pwFiles[,file.exists := file.exists(rdsFile)]

if (pwFiles[file.exists == F, .N] > 0) {
  flog.error("Taxonomically annotated files from Step 3 not found",
             pwFiles[file.exists == F, unique(rdsFile)],
             capture = T)
  stop()
}


pwFiles <- pwFiles[file.exists == T, unique(rdsFile)]
flog.info("Num alignments files to merge: %d", length(pwFiles))

if (length(pwFiles) == 0) {
  flog.error("There are no files to process, stopped.")
  quit()
}


## concatenate all tables
pwAlignments <- rbindlist(lapply(pwFiles, function(fn) { 
  flog.debug("Parsing: %s", fn)
  readRDS(fn)
}), use.names = T, fill = T)


## extract pairwise alignments for species_of_interest
interest_pwAlign <- pwAlignments[query.species == SPECIES_OF_INTEREST]
if (interest_pwAlign[, .N] == 0) {
  flog.error("No alignments for [%s] found", SPECIES_OF_INTEREST)
  stop()
}

flog.info("Loaded: %d records belonging to %d query species",
          interest_pwAlign[,.N],
          interest_pwAlign[,uniqueN(query.species)])



## ****************************************************************************
## prepare heatmap data ----
##
##    * columns are sequences belonging to the species of interest (query)
##    * rows are target sequences to which the query species is aligned
##    * heatmaps can be very large, so we limit the number of target rows by
##      selecting top N matches per taxonomic rank

topN <- 10
topN_alignments <- interest_pwAlign[order(-pident), .SD[1:10], by=same.rank]

hmap_data <- interest_pwAlign[same.rank == 'species' |
                                target %in% topN_alignments$target]
hmap_data[, query := as.character(query)]
hmap_data[, target := as.character(target)]

## add self-matches to heatmap data.table, self-matches will be 100% similar with
## 100% coverage, the same target.species and same.taxa
self_matches <- unique(interest_pwAlign[, .(query, 
                                            query.species,
                                            target = query, 
                                            target.species = query.species,
                                            same.taxa = query.species,
                                            pident = 100, 
                                            coverage = 100,
                                            query.variant.type,
                                            same.rank = 'species')])
hmap_data <- rbind(hmap_data, self_matches, use.names=T, fill=T)
rm(self_matches)

hmap_data[, query.acc := tstrsplit(query, '_q', keep=1)]
hmap_data[, target.acc := tstrsplit(target, '_q', keep=1)]
hmap_data[, ylbl := paste(target.species, target.acc)]

## order rows by decreasing p-ident of the odd one out
odd_one_out <- hmap_data[same.rank == 'species', 
                         .SD[which.min(pident)], query.acc][,.N, target.acc][which.max(N), target.acc]
yorder <- hmap_data[query.acc == odd_one_out][order(pident), unique(ylbl)]
hmap_data[, ylbl := factor(ylbl, levels=yorder)]
xorder <- c(odd_one_out, setdiff(hmap_data[,unique(query.acc)], odd_one_out))
hmap_data[, query.acc := factor(query.acc, levels = xorder)]

hmap_data[, same.rank := droplevels(same.rank)]
hmap_data[, same.rank := factor(same.rank, levels=rev(TAXA_RANKS))]



## ****************************************************************************
## Heatmap ----
##


## find the gap in the histogram for colouring fonts
pcuts <- rev(table(cut(hmap_data$pident, 
                       breaks = seq(0,100,5),
                       labels = seq(5,100,5))))
MIDPOINT <- as.integer(names(which(pcuts == 0)[1]))

p_heat <- ggplot(hmap_data, aes(x=query.acc, y=ylbl)) +
  geom_tile(aes(fill = pident), col = 'white', linewidth=.1) +
  geom_text(aes(label = pident, col = pident < MIDPOINT), size=2.5) +
  scale_x_discrete(position = 'top', 
                   name = '',
                   expand = expansion(mult=c(0,0))) +
  scale_y_discrete(position = 'right',
                   name = '',
                   expand = expansion(mult=c(0,0))) +
  scale_fill_gradient(low='yellow', 
                      high='darkgreen', 
                      name='% identity') +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "white"),
                      guide = 'none') +
  facet_grid(same.rank ~ ., 
             scales = 'free_y',
             space = 'free',
             switch = 'y') +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size=7),
        axis.text.x = element_text(angle=30, hjust=.09),
        legend.position = 'top',
        legend.direction = 'horizontal',
        legend.key.height = unit(3, 'mm'),
        legend.key.width = unit(7, 'mm'),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        panel.spacing.y = unit(3, 'mm'),
        panel.background = element_rect(fill=NA, colour='grey20', linewidth=1),
        strip.background = element_rect(fill=NA, colour='grey20', linewidth=.3),
        strip.text.y.left = element_text(angle=0, hjust=0, vjust=1, face='bold'))


## save figure
flog.info("Plotting: %s", OUTPUT_FIG)
ggsave(p_heat,
       filename = OUTPUT_FIG,
       bg = "white",
       height = FIG_HEIGHT,
       width = FIG_WIDTH)

