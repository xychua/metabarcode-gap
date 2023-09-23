## ****************************************************************************
## Plot pairwise histogram
##
## Author: Xin-Yi Chua (x.chua@connect.qut.edu.au)
##
## This script generates a histogram of pairwise alignments for a SINGLE
## query species of interest, from data of a single primer pair.
##
## The plot is further segregated by the taxonomic ranks of the lowest common
## taxonomy between each query/target pair of sequences.
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
OUTPUT_FIG <- sprintf("%s_MG-histogram.png", PRIMER_PAIR)



## plotting parameters
taxaColours <- c('black', RColorBrewer::brewer.pal(6,name='Dark2'))
names(taxaColours) <- TAXA_RANKS

FIG_HEIGHT <- 10
FIG_WIDTH <- 5.5


## ****************************************************************************
## load data ----
##
##    * specify a query species of interest
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
## Histogram -----
##
##  * drop factor levels, only use the 7 taxonomy ranks
##  * reverse taxonomy rank order for facet_grid

interest_pwAlign[, same.rank := droplevels(same.rank)]
interest_pwAlign[, same.rank := factor(same.rank, levels=rev(TAXA_RANKS))]

p_hist <- ggplot(interest_pwAlign, aes(x=pident, fill=same.rank, col=same.rank)) +
  geom_histogram(binwidth=1, linewidth=.1, alpha=.6) +
  stat_bin(binwidth=1, geom="text", 
           aes(label = ifelse(after_stat(count) %between% c(1,5), after_stat(count), NA)), 
           size = 3, 
           col = "black",
           vjust = -0.25) +
  scale_x_continuous(limits = c(58, 100), name='% identity') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                     expand = expansion(mult=c(0,.15))) + 
  scale_fill_manual(values = taxaColours,
                    name = 'taxonomy rank',
                    aesthetics = c("fill", "col"), 
                    drop = F,
                    guide = "none") +
  labs(title = "Distribution of pairwise alignments for",
       subtitle = SPECIES_OF_INTEREST,
       y = 'Num alignments') +
  facet_grid(same.rank ~ ., 
             drop = F,
             scales = "free_y",
             switch = "y") +
  theme_minimal_hgrid(font_size = 12) +
  theme(plot.title = element_text(size = 12, face = "plain"),
        plot.subtitle = element_text(size = 14, face = "bold.italic"),
        panel.spacing.y = unit(5, 'mm'),
        panel.border = element_rect(colour='grey30'),
        strip.background = element_rect(fill='grey85'))



## save figure
flog.info("Plotting: %s", OUTPUT_FIG)
ggsave(p_hist,
       filename = OUTPUT_FIG,
       bg = "white",
       height = FIG_HEIGHT,
       width = FIG_WIDTH)

