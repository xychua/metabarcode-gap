## ****************************************************************************
## Load metadata data for all primer pairs and combine into one table
##
## __________
##   Author: Xin-Yi Chua
##
## ****************************************************************************

require(futile.logger)
require(data.table)

PRIMER_LABELS <- c(Fish16S = 'Fish16S',
                     teleo = 'Teleo',
                    MiFish = 'MiFish',
                     miniE = 'COI_E',
                     miniL = 'COI_L')

metadata <- rbindlist(lapply(names(PRIMER_LABELS), function(pr) {
  flog.info("Parsing: %s", pr)
  fread(sprintf('../data/nt.201905__%s__taxaMetadata.tsv', pr))
}), use.names=T, fill=T)

## set factor labels
metadata[,primer := factor(primer, levels=names(PRIMER_LABELS))]

flog.info("MG limits loaded: ", metadata[,.N, primer], capture = T)

