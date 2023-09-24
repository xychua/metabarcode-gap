## *****************************************************************************
## Summarise MGboundary
##
## Author: Xin-Yi Chua (x.chua@connect.qut.edu.au)
##
## This script concatenates all the outputs from Step 4 - 004_calc_MGboundary.R
## script and then summarise the information for each species.
##
## For each species, there is only one intra-species min/max boundary and one
## min/max boundary for the closests inter-species.
##
## OUTPUT
##    * two files RDS and TSV format
##    * a table with one row per species
##    * this output is then used as input for the metabarcoding gap scatter plot
##
## *****************************************************************************



## ****************************************************************************
## setup ----
##

library(argparser)
library(data.table)
library(futile.logger)

TAXA_RANKS <- c('superkingdom', 'kingdom','phylum', 'class', 'order', 'family','genus', 'species')



## *****************************************************************************
## parameters ----

parser <- arg_parser("Summarise the metabarcoding gap boundaries for each species (one row per species)", 
                     name = '005_sumamrise_MGboundary.R',
                     hide.opts = T)

parser <- add_argument(parser, 'input-dir', 
                       help = "Directory path that contains aggregated alignments
                       grouped by species.")

parser <- add_argument(parser, 'modified-metadata',
                       help = 'modified metadata file from Step 3 containing the
                       derepID mapped to taxonomy hierarchy. Only the 7 common ranks 
                       are used {species, genus, family, order, class, phylum, 
                       kingdom} and superkingdom if present. All other columns
                       the metadata file are ignored.')

parser <- add_argument(parser, '--output-prefix', 
                       default = 'MGboundary',
                       help = "Output file (no extension required) after 
                       summarising metabarcoding boundaries.")

parser <- add_argument(parser, '--input-ext', 
                       default = '.rds',
                       help = 'file extensions of input files in {input-dir}')

parser <- add_argument(parser, '--no-log', 
                       flag = T,
                       nargs = 0,
                       help = 'will not generate a separate log file')

parser <- add_argument(parser, '--logDir', 
                       default = 'logs', 
                       help = "specify the log file directory")

parser <- add_argument(parser, '--debugOn', 
                       flag = T, 
                       nargs = 0, 
                       help = "enables debugging messages")

## parse arguments ----
args <- parse_args(parser)


## TEST WITHIN R: uncomment to test within R using demo data provided
# args <- parse_args(parser,
#                    c("04-MGboundary",
#                      "metadata-variant-type.rds",
#                      "--output-prefix", "teleo-MGboundary",
#                      "--debugOn",
#                      "--no-log"))



## ****************************************************************************
## setup logfile ----

source("utils/set_logger.R", local = T)



## *****************************************************************************
## parser parameters ----
##

PREFIX    <- args$output_prefix
EXT       <- ifelse(startsWith(args$input_ext, '*'),
                    args$input_ext, 
                    paste0('*', args$input_ext))
isRDS <- grepl('rds', EXT)
if (!isRDS) {
  flog.error("Non *.rds file not yet supported")
  stop()
}

if (!dir.exists(args$input_dir)) {
  flog.error("Input directory not found: ", args$input_dir)
  stop()
}

if (!file.exists(args$modified_metadata)) {
  flog.error("Taxonomy file not found: %s", args$modified_metadata )
  stop()
}

OUTPUT <- list(rds = sprintf("%s.rds", PREFIX),
               tsv = sprintf("%s.tsv", PREFIX))

if (file.exists(OUTPUT$rds)) {
  flog.info("Output file [%s] already exists, process stopped.",
            OUTPUT$rds)
  quit(status = 0)
}

if (file.exists(OUTPUT$tsv)) {
  flog.info("Output file [%s] already exists, process stopped.",
            OUTPUT$tsv)
  quit(status = 0)
}  


flog.info("
###############################################################################
#
# %s
# %s
#
# PARAMETERS:
#              Input directory: %s
#       Taxonomy meatdata file: %s
#               File extension: %s
#
#                      Logfile: %s
#                        Debug: %s
#
# OUTPUT FILES:
#            Rdata object file: %s
#           Tab-separated file: %s
#
###############################################################################
", 
          parser$name,
          parser$description,
          args$input_dir,
          args$modified_metadata,
          EXT,
          logFile,
          args$debug,
          OUTPUT$rds,
          OUTPUT$tsv)



## ****************************************************************************
## concat all tables ----

MGfiles <- list.files(args$input_dir, 
                      pattern = EXT,
                      recursive = T,
                      full.names = T)
flog.info("Num files: %s", length(MGfiles))

dat <- rbindlist(lapply(MGfiles, function(fn) {
  if (file.exists(fn)) {
    flog.info("Parsing: %s", basename(fn))
    tmp <- readRDS(fn)
    tmp
  } else {
    flog.error("File not exists: %s", fn)
    stop()
  }
}), use.names=T, fill=T)

flog.info("Loaded: %d total records belonging to %d query species", 
          dat[,.N],
          dat[,uniqueN(query.species)])



## ****************************************************************************
## summarise per species ----
##
##    for each species, keep the
##    * minimum and maximum intra.species similarity value
##    * minimum and maximum inter.species similarity of the next closest species


if (hasName(dat, 'closest.interpecies')) {
  dat$closest.interspecies <- NULL
}
dat[, closest.interspecies:=max(inter.maxP), query.species]

absLimit <- dat[closest.interspecies == inter.maxP]
absLimit <- absLimit[, .SD[which.min(inter.rank)], by=query.species]

## calculate whether query species is above or below the diagonal line
absLimit[,MG.type := ifelse(intra.minP-inter.maxP == 0, 
                            'on',
                            ifelse(intra.minP - inter.maxP < 0, 
                                   'overlap',
                                   'gap'))]
absLimit[,MG.type:=factor(MG.type, levels =c('overlap','on','gap'))]


flog.info("Overview of metabarcoding gaps",
          absLimit[, .(nSpecies = uniqueN(query.species)), MG.type],
          capture = T)


## ****************************************************************************
## merge with taxonomy lineage ----
##

## load metadata
flog.info("Loading: %s", args$modified_metadata)

if (endsWith(args$modified_metadata, ".rds")) {
  metadata <- readRDS(args$modified_metadata)
} else if (endsWith(args$modified_metadata, ".tsv")) {
  metadata <- fread(args$modified_metadata)
}

if (!hasName(metadata, "variant.type")) {
  flog.error("Metadata file does not have the {variant.type} column from Step 3 - please check")
  stop()
}

## update column names to match format in absLimit {query.species}
cn <- paste0('query.', TAXA_RANKS)
setnames(metadata, TAXA_RANKS, cn)
setnames(metadata,
         c('group', 'group.rank'), 
         c('query.group', 'query.group.rank'))
cn <- c(cn, 'query.group', 'query.group.rank', 'nVariants', 'variant.type')  

tmp <- unique(metadata[,cn, with=F])[absLimit, on=.(query.species)]

if (tmp[,.N] != absLimit[,.N]) {
  ## should try to avoid this situation, for now throw warning message
  flog.warning("Number of query.species is not matching between alignment data and
               taxonomy metadata. There might be duplicate species names with
               different lineages. Check the taxonomy file.")
}
absLimit <- tmp; rm(tmp)


flog.info("Saving MG boundary for [%d query species] as two output files:
1) %s
2) %s", 
          absLimit[,uniqueN(query.species)],
          OUTPUT$rds,
          OUTPUT$tsv)

saveRDS(absLimit, OUTPUT$rds)
fwrite(absLimit, 
       file = OUTPUT$tsv,
       sep = "\t",
       row.names = F)


## *****************************************************************************
## finished ----
##

flog.info("FINISHED", sessionInfo(), capture = T)
invisible(flog.appender(appender.console()))
