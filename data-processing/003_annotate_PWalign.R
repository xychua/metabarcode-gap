## ****************************************************************************
## Annotate PW alignments with taxonomy information
##
## Author: Xin-Yi Chua (x.chua@connect-qut.edu.au)
##
##
## This script parses all files (recursively) given an input directory of
## spliced Rdata objects following the taxonomy hierarchy structure from 
## Step 2 - 002_merge_taxa_groups.R script. This step annotates the query 
## and target amplicon sequences with their taxonomy information.
##
## This requires a metadata file that maps the {queryID/derepID} to a taxonomy
## lineage. The script expects the 8 common taxonomy ranks:
##    superkingdom, kingdom, phylum, class, order, family, genus, species
##
## The script also:
##    * checks that the number of pairwise alignments are correct
##    * calculates the alignment coverage
##    * assigns the lowest taxonomy rank of query/target pair (see OUTPUT)
##
## OUTPUT
##    * New directory 03-taxaGroups by default is created with directory structure
##      mirroring that of the input directory.
##    * Rdata objects (*.rds) output in corresponding sub-directory
##    * Data objects are data.tables with 9 columns:
##
##      COLUMN              DESCRIPTION
##      query               query ID (derepID) of the representative amplicon
##                          sequence after dereplication
##      query.species       name of query species
##      target              target ID (derepID) of the representative amplicon
##                          sequence after dereplication
##      target.species      name of target species
##      pident              percent identity from pairwise alignment range {0..100}
##      coverage            alignment coverage range {0..1}
##      same.rank           lowest taxonomy rank where query/target is the same
##                          e.g., they are the same familiy
##      same.taxa           taxonomy name of the same rank
##      query.variant.type  whether the query species has only one amplicon 
##                          variant (single) or multiple variants (multi)
##
## *****************************************************************************




## *****************************************************************************
## setup ----

library(argparser)
library(data.table)
library(futile.logger)
library(doParallel)
library(parallel)

## only use these columns if they exists in the group metadata file
TAXA_RANKS <- c('superkingdom', 'kingdom','phylum', 'class', 'order', 
                'family','genus', 'species')


## load helper functions
source("utils/qc_alignments.R")
source("utils/annot_taxonomy.R")
source("utils/calc_coverage.R")



## *****************************************************************************
## parameters ----
##

parser <- arg_parser('Parse global PW alignments and returns summary statistics.', 
                     name = '003_annotate_PWalign.R', 
                     hide.opts = T)

parser <- add_argument(parser, 'input-dir', 
                       help = 'Global pairwise alignment text file, the output from
                       VSearch/USearch.')

parser <- add_argument(parser, 'metadata',
                       help = 'metadata file containing the derepID mapped to taxonomy
                       hierarchy. Only the 7 common ranks are used {species, genus,
                       family, order, class, phylum, kingdom} and superkingdom if
                       present. All other columns the metadata file are ignored.')

parser <- add_argument(parser, '--outdir', 
                       default='03-taxaAnnot',
                       help = 'Output directory where the directory hierarchy 
                       structure in the --input-dir will be mirrored.')

parser <- add_argument(parser, '--input-ext', 
                       default='.rds',
                       help = 'file extension of files in the --input-dir')

parser <- add_argument(parser, '--annot-type', 
                       default = 'all',
                       help = "annotate which type of query variant, 
                       choose from {'multi', 'single', 'all'}")

parser <- add_argument(parser, '--file-start',
                       default = 1,
                       help = "")

parser <- add_argument(parser, '--file-num-length',
                       default = -1,
                       help = "")

parser <- add_argument(parser, '--prefix',
                       default='nt.201905__teleo__', 
                       help = "file prefix")

parser <- add_argument(parser, '--check-species-name',
                       flag = T, 
                       nargs = 0,
                       help = sprintf('turn on to check species names are unique 
                       and have no lineage conflicts. Assumes existence of
                       columns {%s} in the [metadata] file.', paste(TAXA_RANKS, collapse=',')))

parser <- add_argument(parser, '--num-clusters',
                       default = 2,
                       help = 'set the number of clusters for parallel processing')

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

## parse arguments from command line ----
args <- parse_args(parser)


## TEST WITHIN R: uncomment to test within R using demo data provided
# args <- parse_args(parser,
#                    c('data/02-taxaGroup',
#                      'data/nt.201905__teleo__taxaMetadata.tsv',
#                      '--annot-type', 'all'))



## *****************************************************************************
## setup logfile ----
##

source("utils/set_logger.R", local = T)



## *****************************************************************************
## parse parameters ----
##

if (!args$annot_type %in% c('multi', 'single', 'all')) {
  stop("--annot-type must be either: 'multi', 'single' or 'all'")
}

PREFIX    <- args$prefix
EXT       <- ifelse(startsWith(args$extension, '*'),
                    args$extension, 
                    paste0('*', args$extension))
START     <- ifelse(is.na(args$file_start) || args$file_start == -1, 1, args$file_start)

isRDS <- grepl('rds', EXT)
if (!isRDS) {
  flog.error("Non *.rds file not yet supported")
  stop()
}

if (!dir.exists(args$input_dir)) {
  flog.error("Input directory not found: %s", args$input_dir)
  stop()
}

if (!file.exists(args$metadata)) {
  flog.info("Metadata file not found: %s", args$metadata)
  stop()
}

if (!dir.exists(args$outdir)) {
  flog.debug("Create output directory: %s", args$outdir)
  dir.create(args$outdir)
}


flog.info("
###############################################################################
#
# %s
# %s
#
# PARAMETERS
#       Input global alignment: %s
#       Taxonomy metadata file: %s
#              Annotation type: %s
#
#               File extension: %s
#             From file[start]: %s 
#                    to length: %s (if -1, means all files)
#                       Prefix: %s
#                      Logfile: %s
#                        Debug: %s
#
#             OUTPUT directory: %s
#
###############################################################################
", 
          parser$name,
          parser$description,
          args$input_dir,
          args$metadata,
          args$annot_type,
          EXT,
          START, args$file_num_length,
          PREFIX,
          logFile,
          args$debug,
          args$outdir)


## ****************************************************************************
## parallel processing ----
##
##    * set the number of clusters

flog.debug("Registering num clusters: %d", args$num_clusters)
cl <- makeCluster(args$num_clusters)
registerDoParallel(cl)




## *****************************************************************************
## load metadata ----
##

flog.info("Loading: %s", args$metadata)
metadata <- fread(args$metadata)
setkey(metadata, derepID)
flog.debug("-- loaded: %d queries (derepID)", metadata[,uniqueN(derepID)])

missing.names <- which(!hasName(metadata, TAXA_RANKS))
if (length(missing.names) > 0) {
  flog.error("There are missing columns in the metadata file: ", missing.names, capture = T)
  stop()
}


if (args$check_species_name) {
  ## check taxonomy, that we can use species name for grouping, as sometimes
  ## there are species with the same name but different lineage
  check.taxa <- metadata[,lapply(.SD, uniqueN), .SDcols=head(TAXA_RANKS, -1), species]
  if (all(check.taxa[,-1] == 1)) {
    flog.info("Check species name: PASSED")
  } else {
    flog.error("Check species name: FAILED")
    stop("Check species name: FAILED")
  }
  rm(check.taxa)
}


## calculate the number of amplicon variants per species
metadata[,nVariants:=uniqueN(derepID), by=TAXA_RANKS]
metadata[,variant.type:=ifelse(nVariants==1, 'single', 'multi')]
flog.debug("Num species with SINGLE amplicon variant: %d (%.2f%%)", 
           metadata[nVariants==1,uniqueN(species)],
           metadata[nVariants==1,uniqueN(species)]/metadata[,uniqueN(species)]*100)
flog.debug("Num species with MULTI amplicon variants: %d (%.2f%%)", 
           metadata[nVariants >1,uniqueN(species)],
           metadata[nVariants >1,uniqueN(species)]/metadata[,uniqueN(species)]*100)



## ****************************************************************************
## prepare input files ----

flog.info("Load pairwise alignments...")
pwFiles <- list.files(args$input_dir, pattern = EXT, recursive = T, full.names = T)
pwFiles <- grep('*_(part|rev)[0-9]+.rds', pwFiles, invert=T, value=T)
flog.debug("Num total files: %d", length(pwFiles))

## find all expected output files and operate on those that do not yet exists
iterFiles <- gsub('.rds','',gsub(paste0(args$input_dir, '/'), '', pwFiles))

if (args$annot_type == "all") {
  iterFiles <- metadata[group %in% iterFiles,
                        .N, 
                        .(output = sprintf('%s/%s_%s.rds', args$outdir, group, variant.type),
                          input = sprintf('%s/%s.rds', args$input_dir, group))]
} else {
  iterFiles <- metadata[group %in% iterFiles & variant.type == ANNOT.TYPE,
                        .N, 
                        .(output = sprintf('%s/%s_%s.rds', args$outdir, group, variant.type),
                          input = sprintf('%s/%s.rds', args$input_dir, group))]
}

iterFiles[, exists := file.exists(output)]
flog.debug("Num expected outputs: %d", iterFiles[,.N])
flog.debug("...num already processed: %d", iterFiles[exists == T, .N])

iterFiles <- iterFiles[exists == F]
if (nrow(iterFiles) == 0) {
  flog.info("No input files to process, process stopped.")
  quit(status = 0)
}


## only work with non existing output
if (args$file_num_length == -1) {
  END <- nrow(iterFiles)
} else {
  END <- min(START + args$file_num_length-1, nrow(iterFiles))
}
if (START <= END) {
  iterFiles <- iterFiles[START:END]
} else {
  flog.error("START > END iterFiles: %d > %d", START, END)
  stop()
}



## ****************************************************************************
## annotate pairwise alignments ----
##
##    * load in all pairwise alignments objects recursively given an input directory
##    * annotate the source amplicon and target amplicon

flog.debug("Num files to iterate: %d", nrow(iterFiles))
z <- apply(iterFiles, MARGIN=1, function(fn) {
  rdsFile <- fn['output']
  if (!file.exists(rdsFile)) {
    flog.debug("Processing: %s", fn['input'])
    pairwise <- readRDS(fn['input'])
    setkey(pairwise, query)
    
    ## check output sub-directory
    subdir <- dirname(rdsFile)
    if (!dir.exists(subdir)) {
      flog.info("Creating: %s", subdir)
      dir.create(subdir, recursive = T)
    }
    
    if (!check_num_alignments_ok(pairwise)) {
      missing.targets <- setdiff(metadata$derepID, pairwise$target)
      missing.targets <- setdiff(missing.targets, pairwise$query)
      if (length(missing.targets) > 0) {
        flog.error("There are missing targets")
        stop("Not yet implemented")
      }
    }
    
    if (!check_strand_alignments_ok(pairwise)) {
      stop("Strand alignment not OK - error handling not yet implemented")
    }
    
    ## merge with metadata get variant.type for query
    tmp <- metadata[,.(query=derepID, variant.type)][pairwise, on=.(query)]
    stopifnot(tmp[,.N] == pairwise[,.N])
    pairwise <- tmp; rm(tmp)
    
    
    ## split into single and multiple variants, save single variant output
    ## do not process the single variant output here as the best intra-species
    ## match would be to itself. The best inter-species match would depend
    ## on what taxonomic rank we are assessing and this can be processed
    ## separately. This is a sanity check that should never be TRUE as iterFiles
    ## is already check above.
    if (args$annot_type != 'all') {
      variantPW <- pairwise[variant.type == args$annot_type]
    } else { 
      variantPW <- pairwise
    }
    if (variantPW[,.N] == 0) {
      flog.warn("No [%s] amplicons in data: %s", args$annot_type, basename(fn))
      return(NA)
    }
    
    
    ## split by query derepID and annotate with taxonomy info to get matching rank
    ## remove unneeded variables to free up memory for parallel processing
    groups <- split(variantPW, variantPW$query)
    rm(pairwise, variantPW); gc()
    
    flog.debug("Num groups: %d", length(groups))
    ## switch to mclapply on Linux env
    pwAlignTaxa <- pblapply(groups,
                            annot_taxonomy,
                            taxonomy = metadata,
                            debug = args$debug)
    
    ## convert character columns to factors to save on file space
    dat <- rbindlist(pwAlignTaxa, use.names=T, fill=T)
    if (dat[,.N] != sum(sapply(pwAlignTaxa, nrow))) {
      flog.error("Error combining list to data.table, number of rows don't match")
      return(NA)
    }
    char.cols <- names(which(lapply(dat, class) == 'character'))
    dat[,(char.cols):=lapply(.SD, as.factor), .SDcols=char.cols]
    
    flog.info("Saving: %-50s [%d records]", rdsFile, dat[,.N])
    saveRDS(dat, file=rdsFile)
    
    rm(dat, groups, pwAlignTaxa)
    gc()
  } else {
    flog.info("[%s] - annotated RDS present, skipped", basename(fn))
  }
})


## *****************************************************************************
## finished ----
##

flog.info("FINISHED", sessionInfo(), capture = T)
invisible(flog.appender(appender.console()))
