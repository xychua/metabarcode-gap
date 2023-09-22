## ****************************************************************************
## README ----
## Author: Xin-Yi Chua (x.chua@connect-qut.edu.au)
##
##
## This script parses all files (recursively) in an input directory 
## and annotates the query and target amplicons with their taxonomy information. 
##
## This requires a metadata file that maps the {derepID/queryID} to a taxonomy
## hierarchy. The script expects the 8 common taxonomy ranks:
##    superkingdom, kingdom, phylum, class, order, family, genus, species
##
## The script also performs:
##    * check the number of pairwise alignments are correct
##    * calculates the alignment coverage
##    * assigns the lowest taxonomy rank of query/target pair (see OUTPUT)
##
## OUTPUT
##
##    * Rdata objects (*.dat)
##    * Tables with 9 columns:
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
##
## ____________________
##   Author: Xin-Yi Chua (x.chua@connect-qut.edu.au)
##  Created: 2022-05-19
## Modified: 2023-04-09
##
# 2023-04-09: - update checking output files
# 2023-04-08: - rename file to help remind run first
#             - update so that --file-num-length is referring to number of inputs
#               with no outputs yet. So if --file-num-length=5, the script will
#               first find which inputs are missing outputs and then operate on 
#               the first 5 of those inputs. If there are 100 inputs of which 
#               the first 20 already have outputs, then it will skip the first 20
#               and operate on inputs #21 to 25.
# 2022-02-09: - fixed if variantPW is empty then return error message
# 2022-06-26: - convert character columns to factors to save space
# 2022-06-23: - change default --file-num-length to allow running from 
#               command line for small data sets
# 2022-06-17: - removed the workaround for --annot-type, keep it simple
# 2022-06-15: - add new parameter '--annot-type' for annotated single pwAlignments
# 2022-06-14: - add in if...else to check if 'multi' and 'single' exists 
#               after split by group. rerun for some array jobs that failed
# 2022-06-13: - change description to differentiate between summarise_PWalign.R script
#             - parallelise annotation processing
#
## *****************************************************************************


## *****************************************************************************
## SETUP ----

library(argparser)
library(data.table)
library(futile.logger)
library(pbapply)
# library(parallel) ## for Linux env

## only use these columns if they exists in the group metadata file
TAXA_RANKS <- c('superkingdom', 'kingdom','phylum', 'class', 'order', 
                'family','genus', 'species')


## load helper functions
source("utils/qc_alignments.R")
source("utils/annot_taxonomy.R")
source("utils/calc_coverage.R")



## *****************************************************************************
## PARAMETERS ----
##

parser <- arg_parser('Parse global PW alignments and returns summary statistics.', 
                     name = '003_annotate_PWalign.R', 
                     hide.opts = T)

parser <- add_argument(parser, 'pwAlignment_path', 
                       help = 'Global pairwise alignment text file, the output from
                       VSearch/USearch. If [--recursive] is set then, will assume
                       this is a directory path containing already spliced pairwise
                       alignments files with [--extension].')

parser <- add_argument(parser, 'metadata',
                       help = 'metadata file containing the derepID mapped to taxonomy
                       hierarchy. Only the 7 common ranks are used {species, genus,
                       family, order, class, phylum, kingdom} and superkingdom if
                       present. All other columns the metadata file are ignored.')

parser <- add_argument(parser, '--outdir', 
                       default='03-taxaAnnot',
                       help = 'Output directory. 
                       If [--recursive] is enabled, then the directory 
                       hierarchy structure in [pwAlignment] will be mirrored 
                       in the output directory')

parser <- add_argument(parser, '--extension', 
                       default='.rds',
                       help = 'if --recursive is enabled then the file extensions
                       of the split pairwise alignment files, [default=*.rds]')

parser <- add_argument(parser, '--file-start',
                       default = 1,
                       help = "for PBS/SLURM or running in parallel mode, to specify 
                       which file to start from. Only used when [--recursive] 
                       is enabled")

parser <- add_argument(parser, '--file-num-length',
                       default = -1,
                       help = 'for PBS/SLURM or running in parallel mode, specify the 
                       number of files to iterate. Only used when [--recursive] 
                       is enabled, otherwise is (default=-1) which will process
                       all files.')

parser <- add_argument(parser, '--annot-type', 
                       default = 'all',
                       help = "annotate which type of query variant, 
                       choose from {'multi', 'single', 'all'}")

parser <- add_argument(parser, '--prefix',
                       default='nt.201905__teleo__', 
                       help = "file prefix")

parser <- add_argument(parser, '--check-species-name',
                       flag = T, 
                       nargs = 0,
                       help = sprintf('turn on to check species names are unique 
                       and have no lineage conflicts. Assumes existence of
                       columns {%s} in the [metadata] file.', paste(TAXA_RANKS, collapse=',')))

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


## TEST WITHIN R: uncomment for testing within the script
# args <- parse_args(parser,
#                    c('02-taxaGroups',
#                      'data/nt.201905__teleo__taxaMetadata.tsv',
#                      '--annot-type', 'all'))



## *****************************************************************************
## LOG FILE ----
##

logFile <- 'no-log-file'
if (!args$no_log) {
  timestamp <- format(Sys.time(), format='%Y%m%d_%H%M%S')
  if (!dir.exists(args$logDir)) {
    dir.create(args$logDir)
  }
  logFile <- file.path(args$logDir, gsub('.R', paste0('_', timestamp,'.log'), parser$name, fixed=T))
  flog.info("Generating log file: %s", logFile)
  
  ## append log file
  invisible(flog.appender(appender.tee(logFile)))
  
  if (args$debugOn) {
    invisible(flog.threshold(DEBUG))
    pboptions(type='none')
  } else {
    invisible(flog.threshold(INFO))
    pboptions(type='txt')
  }
}



## *****************************************************************************
## CHECK PARAMETERS ----
##

if (!args$annot_type %in% c('multi', 'single', 'all')) {
  stop("--annot-type must be either: 'multi', 'single' or 'all'")
}

PREFIX    <- args$prefix
IN        <- list(pwAlignment = args$pwAlignment, 
                  metadata = args$metadata)
EXT       <- ifelse(startsWith(args$extension, '*'),
                    args$extension, 
                    paste0('*', args$extension))
START     <- ifelse(is.na(args$file_start) || args$file_start == -1, 1, args$file_start)

isRDS <- grepl('rds', EXT)
if (!isRDS) {
  flog.error("Non *.rds file not yet supported")
  stop()
}

missing_files <- IN[which(!sapply(IN, file.exists))]
if (length(missing_files) > 0) {
  stop("Missing input file(s):\n", missing_files)
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
#          Group metadata file: %s
#
#                    Recursive: %s
#               File extension: %s
# From file[start] to [length]: %s to %s
#              Annotation type: %s
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
          IN$pwAlignment,
          IN$groupFile,
          args$recursive,
          EXT,
          START, args$file_num_length,
          args$annot_type,
          PREFIX,
          logFile,
          args$debug,
          args$outdir)




## *****************************************************************************
## LOAD METADATA ----
##

flog.info("Loading: %s", IN$metadata)
metadata <- fread(IN$metadata)
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
## ANNOTATE PAIRWISE ALIGNMENTS -----
##
##    * load in all pairwise alignments objects recursively given an input directory
##    * annotate the source amplicon and target amplicon


flog.info("Load pairwise alignments...")
pwFiles <- list.files(IN$pwAlignment, pattern = EXT, recursive = T, full.names = T)
pwFiles <- grep('*_(part|rev)[0-9]+.rds', pwFiles, invert=T, value=T)
flog.debug("Num total files: %d", length(pwFiles))

## find all expected output files and operate on those that do not yet exists
iterFiles <- gsub('.rds','',gsub(paste0(IN$pwAlignment, '/'), '', pwFiles))

if (args$annot_type == "all") {
  iterFiles <- metadata[group %in% iterFiles,
                        .N, 
                        .(output = sprintf('%s/%s_%s.rds', args$outdir, group, variant.type),
                          input = sprintf('%s/%s.rds', IN$pwAlignment, group))]
} else {
  iterFiles <- metadata[group %in% iterFiles & variant.type == ANNOT.TYPE,
                        .N, 
                        .(output = sprintf('%s/%s_%s.rds', args$outdir, group, variant.type),
                          input = sprintf('%s/%s.rds', IN$pwAlignment, group))]
}

iterFiles[, exists := file.exists(output)]
flog.debug("Num expected outputs: %d", iterFiles[,.N])
flog.debug("...num already processed: %d", iterFiles[exists == T, .N])


iterFiles <- iterFiles[exists == F]
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


## iterate through files
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
    
    flog.info("Saving: %s [%d records]", rdsFile, dat[,.N])
    saveRDS(dat, file=rdsFile)
    
    rm(dat, groups, pwAlignTaxa)
    gc()
  } else {
    flog.info("[%s] - annotated RDS present, skipped", basename(fn))
  }
})


## *****************************************************************************
## FINISHED ----
##

flog.info("FINISHED", sessionInfo(), capture = T)
flog.appender(appender.console())
