## ****************************************************************************
## Calculate MG boundary
##
## Author: Xin-Yi Chua (x.chua@connect-qut.edu.au)
##
##
## This script parses all files (recursively) given an input directory following
## Step 3 - 003_annotate_PWalign.R script. The input files should have taxonomy
## information for the query/target amplicon sequences. 
##
## The pairwise sequence alignments are grouped by query.species and this script
## calculates the metabarcoding gap boundaries for each intra- and inter-rank.
##
## OUTPUT
##    * New directory 04-MGboundary by default is created with directory structure
##      mirroring that of the input directory.
##    * Rdata objects (*.rds) output in corresponding sub-directory
##    * Data objects are data.tables with 10 columns:
##
##      COLUMN              DESCRIPTION
##      query.species       name of query species
##      query.variant.type  whether the query species has only one amplicon 
##                          variant (single) or multiple variants (multi)
##      intra.minP          minimum similarity value in intra-species pairwise distribution
##      intra.maxP          maximum similarity value in intra-species pairwise distribution
##      intra.nPairs        number of sequence pairs in intra-species pairwise distribution
##      inter.rank          taxonomy rank to the next closest inter-species
##                          e.g. species from the same genus, family, ..., phylum
##      inter.taxa          taxonomy name of inter-rank
##      inter.minP          minimum similarity value in inter-rank pairwise distribution
##      inter.maxP          maximum similarity value in inter-rank pairwise distribution
##      inter.nPairs        number of sequence pairs in inter-rank pairwise distribution
##      nTargetSpecies      number of target species with same inter-rank
##
##
## *****************************************************************************



## *****************************************************************************
## setup ----

library(argparser)
library(data.table)
library(futile.logger)
library(doParallel)
library(parallel)

source("utils/get_MG_boundary.R")

TAXA_RANKS <- c('superkingdom', 'kingdom','phylum', 'class', 'order', 'family','genus', 'species')



## *****************************************************************************
## parameters ----

parser <- arg_parser('Parse global PW alignments and returns summary statistics.', 
                     name = '004_MGboundary_PWalign.R',
                     hide.opts = T)

parser <- add_argument(parser, 'input-dir', 
                       help = "Directory path that contains pairwise alignments
                       that have been taxonomically annotated.")

parser <- add_argument(parser, '--outdir', 
                       default = '04-MGboundary',
                       help = "Output directory, where the same taxonomy hierarchical
                       directory structure in {input_dir} will be mirrored.")

parser <- add_argument(parser, '--input-ext', 
                       default = '.rds',
                       help = 'file extensions of input files in {input-dir}')

parser <- add_argument(parser, '--file-start',
                       default = 1,
                       help = "to specify which file to start from. Primarily
                       for use in HPC environments with job queues (e.g. PBS/SLURM).
                       Allows for processing subset of input files or to resume
                       when job files.")

parser <- add_argument(parser, '--file-num-length',
                       default = -1,
                       help = 'used with --file-start, to specify the number of
                       files to process from the start index. If -1, means 
                       to process all files.')

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

## parse arguments ----
args <- parse_args(parser)


## TEST WITHIN R: uncomment to test within R using demo data provided
# args <- parse_args(parser,
#                    c("03-taxaAnnot",
#                      '--debugOn',
#                      '--no-log'))




## ****************************************************************************
## setup logfile ----

source("utils/set_logger.R", local = T)



## *****************************************************************************
## parse parameters ----
##

EXT       <- ifelse(startsWith(args$input_ext, '*'),
                    args$input_ext, 
                    paste0('*', args$input_ext))
START     <- ifelse(is.na(args$file_start) || args$file_start == -1, 1, args$file_start)

isRDS <- grepl('rds', EXT)
if (!isRDS) {
  flog.error("Non *.rds file not yet supported")
  stop()
}

if (!dir.exists(args$input_dir)) {
  stop("Input directory doesn't exists: ", args$input_dir)
}

if (!dir.exists(args$outdir)) {
  flog.debug("Creating: %s", args$outdir)
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
#
#               File extension: %s
#                   From start: %s
#                    to length: %s (if -1, means all files)
#
#                      Logfile: %s
#                        Debug: %s
#
# OUTPUT FILES:
#             Output directory: %s
#
###############################################################################
", 
          parser$name,
          parser$description,
          args$input_dir,
          EXT,
          START, 
          args$file_num_length,
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
## prepare files ----
##

flog.info("Load input files ...")
pwFiles <- list.files(args$input_dir, pattern=EXT, recursive=T, full.names = T)
flog.debug("Num total files: %d", length(pwFiles))

## check if file already exists
flog.info("Checking expected output files ...")
expected <- data.table(input=pwFiles, 
                       output=gsub(args$input_dir,
                                   args$outdir,
                                   pwFiles))
expected[, exists:=file.exists(output)]
flog.debug("Num expected files: %d", expected[,.N])
flog.debug("...num already processed files: %d", expected[exists == T, .N])
iterFiles <- expected[exists == F]

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
  flog.error("START file [%d] is greater than END [%d]", START, END)
  stop()
}


## *****************************************************************************
## processing ----
##

flog.debug("Num files to iterate: %d", nrow(iterFiles))
z <- apply(iterFiles, MARGIN = 1, function(fn) {
  rdsFile <- fn[['output']]
  inFile <- fn[['input']]
  
  if (!file.exists(rdsFile)) {
    pwAlignTaxa <- readRDS(inFile)
    
    ## convert factors to characters
    factor.cols <- c('query.species', 'target.species','same.taxa')
    if (!all(hasName(pwAlignTaxa, factor.cols))) {
      flog.error("Expected column names not found: ", factor.cols, capture=T)
      stop()
    }
    pwAlignTaxa[, (factor.cols):=lapply(.SD, as.character), .SDcols = factor.cols]
    
    flog.debug("Processing: %s - %d records", inFile, pwAlignTaxa[,.N])
    
    if ('query.species' %in% names(pwAlignTaxa)) {
      flog.debug("Num query species: %d", pwAlignTaxa[,uniqueN(query.species)])
    } else {
      flog.error("Missing {query.species} column ...")
      print(names(pwAlignTaxa))
    }
    pwSummary <- get_MG_boundary(pwAlignTaxa, debug=args$debug)
    
    ## mirror the taxonomy sub-directories for output
    subdir <- dirname(rdsFile)
    if (!dir.exists(subdir)) {
      flog.debug("Create: %s", subdir)
      dir.create(subdir, showWarnings = F, recursive=T)
    }
    
    flog.info("Saving: %-50s [%d records]", basename(rdsFile), pwSummary[,.N])
    saveRDS(pwSummary, file=rdsFile)
  } else {
    flog.warn("Already process -> skipped: [%s]", basename(inFile))
  }
})

  
  
## *****************************************************************************
## finished ----
##

flog.info("FINISHED", sessionInfo(), capture = T)
invisible(flog.appender(appender.console()))
