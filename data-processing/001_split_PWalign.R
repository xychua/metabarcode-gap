## ****************************************************************************
## Split and sort pairwise alignment file
##
## Author: Xin-Yi Chua (x.chua@connect.qut.edu.au)
##
##
## This script splits a single pairwise alignment text file from VSEARCH/USEARCH
## output format into files based on the dereplicated amplicon ID (derepID/queryID).
##
## As a single pairwise alignment file can be GB in size, the split process is
## carried out in 2-steps to allow support calling from within PBS scripts.
##
##  Step 1) Splits a single input into parts with L number of lines using 
##          Linux based command 'split'.
##          This step can be skipped if the input file is already split by 
##          lines outside of this script. 
##          The script ASSUMES the first input is an uncompressed text file.
##          If input is gzipped, then split the input file first outside 
##          of this script, then run step 2 of this script.
##
##  Step 2) Expects a directory of spliced files and parses each file, sorting
##          the alignments based on the input derepID (queryID).
##          If step 1 is skipped, then this step assumes the {input_pwAlign} 
##          parameter points to a directory of multiple spliced files.
## 
##
## NOTES:
## Requirement for full pairwise alignment: usually with pairwise alignments 
## being symmetrical such that the value for $(i,j) == (j,i)$ we would only 
## need to keep either the upper or lower triangle of a pairwise matrix, that 
## is, only retain (i,j). This is the usual default behaviour from pairwise 
## alignments like VSEARCH/USEARCH which only provides matches for one
## direction. However, one usage of this data is to allow a user to query the
## data to focus on a (or a set of) species of interests, and extract all
## sequences belonging to that set. This means only extracting rows from the
## pairwise matrix of sequences belonging to the species of interest. If we only
## retained half the matrix, then there would be missing values when sub-setting
## by rows. As such, during step 1 we perform the reciprocal alignments when
## sorting the input.
##
##
## ****************************************************************************





## ****************************************************************************
## setup ----
##

library(argparser)
library(data.table)
library(futile.logger)
library(doParallel)
library(parallel)



## ****************************************************************************
## parameters ----
##

parser <- arg_parser("Split and sort pairwise global alignment by queryID",
                     name = '001_split_PWalign.R', 
                     hide.opts = T)

parser <- add_argument(parser, 'input_pwAlign', 
                       help = 'Pairwise global alignment text file, output from 
                       VSEARCH/USEARCH (e.g. nt.201905__teleo__PWaln.gz)')

parser <- add_argument(parser, '--outdir', 
                       default = '01-splitPWaln', 
                       help = 'Output directory')

parser <- add_argument(parser, '--prefix', 
                       default = 'part', 
                       help = "file prefix after splitting files by lines")

parser <- add_argument(parser, '--lines', 
                       short = '-L', 
                       default = as.integer(1000000),
                       help = "put NUMBER lines per output file using Bash `split` command")

parser <- add_argument(parser, '--skip-split-lines', 
                       short = '-P', 
                       flag = T, 
                       nargs = 0,  
                       help = "turn on to skip splitting the input into parts of
                       [L] lines if already split. The script will fail if there are no parts.")

parser <- add_argument(parser, '--skip-split-query', 
                       short = '-Q', 
                       flag = T, 
                       nargs = 0, 
                       help = "turn on to skip splitting input by queryID")

parser <- add_argument(parser, '--file-start',
                       default = 1,
                       help = "only appilicable for Step 2 when splitting by queryID.
                       This is the file index position from which to start process.")

parser <- add_argument(parser, '--file-num-length',
                       default = -1,
                       help = "only applicable for Step 2 when splitting by queryID.
                       This is the number of files to process from --file-start.
                       If -1, means to process all files.")

parser <- add_argument(parser, '--num-clusters',
                       default = 2,
                       help = 'set the number of clusters for parallel processing')

parser <- add_argument(parser, '--no-log', 
                       flag = T, 
                       nargs = 0,
                       help='will not generate a separate log file')

parser <- add_argument(parser, '--logDir', 
                       default = 'logs', 
                       help="specify the existing log directory")

parser <- add_argument(parser, '--debugOn', 
                       flag = T, 
                       nargs = 0, 
                       help = "to turn ON debug mode")

## get parameters from command line ----
args <- parse_args(parser)

## TEST WITHIN R: uncomment to test within R using demo data that is already spliced into parts
# args <- parse_args(parser, c('example-data/00_parts',
#                              '--prefix', 'nt.201905__teleo__part',
#                              '--skip-split-lines',
#                              '--debugOn'))



## ****************************************************************************
## detect operating system ----

if (Sys.info()['sysname'] == 'Windows') {
  warning("
!! NOTE: this script relies on LINUX bash 'split' commandline to perform Step 1. 
You might need to split the input file manually before using this script or 
you can use the {--skip-split-part} parameter.
")
}


## ****************************************************************************
## setup logfile ----

source("utils/set_logger.R", local = T)


## ****************************************************************************
## parser parameters ----
##

PREFIX    <- args$prefix
IN        <- args$input_pwAlign
IN.DIR    <- dirname(IN)
OUT.DIR   <- args$outdir

## check input filetype is not compressed
if (summary(file(IN))$class == 'gzfile') {
  stop("Detected the input file is GZIPPED, this script expects input text file or split the file externally into parts before using this script")
}


if (!file.exists(IN)) {
  flog.error("Input file not found: %s", IN)
  stop()
}

if (!dir.exists(OUT.DIR)) {
  flog.info("Creating output directory: %s", OUT.DIR)
  dir.create(OUT.DIR)
}

flog.info("
###############################################################################
#
# %s
# %s
#
# PARAMETERS:
#       Input global alignment: %s
#            Split file Prefix: %s
#
#    Skip splitting into parts: %s
#      Num lines per tmp parts: %d
#  Skip splitting into queries: %s
#             From file[start]: %s 
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
          IN,
          PREFIX,
          args$skip_split_lines,
          args$lines,
          args$skip_split_query,
          args$file_start,
          args$file_num_length,
          logFile,
          args$debugOn,
          OUT.DIR)


## ****************************************************************************
## parallel processing ----
##
##    * set the number of clusters

flog.debug("Registering num clusters: %d", args$num_clusters)
cl <- makeCluster(args$num_clusters)
registerDoParallel(cl)



## ****************************************************************************
## SPLIT-BY-LINES ----
##
##    * split the input file into parts of [L] lines
##    * this is required before splitting the files by queryID


if (!args$skip_split_lines) {
  ## check input is a file and not a directory
  if (dir.exists(IN)) {
    flog.error("Input path is a directory, use --skip-split-lines: %s", IN)
    stop()
  }
  
  flog.info("Split into parts with %d lines", args$lines)
  partDir <- sprintf("%s/tmp", OUT.DIR)
  flog.debug("Temp directory: %s", partDir)
  if (!dir.exists(partDir)) {
    dir.create(partDir)
  }
  
  outPREFIX <- sprintf('%s/%s', partDir, PREFIX)
  flog.info(cmd <- sprintf('split -d -a5 -l %d %s %s', args$lines, IN, outPREFIX))
  exit.status <- system(cmd)
  if (exit.status != 0) {
    flog.error("Split command exit: %d", exit.status)
  } else {
    split_files <- list.files(partDir)
    flog.info("Split completed: %d files generated", length(split_files))
  }
}


## ****************************************************************************
## SPLIT-BY-QUERY ----
##
##    * expected input is a directory of split text files
##    * this section will then parse each file sort the alignments based on
##      the queryID (derepID) column
##

if (!args$skip_split_query) {
  if (!dir.exists(IN)) {
    stop('Missing input directory: ', IN)
  }
  partDir <- IN
  
  flog.info("Parsing parts and split by queryID")
  iterFiles <- list.files(partDir, pattern=PREFIX, full.names = T)
  if (length(iterFiles) <= 0) {
    stop("No files found in: ", partDir)
  }
  flog.info("Total number of files found: %d", length(iterFiles))
  
  ## support resuming from different file positions
  START     <- ifelse(is.na(args$file_start) || args$file_start == -1, 1, args$file_start)
  ## only work with non existing output
  if (args$file_num_length == -1) {
    END <- length(iterFiles)
  } else {
    END <- min(START + args$file_num_length-1, length(iterFiles))
  }
  if (START <= END) {
    iterFiles <- iterFiles[START:END]
  } else {
    flog.error("START > END iterFiles: %d > %d", START, END)
    stop()
  }
  
  if (length(iterFiles) == 0) {
    flog.error("There are no part files found in: %s ", partDir)
    stop("There are no part files found in: ", partDir)
  }

  ## switch to mclapply on Linux env
  flog.info("Num files to iterate: %d", length(iterFiles))
  tmp <- lapply(iterFiles, function(fn) {
    flog.debug("Parsing: %s", basename(fn))
    part <- gsub('.*(part_[0-9]{5})','\\1',fn)
    pw <- fread(fn, sep='\t', header=F, showProgress=F)
    
    ## remove if first line contains header
    if (pw[1,1] == 'query') {
      pw <- pw[-1,]
    }
    
    if (ncol(pw)==12) {
      names(pw) <- c('query','target','pident','aln_len','mismatch','gaps',
                     'qstart','qend','tstart','tend','evalue','bit_score')
    } else {
      names(pw) <- c('query','target','pident','aln_len','mismatch','gaps',
                     'qstart','qend','tstart','tend','evalue','bit_score','identities')
    }
    
    ## create the reverse pairwise alignment (see details in above README why)
    ## switch the column names so that the 'query' is always the group-by-factor
    pw.rev <- data.table(pw)
    setnames(pw.rev, 
             c('qend','qstart','query','target','tend','tstart'),
             c('tend','tstart','target','query','qend','qstart'))
    pw.bidirect <- rbind(pw, pw.rev, use.names=T, fill=T)
    
    groups <- split(pw.bidirect, pw.bidirect$query)
    flog.info("Splitting files into [%d] derep amplicons", length(groups))
    tmp <- lapply(groups, function(g){
      out <- strsplit(unique(g$query),";",fixed=T)[[1]][1]
      out <- gsub('\\|','__', gsub("/","-",out))
      outFN <- file.path(OUT.DIR, out)
      flog.trace("\t%10s %5d %s", part, nrow(g), gsub("__","\t",out))
      fwrite(g, file=outFN, sep='\t', append=T, quote=F, row.names=F)
    })
    rm(pw, pw.rev, pw.bidirect, groups)
    gc()
  })
}


## ****************************************************************************
## finished ----
##

flog.info("FINISHED", sessionInfo(), capture = T)
invisible(flog.appender(appender.console()))
