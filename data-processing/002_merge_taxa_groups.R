## *****************************************************************************
## README ----
## Author: Xin-Yi Chua (x.chua@connect-qut.edu.au)
##
##
## This script takes the pairwise global alignments after being split into
## groups based on their derepID/queryID, and reorganises the data based on
## the taxonomy lineage. 
##
## This requires a metadata file that maps the {derepID/queryID} to a
## {group} following the taxonomy hierarchy specified in file path format,
## for example:
##
##  DEREPID               GROUP
##  KF612341.1_q914_a100  Eukaryota/Metazoa/Chordata/Chordata-subgroup-1
##  KT203434.1_q920_a99   Eukaryota/Metazoa/Chordata/Chordata-subgroup-1
##  MH085614.1_q237_a99   Eukaryota/Metazoa/Chordata/Actinopteri/Actinopteri-subgroup-6
##  AP009147.1_q908_a100  Eukaryota/Metazoa/Chordata/Actinopteri/Cypriniformes/Cypriniformes-subgroup-2
##
## * first 2 dereplicated amplicon will be grouped into the Chordata-subgroup-1
## * third one is grouped into Actinopteri-subgroup-6
## * forth one is grouped into Cypriniformes-subgroup-2 etc
##
## The output creates a new directory that follows the hierarchical filepath
## format as per the metadata file. Binary Rdata objects (*.rds) are saved 
## in their corresponding sub-folders.
##
##
## NOTE: After this step, the input files can be deleted.
##
## ____________________
##   Author: Xin-Yi Chua (x.chua@connect-qut.edu.au)
##  Created: 2021-05-18
## Modified: 2022-06-26
##
## 2022-06-26: remove START and END parameter settings
## 2022-06-23: will go ahead and remove the sub-files after combining
## 2022-05-25: new smaller group size (500), output directory nt.20190522.61.taxa
## 2022-05-18: max group size 1000, output directory nt.20190522.62.taxa/
##
## *****************************************************************************




## *****************************************************************************
## SETUP ----
##

library(argparser)
library(data.table)
library(futile.logger)
library(pbapply)
# library(parallel) ## for Linux env


## *****************************************************************************
## PARAMETERS ----
##

parser <- arg_parser('Merge pairwise alignments by taxonomy groups', 
                     name = '002_merge_taxa_groups.R', 
                     hide.opts = T)

parser <- add_argument(parser, 'indir',
                       help = 'input directory with pairwise sequences sorted by queryID')

parser <- add_argument(parser, 'metadata',
                       help = 'metadata file that maps the queryIDs to the taxonomic group')

parser <- add_argument(parser, '--group-col', 
                       default = 'group',
                       help = 'column name in [metadata] file that contains the 
                       groups to which input query files will be assigned')

parser <- add_argument(parser, '--outdir', 
                       default = '02-taxaGroups', 
                       help = 'Output directory')

parser <- add_argument(parser, '--no-log', 
                       flag = T, 
                       nargs = 0,
                       help = 'will not generate a separate log file')
parser <- add_argument(parser, '--logDir',
                       default = 'logs', 
                       help = "specify the existing log directory")
parser <- add_argument(parser, '--debugOn', 
                       flag = T, 
                       nargs = 0,
                       help = "to turn ON debug mode")

## get parameters from command line ----
args <- parse_args(parser)


## TEST WITHIN R: uncomment to test within R using demo data provided
# args <- parse_args(parser, 
#                    c('data/01-split',
#                      'data/nt.201905__teleo__taxaMetadata.tsv',
#                      '--debugOn'))




## ****************************************************************************
## setup logfile ----

source("utils/set_logger.R", local = T)



## ****************************************************************************
## parser parameters ----
##

if (!dir.exists(args$indir)) {
  stop("Missing input directory: ", args$indir)
}

if (!file.exists(args$metadata)) {
  stop("Missing metadata file: ", args$metadata)
}

if (!dir.exists(args$outdir)) {
  flog.info("Creating output directory: %s", args$outdir)
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
#                Metadata file: %s
#                 Group column: %s
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
          args$indir,
          args$metadata,
          args$group_col,
          logFile,
          args$debugOn,
          args$outdir)




## *****************************************************************************
## LOAD METADATA ----
##
##    * user can specify the {group} column using the {--group-col} parameter setting
##    * Check for missing identifiers --> warning message
##

metadata <- fread(args$metadata)
if (!(args$group_col %in% names(metadata))) {
  flog.error("Group column [%s] not found in metadata file", args$group_col)
  stop()
}

pwFiles <- data.table(path=list.files(args$indir, 
                                      full.names = T, 
                                      recursive = F, 
                                      include.dirs = F))
pwFiles[,dir:=dirname(path)]
pwFiles[,derepID:=basename(path)]
pwFiles[,isDir:=file.info(path)$isdir]
pwFiles <- pwFiles[isDir == F]
flog.info("Num files to parse: %d", pwFiles[,.N])


## check if any query files are missing
missing.dereps <- setdiff(metadata$derepID, pwFiles$derepID)
flog.warn("Num missing queries (derepID): %d", length(missing.dereps))
flog.warn("Missing queries (derepIDs): ", missing.dereps, capture = T)



## ****************************************************************************
## MERGE TAXONOMY ----
##
##    * merge taxonomy information with file names for grouping


pwFiles_meta <- metadata[,c('derepID', args$group_col), with=F][pwFiles, on=.(derepID)]
setnames(pwFiles_meta, args$group_col, 'group')

if (pwFiles_meta[is.na(group), .N] > 0) {
  flog.warn("There are queries (derepIDs) with NA assigned group - check metadata")
}



## ****************************************************************************
## MAP QUERIES to GROUPS ----
##
##    * split metadata into groups and iterate through each group to combine input files
##


pwGroups <- split(pwFiles_meta, pwFiles_meta$group, drop = T)
flog.info("Num groups: %d", length(pwGroups))
flog.debug("Num queries per group: ", 
           pwFiles_meta[order(group), .N, group], capture=T)

z <- pblapply(pwGroups, function(grp) {
  grp_name <- grp$group[1]
  rdsFile <- sprintf('%s/%s.rds', args$outdir, grp_name)
  
  if (file.exists(rdsFile)) {
    flog.info("Skipped-already parsed: %s", rdsFile)
  } else if (grp[,.N] > 0) {
    
    ## create sub-directories as required
    subdir <- dirname(rdsFile)
    if (!dir.exists(subdir)) {
      dir.create(subdir, recursive = T)
    }
    
    flog.info("Parsing: %s [%d files]", grp_name, grp[,.N])
    
    dat <- rbindlist(mclapply(grp$path, function(fn) {
      tmp <- fread(fn)
      flog.debug("Parsing: %s [%d records]", basename(fn), tmp[,.N])
      tmp
    }), use.names=T, fill=T)

    flog.info("Saving: %s [%d records]", rdsFile, dat[,.N])
    saveRDS(dat, file=rdsFile)
    
    rm(dat); gc()
  }
})


## ****************************************************************************
## Finished ----
##

flog.info("FINISHED", sessionInfo(), capture = T)
invisible(flog.appender(appender.console()))
