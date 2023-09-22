#' set_logger
#' 
#' helper script that sets the logging file (if enabled)
#' and logger thresholds
#' 
#' 

logFile <- 'no-log-file'
if (!args$no_log) {
  timestamp <- format(Sys.time(), format='%Y%m%d_%H%M%S')
  if (!dir.exists(args$logDir)) {
    dir.create(args$logDir)
  }
  logFile <- file.path(args$logDir, gsub('.R', paste0('-', timestamp,'.log'), parser$name, fixed=T))
  flog.info("Generate log file: %s", logFile)
  
  ## append log file
  invisible(flog.appender(appender.tee(logFile)))
}

if (args$debugOn) {
  invisible(flog.threshold(DEBUG))
} else {
  invisible(flog.threshold(INFO))
}