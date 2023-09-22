#' calc_coverage
#' 
#' This function calculates the coverage percentage of each alignment in a
#' pairwise alignment data.table
#' 
#' @param pwAlign.taxa data.table of pairwise alignments with original amplicon
#'                     lengths for {query} and {target}
#' @return new column attribute {coverage} added
#' 
calc_coverage <- compiler::cmpfun(function(pwAlign.taxa) {
  if ('coverage' %in% pwAlign.taxa) {
    stop("Column {coverage} already exists in pwAlign.taxa data.table")
  }
  pwAlign.taxa[,coverage := 100*(aln_len-gaps)/max(query.amp_length, target.amp_length), .(query, target)]
  
  flog.trace("Coverage range: %.3f - %.3f", 
             pwAlign.taxa[,min(coverage)], 
             pwAlign.taxa[,max(coverage)])
})