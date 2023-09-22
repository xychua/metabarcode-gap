## Alignment QC checks

#' check_num_alignments_ok
#'
#' This function checks that the number of pairwise alignments are
#' as expected. Where the number of expected alignments are
#' numQueries * numTargets - nuQueries (self-matches).
#'
#' @param pairwise data.table of pairwise alignments with query and target
#'                 identifiers. Expect column attributes {query} and {target}
#' @return True/False the number of alignments is the number expected
#' @examples
#' \dontrun{
#' check_num_alignments_ok(pairwise)
#' }
#'
#' @import data.table
#' @import futile.logger

check_num_alignments_ok <- compiler::cmpfun(function(pairwise) {
  cn <- c('query','target')
  if (!all(cn %in% names(pairwise))) {
    stop("Missing column attributes: ", cn)
  }

  nAlign <- pairwise[,.N]
  nQueries <- pairwise[,uniqueN(query)]
  nTargets <- pairwise[,uniqueN(target)]

  flog.debug("Loaded: %d alignments for %d queries %d targets", nAlign, nQueries, nTargets)

  nExpected = (nQueries * nTargets) - nQueries
  if (nAlign > nExpected) {
    flog.warn("Number of alignments is greater than expected")
  }

  if (nAlign < nExpected) {
    flog.warn("Number of alignments is less than expected")
  }

  return(nAlign == nExpected)
})


#' check_strand_alignments_ok
#'
#' This function checks that the {start} position is always smaller than the
#' {end} position in the pairwise alignments. It checks both the
#' {query} and {target}.
#'
#' @param pairwise data.table of pairwise alignments with query and target
#'                 identifiers. Expect column attributes {qstart, qend, tstart}
#'                 and {tend}.
#' @return True/False
#'
check_strand_alignments_ok <- compiler::cmpfun(function(pairwise) {
  if (!all(cn <- c('qstart','qend','tstart','tend') %in% names(pairwise))) {
    stop("Missing column attributes: ", cn)
  }


  if (pairwise[qstart > qend, .N] > 0) {
    flog.warn("Some queries have [start > end], match on the reverse strand")
    return(FALSE)
  }
  if (pairwise[tstart > tend, .N] > 0) {
    flog.warn("Some targets have [start > end], match on the reverse strand")
    return(FALSE)
  }
  return(TRUE)
})
