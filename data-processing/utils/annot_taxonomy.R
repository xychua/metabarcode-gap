#' annot_taxonomy
#'
#' This function annotates the given pairwise alignment with taxonomy lineage
#' metadata. This is repeated for both the \code{query} and \code{target}
#' sequences. The function then adds a new column \code{same.taxa} that denotes
#' at which taxonomic rank the pair of aligned sequences (i, j) are the same.
#'
#' @param pwAlign  data.table of pairwise alignments with \code{query} and
#'                 \code{target} columns
#' @param taxonomy data.table of taxonomy lineage from species to *Superkingdom*
#'                 with \code{derepID} column that is used to match the
#'                 \code{query} and \code{target} columns from \code{pwAlign}
#'                 input
#' @return data.table with assigned taxonomy lineage for the \code{query}
#'         sequence and the \code{same.rank} and \code{same.taxa} columns
#'
#' @import data.table
#' @import futile.logger

annot_taxonomy <- compiler::cmpfun(function(pwAlign,
                                            taxonomy,
                                            debug=F) {
  if (debug) { flog.threshold(DEBUG) }

  flog.debug("Parsing: query %s - %d alignments", pwAlign[,unique(query)], pwAlign[,.N])

  ## get taxonomy for query
  taxa <- taxonomy[derepID %in% pwAlign$query]
  names(taxa) <- paste0('query.',names(taxa))
  setnames(taxa, 'query.derepID', 'query')
  pwAlign.taxa <- taxa[pwAlign, on=.(query)]
  stopifnot(pwAlign.taxa[,.N] == pwAlign[,.N])

  ## get taxonomy for target
  taxa  <- taxonomy[derepID %in% pwAlign$target]
  names(taxa) <- paste0('target.', names(taxa))
  setnames(taxa, 'target.derepID', 'target')
  pwAlign.taxa <- taxa[pwAlign.taxa, on=.(target)]
  stopifnot(pwAlign.taxa[,.N] == pwAlign[,.N])

  ## identify the taxonomic rank for each pair of alignments
  flog.trace("Assigning {same.rank} ...", pwAlign.taxa[,.N])
  pwAlign.taxa[query.species==target.species,
               `:=`(same.rank='species', same.taxa=query.species)]

  pwAlign.taxa[is.na(same.rank) & query.genus==target.genus,
               `:=`(same.rank='genus', same.taxa=query.genus)]

  pwAlign.taxa[is.na(same.rank) & query.family==target.family,
               `:=`(same.rank='family', same.taxa=query.family)]

  pwAlign.taxa[is.na(same.rank) & query.order==target.order,
               `:=`(same.rank='order', same.taxa=query.order)]

  pwAlign.taxa[is.na(same.rank) & query.class==target.class,
               `:=`(same.rank='class', same.taxa=query.class)]

  pwAlign.taxa[is.na(same.rank) & query.phylum==target.phylum,
               `:=`(same.rank='phylum', same.taxa=query.phylum)]

  pwAlign.taxa[is.na(same.rank) & query.kingdom==target.kingdom,
               `:=`(same.rank='kingdom', same.taxa=query.kingdom)]

  pwAlign.taxa[is.na(same.rank) & query.superkingdom==target.superkingdom,
               `:=`(same.rank='superkingdom', same.taxa=query.superkingdom)]

  pwAlign.taxa[,same.rank:=factor(same.rank, levels=TAXA_RANKS, ordered=T)]

  calc_coverage(pwAlign.taxa)

  ## keep minimum columns
  cn <- c('query',
          'query.species',
          'target',
          'target.species',
          'pident',
          'coverage',
          'same.rank',
          'same.taxa',
          'query.variant.type')
  return(pwAlign.taxa[, cn, with=F])
})
