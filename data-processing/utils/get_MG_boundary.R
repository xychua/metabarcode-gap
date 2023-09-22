#' get_MG_boundary
#'
#' This function will summarise and extract the boundaries of intra- and
#' inter-rank pairwise alignments.
#'
#' @param pwAlignTaxa data.table of pairwise alginments that have been
#'                    taxonomically annotated and contains the {same.taxa}
#'                    column from \link{annot_taxonomy}
#' @return MG.boundary data.table with columns:
#'         1. query.species
#'         2. same.rank
#'         3. same.taxa
#'         4. intra.minP
#'         5. intra.maxP
#'         6. intra.nPairs
#'         7. nTargetSpecies
#'         8. inter.minP
#'         9. inter.maxP
#'        10. inter.nPairs

get_MG_boundary <- compiler::cmpfun(function(pwAlignTaxa, debug=F) {

  if (debug) { flog.threshold(DEBUG) }

  if (!'coverage' %in% names(pwAlignTaxa)) {
    flog.error("Missing {coverage} column attribute")
  } else {
    
    ## ****
    ## calculate INTRA-species boundaries
    ##
    ##    - split the data based on the query.variant.type
    ##    - species with SINGLE amplicon variant has highest match to themselves
    
    intra.dist <- data.table()
    
    groups <- split(pwAlignTaxa, pwAlignTaxa$query.variant.type)
    if (groups$multi[same.rank == 'species', .N] > 0) {
      intra.dist <- rbind(intra.dist,
                          groups$multi[same.rank == 'species',
                                       .(intra.minP = min(pident),
                                         intra.maxP = max(pident),
                                         intra.nPairs = .N,
                                         nTargetSpecies = uniqueN(target.species)),
                                       .(query.species, same.rank, same.taxa, query.variant.type)],
                          use.names = T, 
                          fill = T)
    }

    if (groups$single[,.N] > 0) {    
      flog.trace("Creating self-matches for query species with SINGLE amplicon variants")
      ## the self-match is the only possible same-species comparison
      intra.dist <- rbind(intra.dist, 
                          unique(groups$single[,.(query.species,
                                                  same.rank = 'species',
                                                  same.taxa = query.species,
                                                  query.variant.type,
                                                  intra.minP = 100,
                                                  intra.maxP = 100,
                                                  intra.nPairs = 1,
                                                  nTargetSpecies = 1)]),
                          use.names = T,
                          fill = T)
    }
    
    flog.debug('intra-species: %d query.species - %d records',
               intra.dist[,uniqueN(query.species)],
               intra.dist[,.N])
    
    
    ## *****
    ## calculate INTER-species boundaries
    inter.dist <- pwAlignTaxa[same.rank!='species',
                              .(inter.minP=min(pident),
                                inter.maxP=max(pident),
                                inter.nPairs=.N,
                                nTargetSpecies=uniqueN(target.species)),
                              .(query.species, same.rank, same.taxa)]
    flog.debug('inter-species: %d query.species - %d records',
               inter.dist[,uniqueN(query.species)],
               inter.dist[,.N])
    
    
    ## merge tables
    if (exists('intra.dist')) {
      stopifnot(identical(intra.dist$query.species, intra.dist$same.taxa))
      stopifnot(unique(intra.dist$same.rank) == 'species')
      stopifnot(unique(intra.dist$nTargetSpecies) == 1)
      ## remove not required columns before merge
      intra.dist[,c('same.taxa', 'same.rank', 'nTargetSpecies'):=NULL]
      setnames(inter.dist, c('same.taxa','same.rank'), c('inter.taxa','inter.rank'))
      MG.boundary <- intra.dist[inter.dist, on=.(query.species)]
    } else {
      flog.error("Missing {intra.dist} data table")
    }
    flog.debug('Return: %d species - %d records',
               MG.boundary[,uniqueN(query.species)],
               MG.boundary[,.N])
    return(MG.boundary)
  }
})
