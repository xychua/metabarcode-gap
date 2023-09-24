## ****************************************************************************
## plot MG-scatterplot function ----
##    helper function to generate one MG scatterplot for given data
##    data table must have following columns
##      - intra.minP > x-coordinate
##      - inter.maxP > y-coordinate
##      - inter.rank > shape and colour
##      - interest   > will appeared coloured, otherwise appears as background
##      - show.lbl   > only if 'show.lbl = T', the label to display in graph

plotMGscatter <- function(X, 
                          show.lbl=F, 
                          lbl.size=3,
                          nudge_x = -3,
                          nudge_y = -3) {
  minV <- floor(X[,min(intra.minP,inter.maxP)])
  
  p <- ggplot(X, aes(x=intra.minP, 
                     y=inter.maxP, 
                     col=MG.type, 
                     shape=inter.rank, 
                     size=nVariants)) + 
    geom_abline(intercept=0, slope=1, col='grey30', linewidth=.3) +
    geom_point(col='grey85') +
    geom_point(data=X[interest==T], stroke=1) +
    facet_wrap(primer ~ ., ncol=2) +
    scale_colour_manual(values = MGcolours, name='MG type', drop=F) +
    scale_shape_manual(values = taxaShapes,   name='Closest relative', drop=F) +
    scale_size(range=c(.8, 10), name='nVariants') +
    scale_x_continuous(breaks=seq(0,100,10), 
                       expand=expansion(mult=c(0.05, .3))) +
    scale_y_continuous(breaks=seq(0,100,10), 
                       expand=expansion(mult=c(0.05, .3))) +
    facet_wrap(primer ~ ., ncol=2) +
    coord_equal() +
    labs(x="Least alike SAME species",
         y="Most alike DIFFERENT species") +
    guides(colour=guide_legend(order=1),
           shape=guide_legend(order=1, override.aes=list(size=3)),
           size=guide_legend(order=2, override.aes=list(shape=1))) +
    theme(axis.title.x = element_text(size=14, margin=margin(3,1,1,1, 'mm')),
          axis.title.y = element_text(size=14, margin=margin(1,5,1,1, 'mm')),
          panel.background = element_rect(colour='grey30', linewidth=.3),
          panel.spacing = unit(5, 'mm'),
          strip.text.x = element_text(size=14, face='bold', colour='grey30', hjust=0),
          strip.background = element_blank(),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10),
          legend.position = 'right',
          legend.justification = 'top',
          legend.direction ='vertical',
          legend.background = element_blank())
  
  if (show.lbl) {
    vlabels <- unique(X[!is.na(show.lbl) & MG.type %in% c('overlap', 'on'), 
                        .(show.lbl, interest, intra.minP, inter.maxP, primer, inter.rank, MG.type)])
    hlabels <- unique(X[!is.na(show.lbl) & MG.type == 'gap', 
                        .(show.lbl, interest, intra.minP, inter.maxP, primer, inter.rank, MG.type)])
    
    p <- p + ggnewscale::new_scale_colour() +
      geom_text_repel(data=hlabels, aes(label=show.lbl, col=MG.type), 
                      show.legend = F,
                      size = lbl.size, 
                      #key_glyph='point',
                      box.padding = .1, 
                      max.overlaps = 50,
                      segment.size = .3,
                      direction = 'y',
                      hjust = 0,
                      nudge_y = nudge_y,
                      nudge_x = 105 - hlabels$intra.minP,
                      parse = T) +
      geom_text_repel(data=vlabels, aes(label=show.lbl, col=MG.type), 
                      show.legend = F,
                      size = lbl.size, 
                      #key_glyph='point',
                      box.padding = .1, 
                      max.overlaps = 50,
                      segment.size = .3,
                      angle = 90, 
                      direction = 'x', 
                      hjust = 0,
                      vjust = 0,
                      nudge_x = nudge_x,
                      nudge_y = 103 - vlabels$inter.maxP,
                      parse = T) +
      scale_colour_manual(values = MGcolours, name='MG region')
  }
  return(p)
}


## ****************************************************************************
## add labels ----

add_repel_labels <- function(gg_plot, lbl.size) {
  return(gg_plot + ggnewscale::new_scale_colour() +
           geom_text_repel(aes(label=show.lbl, col=MG.type),
                           size = lbl.size, 
                           box.padding = .1, 
                           max.overlaps = 100,
                           segment.size = .3,
                           parse=T) +
           scale_colour_manual(values = MGcolours, name='MG region'))
}
