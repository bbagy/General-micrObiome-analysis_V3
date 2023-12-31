
calRare <- function(psdata, measures, depths, parallel=FALSE) {
  require('plyr') # ldply
  require('reshape2') # melt
  require('doParallel')
  
  # set parallel options if required
  if (parallel) {
    paropts  <- list(.packages=c("phyloseq", "reshape2"))
  } else {
    paropts  <- NULL
  }
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive() && ! parallel, 'text', 'none'), .parallel=parallel, .paropts=paropts)
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}



Go_rare2 <- function(psIN, project, alpha_metrics, color, group, xlimit, plot = TRUE, parallel = FALSE, se = TRUE) {
  ps.rare <- calRare(psIN, alpha_metrics, rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
  summary(ps.rare)
  ps.summary <- ddply(ps.rare, c('Depth', 'Sample', 'Measure'), summarise, 
                      Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
  
  ps.summary$Sample <-  gsub("X", "", ps.summary$Sample); ps.summary$Sample
  ps.summ.verbose <- data.frame(merge(ps.summary, data.frame(sample_data(psIN)),by.x = 'Sample', by.y = 'row.names'))
  
  sample_variables(psIN)
  ps.summ.verbose$Measure
  
  
  # Plot
  plotlist <- list()
  for (am in alpha_metrics){
    ps.summ.verbose.sel <- subset(ps.summ.verbose, Measure == am)
    p <- ggplot(ps.summ.verbose.sel, mapping = aes_string(x = "Depth", y = "Alpha_diversity_mean", 
                                                   ymin = "Alpha_diversity_mean - Alpha_diversity_sd", 
                                                   ymax = "Alpha_diversity_mean + Alpha_diversity_sd", 
                                                   colour = color, group = group))+
      xlim(NA, xlimit) + theme_light() + geom_line() + #geom_pointrange(size = 0.1) +
      theme(legend.position="right", legend.text=element_text(size=8))+ guides(col = guide_legend(ncol = 2)) +
      ggtitle(sprintf("%s-Rarefaction curve", am )) + labs(x = "Sequence Sample Size", y = am)
    # + facet_wrap(facets = ~ StudyID, scales = 'free_y')
    #print(p)
    #plotlist[[length(plotlist)+1]] <- p 
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
  #multiplot(plotlist=plotlist, cols=cols, rows=rows)

}



