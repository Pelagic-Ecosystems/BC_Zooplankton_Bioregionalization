'p_indvalScatter_25' <- function(indspe, plot = TRUE, threshold = 0.25) {
  # Created P.Pata - June 14, 2020
  # Received the output of the indicator species analysis ( indval() ).
  # Revised - February 4, 2021 to only plot and return taxa which are >= 0.25 indval (i.e., 0.5 fidelity and 0.5% specificity)
  require(ggplot2)
  require(tibble)
  require(ggrepel)
  # Apply a correction for multiple testing [http://biol09.biol.umontreal.ca/PLcourses/Indicator_species.pdf]
  prob.corrected = p.adjust(indspe$pval, "holm") # increases number of insignificant species
  
  # Un/comment to see details of which taxa did not make the cut
  # # check which species have very low significance
  # nsigsp <- rownames_to_column(as.data.frame(prob.corrected ),"Species") %>% 
  #   rename(pval = "prob.corrected") %>% filter(pval > 0.05)
  # plot(indspe$pval,prob.corrected)
  
  # Isolate significant indicator value species and frequency
  indspe.iv <- rownames_to_column(as.data.frame( indspe$indcls ),"Species") %>%  
    rename(indval = "indspe$indcls")
  indspe.clust <- rownames_to_column(as.data.frame( indspe$maxcls ),"Species") %>%  
    rename(cluster = "indspe$maxcls")
  sigsp <- rownames_to_column(as.data.frame(prob.corrected ),"Species") %>% 
    rename(pval = "prob.corrected") %>% filter(pval < 0.05) 
  indspe.freq <- rowSums(as.data.frame(indspe$relfrq)) / ncol(indspe$relfrq) * 100
  indspe.freq <- rownames_to_column(as.data.frame( indspe.freq ),"Species")%>%  
    rename(relfreq = "indspe.freq")
  sigsp <- left_join(sigsp, indspe.iv, by = "Species")
  sigsp <- left_join(sigsp, indspe.freq, by = "Species")
  sigsp <- left_join(sigsp, indspe.clust, by = "Species")
  sigsp <- sigsp %>%
    filter(indval >= threshold)
  
  if (plot == TRUE) {
    # Scatterplot of highest indval clusters
    p <- ggplot(sigsp, aes(x=relfreq, y=indval, label = Species)) + 
      geom_point(data = sigsp, size=4, alpha = 0.8, shape = 16, 
                 aes(colour=as.factor(cluster))) +
      theme(legend.title = element_text(color = 'black', size = 12, face = 'bold'),
            legend.text = element_text(size = 12),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 12)) +
      geom_text_repel(data = sigsp, size=4, hjust=0, 
                      nudge_x = 0.5, fontface="italic",
                      aes(color = as.factor(cluster))) +
      xlab("% Presence in samples") + ylab("Indicator Value") +
      scale_colour_discrete(name = "Clusters") +
      theme_bw()
    plot(p)
  }
  
  return(sigsp)
}