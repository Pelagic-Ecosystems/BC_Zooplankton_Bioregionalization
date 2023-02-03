# Goodness-of-clustering evaluation
# Patrick Pata
# Modified October 21, 2021
# 
# Set up goodness of clustering evaluation. This chunk contains a wrapper function 
# for two cluster validation stats from the optimclass and fpc packages. Indices 1-5 
# are evaluated against the species abundance table while indices 6-14 are evaluated 
# against the dissimilarity matrix. Also included with the 8 recommended indices in 
# the optpart package (Roberts et al. 2015) are 2 other indices. The within-cluster 
# sum of squares is the metric optimized in ward hierarchical clustering and k-means
# clustering. The Calinski-Harabaz criterion is a common cluster evaluation metric.
#
# Indices 11-14 are not used in the formal analysis because they have either very 
# similar goals as the other indices or the direction of optimal clustering is 
# unclear to me. These are supplied as reference for exploration and other indices 
# in the fpc package are not included here as well.

require(fpc)
require(optpart)

goc_index <- function(spemat, multiclust, distmat) {
  ceval <- matrix(NA, nrow = ncol(multiclust), ncol = 14)
  # names = list("indval","isamic","optimclass","tabdev","totchi","disdiam","partana","silhouette"))
  # include clusterstats indices
  pval = 0.000001 # Tichy recommends 10^-6
  
  for ( k in seq(ncol(multiclust)) ) {
    clust <- multiclust[,k]
    
    # 1. Indval (higher better) 
    # pobj <- optindval(spemat, clust, maxitr = 10, minsiz = 5) # seems to take longer than indval()
    pobj <- indval(spemat,clust)
    ceval[k,1] <- sum(pobj$indcls[pobj$pval < 0.05] )
    
    # 2. ISAMIC (higher better)
    pobj <- isamic(spemat, clust, sort = FALSE)
    ceval[k,2] <- sum(pobj) / length(pobj)
    
    # 3. Optimclass (higher better)
    tmpx <- mondo.fisher(comm=spemat,clust=clust)
    ceval[k,3] <- sum(tmpx <= pval) # number of significant taxa (across all clusters) - original
    
    # 4. TABDEV (lower better)
    pobj <- tabdev(spemat,clust)
    ceval[k,4] <- pobj$totdev
    
    # 5. TOTCHI (lower better)
    totchi <- 0
    for (c in seq(max(clust))) { # this should theoretically work
      pobj <- cca(spemat[clust==c,])
      totchi <- totchi+ pobj$tot.chi
    }
    ceval[k,5] <- totchi
    
    # Community-based evaluators
    
    # 6. DISDIAM (lower beter)
    pobj <- disdiam(clust, distmat ) # revised Feb 25, 2021 to use the original function
    ceval[k,6] <- pobj$mean
    
    # 7. PARTANA (higher better)
    pobj <- partana(clust, distmat)
    ceval[k,7] <- pobj$ratio
    
    # 8. Silhouette Width (higher better)
    pobj <- gensilwidth(clust, distmat)
    ceval[k,8] <- mean(pobj[,3])
    
    # 9. Add common indices from clusterstats() of fpc package
    cstat <- cluster.stats(distmat, clust, silhouette = FALSE)
    ceval[k,9] <- cstat$within.cluster.ss # lower better
    ceval[k,10] <- cstat$ch # higher better
    
    
    ceval[k,11] <- cstat$dunn # higher better
    ceval[k,12] <- cstat$wb.ratio # average within/average between... like partana so higher better? 
    ceval[k,13] <- cstat$pearsongamma # A version of hubert's gamma coefficient... higher better
    ceval[k,14] <- cstat$entropy # lower better? unclear from Melia 2007
  }
  
  ceval <- as.data.frame( cbind(c(2:max(multiclust)), ceval) )
  colnames(ceval) <- c("Num_Clusters","Indval","ISAMIC","Optimclass","TABDEV",
                       "TOTCHI","DISDIAM","PARTANA","Silhoutte",
                       "WithinClusterSS","CalinskiHarabaz",
                       "DunnIndex","WBRatio","PearsonGamma","Entropy") 
  
  return(ceval)
}

# --- Optimclass is Fisher test with pval used as "measure of positive fidelity of species to cluster"
mondo.fisher <- function(comm, clust) {
  rows <- ncol(comm)
  cols <- length(table(clust))
  res <- matrix(NA, nrow = rows, ncol = cols)
  for (i in 1:rows) {
    for (j in 1:cols) {
      
      # res[i, j] <- fisher.test(comm[, i] > 0, clust == j)$p.val
      x <- comm[,i]>0
      y <- clust == j
      z <- table(x,y)
      if(nrow(z)>1 & ncol(z)>1){
        res[i, j] <- fisher.test(x,y)$p.val
      } else {
        res[i, j] <- 2
      }
      
    }
  }
  res
}
