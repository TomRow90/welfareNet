#' Compare simulated etworks
#'
#' This function computes various performance measures to be used in simulation studies
#'
#'
#'
#' @param true true network
#' @param est estimated network
#' @param metric which performance measures to use
#' @param directed is the network directed or not
#'
#' @details
#' this is a function from https://osf.io/preprints/psyarxiv/4j3hf_v1 added to this package for convenience to be used
#' in simulation studies
#'
#' @export

compare.sim.networks <- function(true, est, metric = NULL, directed = FALSE){

  # Warning if input is not a matrix:
  if (!is.matrix(true) | !is.matrix(est)) stop("Input must be a weight matrix")

  if (is.null(metric)){
    metric <- c("true_positives", "true_negatives","false_positives","false_negatives","sensitivity", "signed_sensitivity", "sensitivity_top50",
                "sensitivity_top25", "sensitivity_top10", "specificity",
                "precision", "precision_top50", "precision_top25",
                "precision_top10", "jaccard_index", "correlation",
                "correlation_abs", "correlation_true", "bias", "bias_true",
                "centrality_Pearson_cor", "centrality_Kendall_cor",
                "centrality_top5", "centrality_top3", "centrality_top1")
  } else if(any(!(metric %in% c("true_positives", "true_negatives","false_positives","false_negatives","sensitivity", "signed_sensitivity", "sensitivity_top50",
                                "sensitivity_top25", "sensitivity_top10", "specificity",
                                "precision", "precision_top50", "precision_top25",
                                "precision_top10", "jaccard_index",
                                "correlation", "correlation_abs", "correlation_true",
                                "bias", "bias_true", "centrality_Pearson_cor",
                                "centrality_Kendall_cor", "centrality_top5",
                                "centrality_top3", "centrality_top1")))) {
    stop("metric invalid; needs to be 'true_positives', 'true_negatives','false_positives','false_negatives','sensitivity', 'signed_sensitivity', 'sensitivity_top50',
    'sensitivity_top25', 'sensitivity_top10', 'specificity', 'precision', 'precision_top50',
    'precision_top25', 'precision_top10', 'jaccard_index',
         'correlation', 'correlation_abs', 'correlation_true', 'bias', 'bias_true',
         'centrality_Pearson_cor', 'centrality_Kendall_cor', 'centrality_top5',
         'centrality_top3', 'centrality_top1'")
  }

  # Check if centrality needs to be computed:
  if(any(metric %in% c("centrality_Pearson_cor", "centrality_Kendall_cor", "centrality_top5", "centrality_top3", "centrality_top1"))){
    centTrue <- qgraph::centrality(true)
    centEst <- qgraph::centrality(est)
  }

  # Check if network is directed:
  if (directed){
    true <- c(true)
    est <- c(est)
  } else {
    true <- true[upper.tri(true, diag = FALSE)]
    est <- est[upper.tri(est, diag = FALSE)]
  }

  # Output list:
  out <- list()

  # True positives:
  truePos <- sum(est != 0 &  true != 0)

  # False pos:
  falsePos <- sum(est != 0 & true == 0)

  # True Neg:
  trueNeg <- sum(est == 0 & true == 0)

  # False Neg:
  falseNeg <- sum(est == 0 & true != 0)

  #true pos
  if("true_positives" %in% metric){
    out$true_positives <- truePos
  }

  #true neg
  if("true_negatives" %in% metric){
    out$true_negatives<- trueNeg
  }

  #False pos
  if("false_positives" %in% metric){
    out$false_positives <- falsePos
  }

  #False neg
  if("false_negatives" %in% metric){
    out$false_negatives <- falseNeg
  }

  # Sensitivity:
  if("sensitivity" %in% metric){
    out$sensitivity <- truePos / (truePos + falseNeg)
  }

  # Signed sensitivity:
  if("signed_sensitivity" %in% metric){
    truePosSigned <- sum(est != 0 &  true != 0 & sign(est) == sign(true))
    out$signed_sensitivity <- truePosSigned / (truePos + falseNeg)
  }

  # Sensitivity top 50%:
  if("sensitivity_top50" %in% metric){
    top50 <- which(abs(true) > median(abs(true[true!=0])))
    out[["sensitivity_top50"]] <- sum(est[top50]!=0 & true[top50] != 0) / sum(true[top50] != 0)
  }

  # Sensitivity top 25%:
  if("sensitivity_top25" %in% metric){
    top25 <- which(abs(true) > quantile(abs(true[true!=0]), 0.75))
    out[["sensitivity_top25"]] <- sum(est[top25]!=0 & true[top25] != 0) / sum(true[top25] != 0)
  }

  # Sensitivity top 10%:
  if("sensitivity_top10" %in% metric){
    top10 <- which(abs(true) > quantile(abs(true[true!=0]), 0.90))
    out[["sensitivity_top10"]] <- sum(est[top10]!=0 & true[top10] != 0) / sum(true[top10] != 0)
  }

  # Specificity:
  if("specificity" %in% metric){
    out$specificity <- trueNeg / (trueNeg + falsePos)
  }

  # Precision (1 - FDR):
  if("precision" %in% metric){
    out$precision <- truePos / (falsePos + truePos)
  }

  # precision top 50% (of estimated edges):
  if("precision_top50" %in% metric){
    top50 <- which(abs(est) > median(abs(est[est!=0])))
    out[["precision_top50"]] <- sum(est[top50]!=0 & true[top50] != 0) / sum(est[top50] != 0)
  }

  # precision top 25%:
  if("precision_top25" %in% metric){
    top25 <- which(abs(est) > quantile(abs(est[est!=0]), 0.75))
    out[["precision_top25"]] <- sum(est[top25]!=0 & true[top25] != 0) / sum(est[top25] != 0)
  }

  # precision top 10%:
  if("precision_top10" %in% metric){
    top10 <- which(abs(est) > quantile(abs(est[est!=0]), 0.90))
    out[["precision_top10"]] <- sum(est[top10]!=0 & true[top10] != 0) / sum(est[top10] != 0)
  }

  # Jaccard index:
  if("jaccard_index" %in% metric){
    inter <- sum(true & est)
    union <- sum(true | est)

    if (union == 0) { # avoid division by zero
      out$jaccard_index <- 1
    } else {
      # Calculate Jaccard Index
      out$jaccard_index <- inter / union
    }
  }

  # correlation:
  if("correlation" %in% metric){
    out$correlation <- cor(est, true)
  }

  # correlation between absolute edges:
  if("correlation_abs" %in% metric){
    out$abs_correlation <- cor(abs(est), abs(true))
  }

  # correlation between true edge weights:
  if("correlation_true" %in% metric){
    if (truePos > 0){
      trueEdges <- est != 0 & true != 0
      out$correlation_true <- cor(est[trueEdges], true[trueEdges])
    } else {
      out$correlation_true <- NA
    }
  }

  # average bias between edge weights:
  if("bias" %in% metric){
    out$bias <- mean(abs(est-true), na.rm = TRUE)
  }

  # average bias between true edge weights:
  if("bias_true" %in% metric){
    if (truePos > 0){
      trueEdges <- est != 0 & true != 0
      out$bias_true <- mean(abs(est[trueEdges] - true[trueEdges]))
    } else {
      out$bias_true <- NA
    }
  }

  # Pearson correlations:
  if("centrality_Pearson_cor" %in% metric){
    if(directed){
      out$in_strength_correlation <- cor(centTrue$InDegree, centEst$InDegree)
      out$out_strength_correlation <- cor(centTrue$OutDegree, centEst$OutDegree)
    } else {
      out$strength_correlation <- cor(centTrue$OutDegree, centEst$OutDegree)
    }
    out$closeness_correlation <- cor(centTrue$Closeness, centEst$Closeness)
    out$betweenness_correlation <- cor(centTrue$Betweenness, centEst$Betweenness)
  }

  # Kendall correlations:
  if("centrality_Kendall_cor" %in% metric){
    if(directed){
      out$in_strength_correlation_kendall <- cor(centTrue$InDegree, centEst$InDegree, method = "kendall")
      out$out_strength_correlation_kendall <- cor(centTrue$OutDegree, centEst$OutDegree, method = "kendall")
    } else{
      out$strength_correlation_kendall <- cor(centTrue$OutDegree, centEst$OutDegree, method = "kendall")
    }
    out$closeness_correlation_kendall <- cor(centTrue$Closeness, centEst$Closeness, method = "kendall")
    out$betweenness_correlation_kendall <- cor(centTrue$Betweenness, centEst$Betweenness, method = "kendall")

  }

  # Centrality top 5:
  if("centrality_top5" %in% metric){
    if(directed){
      out$in_strength_top5 <- mean(order(centEst$InDegree, decreasing = TRUE)[1:5] %in% order(centTrue$InDegree, decreasing = TRUE)[1:5])
      out$out_strength_top5 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1:5] %in% order(centTrue$OutDegree, decreasing = TRUE)[1:5])
    } else{
      out$strength_top5 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1:5] %in% order(centTrue$OutDegree, decreasing = TRUE)[1:5])
    }
    out$closeness_top5 <-  mean(order(centEst$Closeness, decreasing = TRUE)[1:5] %in% order(centTrue$Closeness, decreasing = TRUE)[1:5])
    out$betweenness_top5 <-  mean(order(centEst$Betweenness, decreasing = TRUE)[1:5] %in% order(centTrue$Betweenness, decreasing = TRUE)[1:5])
  }

  # Centrality top 3:
  if("centrality_top3" %in% metric){
    if(directed){
      out$in_strength_top3 <- mean(order(centEst$InDegree, decreasing = TRUE)[1:3] %in% order(centTrue$InDegree, decreasing = TRUE)[1:3])
      out$out_strength_top3 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1:3] %in% order(centTrue$OutDegree, decreasing = TRUE)[1:3])
    } else{
      out$strength_top3 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1:3] %in% order(centTrue$OutDegree, decreasing = TRUE)[1:3])
    }
    out$closeness_top3 <-  mean(order(centEst$Closeness, decreasing = TRUE)[1:3] %in% order(centTrue$Closeness, decreasing = TRUE)[1:3])
    out$betweenness_top3 <-  mean(order(centEst$Betweenness, decreasing = TRUE)[1:3] %in% order(centTrue$Betweenness, decreasing = TRUE)[1:3])

  }

  # Centrality top 1:
  if("centrality_top1" %in% metric){
    if(directed){
      out$in_strength_top1 <- mean(order(centEst$InDegree, decreasing = TRUE)[1] %in% order(centTrue$InDegree, decreasing = TRUE)[1])
      out$out_strength_top1 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1] %in% order(centTrue$OutDegree, decreasing = TRUE)[1])
    }else{
      out$strength_top1 <- mean(order(centEst$OutDegree, decreasing = TRUE)[1] %in% order(centTrue$OutDegree, decreasing = TRUE)[1])
    }
    out$closeness_top1 <- mean(order(centEst$Closeness, decreasing = TRUE)[1] %in% order(centTrue$Closeness, decreasing = TRUE)[1])
    out$betweenness_top1 <- mean(order(centEst$Betweenness, decreasing = TRUE)[1] %in% order(centTrue$Betweenness, decreasing = TRUE)[1])
  }

  # Return evaluation metrics:
  return(out)
}



#' Generate GGM function. From (bootnet: https://github.com/SachaEpskamp/bootnet/tree/master). Used for simulation studies.
#'
#' @param Nvar number of nodes
#' @param p rewiring probability
#' @param nei average neighbourhood of each node
#' @param parRange range of partial correlations
#' @param constant default = 1.5
#' @param propPositive Proportion of positive partial correlations
#' @param clusters number of clusters ig graph = "cluster"
#' @param graph one of smallworld, random, scalefree, hub, or cluster
#'
#' @return partial correlation matrix
#' @export



genGGM <- function(
    Nvar,
    p = 0, # Rewiring probability if graph = "smallworld" or "cluster", or connection probability if graph = "random". If cluster, can add multiple p's for each cluster, e.g., "c(.1, .5)"
    nei = 1,
    parRange = c(0.5,1),
    constant = 1.5,
    propPositive = 0.5,
    clusters = NULL, #number of clusters if graph = "cluster"
    graph = c("smallworld","random", "scalefree", "hub", "cluster")
){
  graph <- match.arg(graph)


  ## Approach from
  # Yin, J., & Li, H. (2011). A sparse conditional gaussian graphical model for analysis of genetical genomics data. The annals of applied statistics, 5(4), 2630.

  # Simulate graph structure:
  if (graph == "smallworld"){
    # Watts Strogatz small-world
    trueKappa <- as.matrix(igraph::as_adjacency_matrix(igraph::sample_smallworld(1,Nvar,nei,p)))
  } else if (graph == "random"){
    # Ranodm network:
    trueKappa <- as.matrix(igraph::as_adjacency_matrix(igraph::erdos.renyi.game(Nvar, p)))
  } else if (graph == "scalefree") {
    if(!requireNamespace("BDgraph")) stop("'BDgraph' package needs to be installed.")

    trueKappa <- BDgraph::bdgraph.sim(p = Nvar, graph = "scale-free")$G
  } else if (graph == "hub") {
    if(!requireNamespace("BDgraph")) stop("'BDgraph' package needs to be installed.")

    trueKappa <- BDgraph::bdgraph.sim(p = Nvar, graph = "hub")$G
    class(trueKappa) <- "matrix"
  } else if (graph == "cluster") {
    if(!requireNamespace("BDgraph")) stop("'BDgraph' package needs to be installed.")

    trueKappa <-  BDgraph::bdgraph.sim(p = Nvar, graph = "cluster", prob = p, class = clusters)$G #can be
    class(trueKappa) <- "matrix"
  }

  # Make edges negative and add weights:
  trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] * sample(c(-1,1),sum(upper.tri(trueKappa)),TRUE,prob=c(propPositive,1-propPositive)) *
    runif(sum(upper.tri(trueKappa)), min(parRange ),max(parRange ))

  # Symmetrize:
  trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]

  # Make pos def:
  diag(trueKappa) <- constant * rowSums(abs(trueKappa))
  diag(trueKappa) <- ifelse(diag(trueKappa)==0,1,diag(trueKappa))
  trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
  trueKappa <- (trueKappa + t(trueKappa)) / 2

  return(as.matrix(qgraph::wi2net(trueKappa)))
}




#' simulate data from a given GGM (again taken from bootnet: https://github.com/SachaEpskamp/bootnet/tree/master)
#'
#' @param trueNet partial correlation matrix
#' @param sampleSize sample size
#'
#' @return data as matrix
#' @export


sim.data <- function(trueNet, sampleSize) {

  #simulate date
  graph <- trueNet
  intercepts <- rep(0,ncol(graph))

  # standardize:
  if (!all(diag(graph) == 0 | diag(graph) == 1)){
    graph <- cov2cor(graph)
  }

  # Remove diag:
  diag(graph) <- 0

  # Generate data:
  # True sigma:
  if (any(eigen(diag(ncol(graph)) - graph)$values < 0)){
    stop("Precision matrix is not positive semi-definite")
  }

  Sigma <- cov2cor(solve(diag(ncol(graph)) - graph))


  # Generate data:
  data <- mvtnorm::rmvnorm(sampleSize, sigma = Sigma)

}




#' find approximate pcor given an adjacency matrix and parameter range (again taken from bootnet: https://github.com/SachaEpskamp/bootnet/tree/master)
#'
#' @param adj the signed adjacency matrix
#' @param parRange the range of absolute partial correlations to start from (actual pcors will likely be lower than this)
#'
#' @return true.net = the partial correlation matrix
#' @export

adj.to.pcor <- function(adj, parRange=c(0.5,1)) {

  trueKappa <- adj

  # Make edges adj# Make edges negative and add weights:
  trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] *
    runif(sum(upper.tri(trueKappa)), min(parRange ),max(parRange ))

  # Symmetrize:
  trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]

  # Make pos def:
  diag(trueKappa) <- 1.5 * rowSums(abs(trueKappa))
  diag(trueKappa) <- ifelse(diag(trueKappa)==0,1,diag(trueKappa))
  trueKappa <- trueKappa/diag(trueKappa)[row(trueKappa)]
  trueKappa <- (trueKappa + t(trueKappa)) / 2

  return(true.net = trueKappa)

}






