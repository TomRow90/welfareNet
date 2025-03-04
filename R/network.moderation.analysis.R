#' Perform network moderation analysis of zero-order networks with a permutation test
#'
#' This function performs network moderation analysis. That is, it performs a permutation test to compare two networks
#' on various global, nodal, and edge level characteristics. It uses zero-order distance correlations.
#'
#'
#'
#' @param data1 the dataset as matrix or data frame for group 1
#' @param data2 the dataset as matrix or data frame for group 2
#' @param net1 the estimate.dcor.network object of group 1
#' @param net2 the estimated.dcor.network object of group 2
#' @param permutations the number of permutations to perform. Default is 1000.
#'
#' @return dissimilarity.results = network.dissimilarity estimates, difference, P and S values
#' @return global.results = global network topology estimates, difference, P and S values
#' @return strength.results = Node strength estimates, difference, P and S values
#' @return edge.differences = Edge weights, difference, P and S values
#' @return perms = A list of all permuted null distributions for each network statistic
#'
#' @details
#' Note that in the returned objects/results, both the P and S value are supplied. For dissimilarity
#' the P value is one-sided, for all other network descriptors it is two-sided. The S value is also
#' known as the surprisal value and is simply calculated as -log2(P value). It converts the P value
#' into binary 'bits' of information against the test hypothesis and can be interpreted in the context
#' of flipping a coin. That is, a P value of 0.05 as an S value is 4.32, rounded to the nearest integer
#' we can say that assuming the null hypothesis (and all other assumptions) are true, then the result
#' is as surprising as getting 4 heads in a row in a series of coin tosses. See Rafi and Greenland (2020)
#' for more information and explanation about the surprisal value (https://doi.org/10.1186/s12874-020-01105-9).
#'
#' @export


network.moderation.analysis <- function(data1, data2, net1, net2, permutations) {

  #get node and edge number
  nodes <- ncol(data1)
  edges <- nodes * (nodes - 1) / 2

  #get empirical networks
  net1 <- net1$saturated.net
  net2 <- net2$saturated.net

  #calculate empirical topological descriptors
  emperical.frob.norm <- sF(net1, net2, nodes)

  emperical.weighted.density.diff <- weighted.density(net1) - weighted.density(net2)

  emperical.weighted.cc.diff <- global.weighted.cc(net1) - global.weighted.cc(net2)

  emperical.strength.diff <- strength(net1) - strength(net2)

  emperical.edge.diff <- net1 - net2
  emperical.edge.diff <- emperical.edge.diff[lower.tri(emperical.edge.diff)]


  #permutation test function. Adapted from Van Borkulo et al (2023)
  perm.test <- function(data1, data2) {

    nodes <- ncol(data1)
    edges <- nodes*(nodes - 1) / 2

    combined.data <- rbind(data1,data2)

    n.obs1 <- nrow(data1)
    n.obs2 <- nrow(data2)

    #generate permuted dataset
    s <- sample(1:(n.obs1 + n.obs2),n.obs1,replace=FALSE)
    x1perm <- combined.data[s,]
    x2perm <- combined.data[-s, ]

    #estimate networks
    perm.net1 <- estimate.dcor.network(x1perm, select = "saturated")
    perm.net1 <- perm.net1$saturated.net

    perm.net2 <- estimate.dcor.network(x2perm, select="saturated")
    perm.net2 <- perm.net2$saturated.net

    #global descriptors
    perm.frob.norm <- sF(perm.net1, perm.net2, nodes)

    perm.weighted.density.diff <- weighted.density(perm.net1) - weighted.density(perm.net2)

    perm.weighted.cc.diff <- global.weighted.cc(perm.net1) - global.weighted.cc(perm.net2)

    #centrality differences
    strength.diff <- strength(perm.net1) - strength(perm.net2)

    #edge differences
    edge.diff <- perm.net1 - perm.net2
    edge.diff <- edge.diff[lower.tri(edge.diff)]


    return(list(frob.norm = perm.frob.norm,
                weighted.density.diff = perm.weighted.density.diff,
                weighted.cc.diff = perm.weighted.cc.diff,
                strength.diff = strength.diff,
                edge.diff = edge.diff))
  }



  perms <- pbapply::pbreplicate(perm.test(data1=data1, data2=data2),
                                n=permutations, simplify=FALSE)

  #extract permutations
  frob.norm.perms <- extract.property(perms, "frob.norm")

  weighted.density.perms <- extract.property(perms, "weighted.density.diff")

  weighted.cc.perms <- extract.property(perms, "weighted.cc.diff")

  strength.perms <- extract.property(perms, "strength.diff")

  edge.diff.perms <- extract.property(perms, "edge.diff")


  #calculate P values for all descriptors
  frob.norm.p <- one.sided.p.function(frob.norm.perms, emperical.frob.norm)

  weighted.density.p <- two.sided.p.function(weighted.density.perms, emperical.weighted.density.diff)

  weighted.cc.p <- two.sided.p.function(weighted.cc.perms, emperical.weighted.cc.diff)

  ####### produce output objects############

  #dissimilarities
  network.dissimilarities <- as.data.frame(matrix(0, nrow=1, ncol=4))
  colnames(network.dissimilarities) <- c("Dissimilarity Measure","Value","P value", "S value")

  network.dissimilarities$'Dissimilarity Measure' <- c("Frobenius Norm Dissimilarity")

  network.dissimilarities$Value <- emperical.frob.norm

  network.dissimilarities$`P value` <- frob.norm.p

  network.dissimilarities$`S value` <- -log2(network.dissimilarities$`P value`)

  #global descriptors
  global.results <- as.data.frame(matrix(0, nrow=2, ncol=6))
  colnames(global.results) <- c("Descriptor", "Net1", "Net2", "Difference",
                                "P value", "S value")

  global.results$Descriptor <- c("Weighted density",
                                 "Weighted Clustering")

  global.results$Net1[1] <- weighted.density(net1)
  global.results$Net2[1] <- weighted.density(net2)

  global.results$Net1[2] <- global.weighted.cc(net1)
  global.results$Net2[2] <- global.weighted.cc(net2)


  global.results$Difference <- c(emperical.weighted.density.diff,
                                 emperical.weighted.cc.diff)

  global.results$`P value` <- c(weighted.density.p,
                                weighted.cc.p)

  global.results$`S value` <- -log2(global.results$`P value`)

  ########. Node level descriptors #########

  #strength
  strength.differences <- as.data.frame(matrix(0, nrow = nodes, ncol = 6))
  colnames(strength.differences) <- c("Node", "Net1","Net2","Difference", "P value", "S value")
  strength.differences$Node <- colnames(data1)

  strength.differences$Net1 <- strength(net1)
  strength.differences$Net2 <- strength(net2)

  strength.differences$`Difference` <- emperical.strength.diff

  for(i in 1:nodes) {

    strength.differences$`P value`[i] <- two.sided.p.function(strength.perms[i,], strength.differences$Difference[i])

  }

  strength.differences$`S value` <- -log2(strength.differences$`P value`)

  ########### edge level ##################

  edge.differences <- as.data.frame(matrix(0, nrow = edges, ncol = 6))
  colnames(edge.differences) <- c("Edge", "Net1", "Net2", "Difference", "P value", "S value")
  pairs <- combn(colnames(data1), 2, FUN = function(x) paste(x, collapse = "-"), simplify = FALSE)
  pairs <- as.character(pairs)
  edge.differences$Edge <- pairs

  edge.differences$Net1 <- net1[lower.tri(net1)]
  edge.differences$Net2 <- net2[lower.tri(net2)]

  edge.differences$Difference <- emperical.edge.diff

  for(i in 1:edges) {

    edge.differences$`P value`[i] <- two.sided.p.function(edge.diff.perms[i,], edge.differences$Difference[i])

  }

  edge.differences$`S value` <- -log2(edge.differences$`P value`)

  #create permutation list
  perms <- list()

  perms[[1]] <- frob.norm.perms
  perms[[2]] <- weighted.density.perms
  perms[[3]] <- weighted.cc.perms
  perms[[4]] <- strength.perms
  perms[[5]] <- edge.diff.perms

  names(perms) <- c("Frobenius Norm Dissimarlity","Weighted Density",
                    "Weighted Clustering Coefficient",
                    "Node Strength","Edge Weight Differences")


  return(list(dissimilarity.results = network.dissimilarities,
              global.results = global.results,
              strength.results = strength.differences,
              edge.differences = edge.differences,
              perms = perms))

}
