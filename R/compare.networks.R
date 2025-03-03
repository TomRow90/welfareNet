#' Perform network moderation analysis with a permutation test
#'
#' This function performs network moderation analysis. That is, it performs a permutation test to compare two networks
#' on various global, nodal, and edge level characteristics.
#'
#'
#'
#' @param data1 the dataset as matrix or data frame for group 1
#' @param data2 the dataset as matrix or data frame for group 2
#' @param net1 the estimate.network object of group 1 from the estimate.network function
#' @param net2 the estimated.network object of group 2 from the estimate.network function
#' @param permutations the number of permutations to perform. Default is 1000.
#'
#' @return dissimilarity.results = network.dissimilarity estimates, difference, P and S values
#' @return global.results = global network topology estimates, difference, P and S values
#' @return degree.results = Node degree estimates, difference, P and S values
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


compare.networks <- function(data1, data2, net1, net2, permutations) {

  #get node and edge number
  nodes <- ncol(data1)
  edges <- nodes * (nodes - 1) / 2

  #calculate emperical topological descriptors
  emperical.ham.d <- hamming.distance(net1, net2)

  emperical.jaccard <- jaccard.similarity(net1, net2)

  emperical.frob.norm <- sF(net1, net2, nodes)

  emperical.density.diff <- unweighted.density(net1) - unweighted.density(net2)

  emperical.weighted.density.diff <- weighted.density(net1) - weighted.density(net2)

  emperical.pos.diff <- proportion.pos(net1) - proportion.pos(net2)

  emperical.neg.diff <-  proportion.neg(net1) - proportion.neg(net2)

  emperical.cc.diff <- global.cc(net1) - global.cc(net2)

  emperical.weighted.cc.diff <- global.weighted.cc(net1) - global.weighted.cc(net2)

  emperical.degree.diff <- degree(net1) - degree(net2)

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
    perm.net1 <- estimate.network(x1perm, select = "saturated")
    perm.net1 <- perm.net1$selected.net

    perm.net2 <- estimate.network(x2perm, select="saturated")
    perm.net2 <- perm.net2$selected.net

    #global descriptors
    perm.ham.d <- hamming.distance(perm.net1, perm.net2)

    perm.jaccard <- jaccard.similarity(perm.net1, perm.net2)

    perm.frob.norm <- sF(perm.net1, perm.net2, nodes)

    perm.density.diff <- unweighted.density(perm.net1) - unweighted.density(perm.net2)

    perm.weighted.density.diff <- weighted.density(perm.net1) - weighted.density(perm.net2)

    perm.pos.diff <- proportion.pos(perm.net1) - proportion.pos(perm.net2)

    perm.neg.diff <-  proportion.neg(perm.net1) - proportion.neg(perm.net2)

    perm.cc.diff <- global.cc(perm.net1) - global.cc(perm.net2)

    perm.weighted.cc.diff <- global.weighted.cc(perm.net1) - global.weighted.cc(perm.net2)

    #centrality differences
    degree.diff <- degree(perm.net1) - degree(perm.net2)

    strength.diff <- strength(perm.net1) - strength(perm.net2)

    #edge differences
    edge.diff <- perm.net1 - perm.net2
    edge.diff <- edge.diff[lower.tri(edge.diff)]


    return(list(ham.d = perm.ham.d,
                frob.norm = perm.frob.norm,
                jaccard = perm.jaccard,
                density.diff = perm.density.diff,
                weighted.density.diff = perm.weighted.density.diff,
                positive.diff = perm.pos.diff,
                negative.diff = perm.neg.diff,
                cc.diff = perm.cc.diff,
                weighted.cc.diff = perm.weighted.cc.diff,
                degree.diff = degree.diff,
                strength.diff = strength.diff,
                edge.diff = edge.diff))

  }



  perms <- pbapply::pbreplicate(perm.test(data1=data1, data2=data2),
                               n=permutations, simplify=FALSE)

  #extract permutations
  ham.d.perms <- extract.property(perms, "ham.d")

  jaccard.perms <- extract.property(perms, "jaccard")

  frob.norm.perms <- extract.property(perms, "frob.norm")

  density.perms <- extract.property(perms, "density.diff")

  weighted.density.perms <- extract.property(perms, "weighted.density.diff")

  positive.perms <- extract.property(perms, "positive.diff")

  negative.perms <- extract.property(perms, "negative.diff")

  cc.perms <- extract.property(perms, "cc.diff")

  weighted.cc.perms <- extract.property(perms, "weighted.cc.diff")

  degree.perms <- extract.property(perms, "degree.diff")

  strength.perms <- extract.property(perms, "strength.diff")

  edge.diff.perms <- extract.property(perms, "edge.diff")


  #calculate P values for all descriptors

  ham.d.p <- one.sided.p.function(ham.d.perms, emperical.ham.d)

  jaccard.p <- one.sided.p.function(jaccard.perms, emperical.jaccard)

  frob.norm.p <- one.sided.p.function(frob.norm.perms, emperical.frob.norm)

  density.p <- two.sided.p.function(density.perms, emperical.density.diff)

  weighted.density.p <- two.sided.p.function(weighted.density.perms, emperical.weighted.density.diff)

  positive.p <- two.sided.p.function(positive.perms, emperical.pos.diff)

  negative.p <- two.sided.p.function(negative.perms, emperical.neg.diff)

  cc.p <- two.sided.p.function(cc.perms, emperical.cc.diff)

  weighted.cc.p <- two.sided.p.function(weighted.cc.perms, emperical.weighted.cc.diff)

  ####### produce output objects############

  #dissimilarities
  network.dissimilarities <- as.data.frame(matrix(0, nrow=3, ncol=4))
  colnames(network.dissimilarities) <- c("Dissimilarity Measure","Value","P value", "S value")

  network.dissimilarities$'Dissimilarity Measure' <- c("Hamming Dissimilarity", "Jaccard Dissimilarity","Frobenius Norm Dissimilarity")

  network.dissimilarities$Value <- c(emperical.ham.d, emperical.jaccard, emperical.frob.norm)

  network.dissimilarities$`P value` <- c(ham.d.p, jaccard.p, frob.norm.p)

  network.dissimilarities$`S value` <- -log2(network.dissimilarities$`P value`)

  #global descriptors
  global.results <- as.data.frame(matrix(0, nrow=6, ncol=6))
  colnames(global.results) <- c("Descriptor", "Net1", "Net2", "Difference",
                                "P value", "S value")

  global.results$Descriptor <- c("Density", "Weighted density",
                                 "Positive density", "Negative density",
                                 "Clustering", "Weighted Clustering")

  global.results$Net1[1] <- unweighted.density(net1)
  global.results$Net2[1] <- unweighted.density(net2)

  global.results$Net1[2] <- weighted.density(net1)
  global.results$Net2[2] <- weighted.density(net2)

  global.results$Net1[3] <- proportion.pos(net1)
  global.results$Net2[3] <- proportion.pos(net2)

  global.results$Net1[4] <- proportion.neg(net1)
  global.results$Net2[4] <- proportion.neg(net2)

  global.results$Net1[5] <- global.cc(net1)
  global.results$Net2[5] <- global.cc(net2)

  global.results$Net1[6] <- global.weighted.cc(net1)
  global.results$Net2[6] <- global.weighted.cc(net2)


  global.results$Difference <- c(emperical.density.diff, emperical.weighted.density.diff,
                                 emperical.pos.diff, emperical.neg.diff,
                                 emperical.cc.diff, emperical.weighted.cc.diff)

  global.results$`P value` <- c(density.p, weighted.density.p,
                                positive.p, negative.p,
                                cc.p, weighted.cc.p)

  global.results$`S value` <- -log2(global.results$`P value`)

  ########. Node level descriptors #########

  #degree
  degree.differences <- as.data.frame(matrix(0, nrow = nodes, ncol = 6))
  colnames(degree.differences) <- c("Node", "Net1","Net2","Difference", "P value", "S value")
  degree.differences$Node <- colnames(data1)

  degree.differences$Net1 <- degree(net1)
  degree.differences$Net2 <- degree(net2)

  degree.differences$`Difference` <- emperical.degree.diff

  for(i in 1:nodes) {

    degree.differences$`P value`[i] <- two.sided.p.function(degree.perms[i,], degree.differences$`Difference`[i])

  }

  degree.differences$`S value` <- -log2(degree.differences$`P value`)

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

  #signed weighted edges

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

  perms[[1]] <- ham.d.perms
  perms[[2]] <- jaccard.perms
  perms[[3]] <- frob.norm.perms
  perms[[4]] <- density.perms
  perms[[5]] <- weighted.density.perms
  perms[[6]] <- positive.perms
  perms[[7]] <- negative.perms
  perms[[8]] <- cc.perms
  perms[[9]] <- weighted.cc.perms
  perms[[10]] <- degree.perms
  perms[[11]] <- strength.perms
  perms[[12]] <- edge.diff.perms

  names(perms) <- c("Hamming Dissimarlity","Jaccard Dissimilarity","Frobenius Norm Dissimarlity",
                    "Unweighted Density","Weighted Density","Proportion Positive Edges",
                    "Proportion Negative Edges","Clustering Coefficient","Weighted Clustering Coefficient",
                    "Node Degree","Node Strength","Edge Weight Differences")




  return(list(dissimilarity.results = network.dissimilarities,
              global.results = global.results,
              degree.results = degree.differences,
              strength.results = strength.differences,
              edge.differences = edge.differences,
              perms = perms))

}
