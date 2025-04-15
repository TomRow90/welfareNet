#' Perform network comparison analysis of low order conditional unweighted networks with a permutation test
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


qp.mgm.comparison <- function(data1, data2, net1, net2, type, permutations) {

  #get node and edge number
  nodes <- ncol(data1)
  edges <- nodes * (nodes - 1) / 2

  #get empirical networks
  net1 <- net1$structure
  net2 <- net2$structure

  #calculate empirical topological descriptors
  emperical.ham.d <- hamming.distance(net1, net2)

  emperical.density.diff <- unweighted.density(net1) - unweighted.density(net2)

  emperical.cc.diff <- global.cc(net1) - global.cc(net2)

  emperical.degree.diff <- degree(net1) - degree(net2)


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
    perm.net1 <- estimate.q1.mgm(x1perm, alpha=0.1, type=type)
    perm.net1 <- perm.net1$structure

    perm.net2 <- estimate.q1.mgm(x2perm, alpha=0.1, type=type)
    perm.net2 <- perm.net2$structure

    #global descriptors
    perm.ham.d <- hamming.distance(perm.net1, perm.net2)

    perm.density.diff <- unweighted.density(perm.net1) - unweighted.density(perm.net2)

    perm.cc.diff <- global.cc(perm.net1) - global.cc(perm.net2)

    #centrality differences
    degree.diff <- degree(perm.net1) - degree(perm.net2)


    return(list(ham.d = perm.ham.d,
                density.diff = perm.density.diff,
                cc.diff = perm.cc.diff,
                degree.diff = degree.diff))
  }



  perms <- pbapply::pbreplicate(perm.test(data1=data1, data2=data2),
                                n=permutations, simplify=FALSE)

  #extract permutations
  ham.d.perms <- extract.property(perms, "ham.d")

  density.perms <- extract.property(perms, "density.diff")

  cc.perms <- extract.property(perms, "cc.diff")

  degree.perms <- extract.property(perms, "degree.diff")


  #calculate P values for all descriptors
  ham.d.p <- one.sided.p.function(ham.d.perms, emperical.ham.d)

  density.p <- two.sided.p.function(density.perms, emperical.density.diff)

  cc.p <- two.sided.p.function(cc.perms, emperical.cc.diff)

  ####### produce output objects############

  #dissimilarities
  network.dissimilarities <- as.data.frame(matrix(0, nrow=1, ncol=4))
  colnames(network.dissimilarities) <- c("Dissimilarity Measure","Value","P value", "S value")

  network.dissimilarities$'Dissimilarity Measure' <- c("Hamming Distance")

  network.dissimilarities$Value <- emperical.ham.d

  network.dissimilarities$`P value` <- ham.d.p

  network.dissimilarities$`S value` <- -log2(network.dissimilarities$`P value`)

  #global descriptors
  global.results <- as.data.frame(matrix(0, nrow=2, ncol=6))
  colnames(global.results) <- c("Descriptor", "Net1", "Net2", "Difference",
                                "P value", "S value")

  global.results$Descriptor <- c("Density",
                                 "Clustering")

  global.results$Net1[1] <- unweighted.density(net1)
  global.results$Net2[1] <- unweighted.density(net2)

  global.results$Net1[2] <- global.cc(net1)
  global.results$Net2[2] <- global.cc(net2)


  global.results$Difference <- c(emperical.density.diff,
                                 emperical.cc.diff)

  global.results$`P value` <- c(density.p,
                                cc.p)

  global.results$`S value` <- -log2(global.results$`P value`)

  ########. Node level descriptors #########

  #strength
  degree.differences <- as.data.frame(matrix(0, nrow = nodes, ncol = 6))
  colnames(degree.differences) <- c("Node", "Net1","Net2","Difference", "P value", "S value")
  degree.differences$Node <- colnames(data1)

  degree.differences$Net1 <- degree(net1)
  degree.differences$Net2 <- degree(net2)

  degree.differences$`Difference` <- emperical.degree.diff

  for(i in 1:nodes) {

    degree.differences$`P value`[i] <- two.sided.p.function(degree.perms[i,], degree.differences$Difference[i])

  }

  degree.differences$`S value` <- -log2(degree.differences$`P value`)



  #create permutation list
  perms <- list()

  perms[[1]] <- ham.d.perms
  perms[[2]] <- density.perms
  perms[[3]] <- cc.perms
  perms[[4]] <- degree.perms

  names(perms) <- c("Hamming Distance Dissimarlity","Density",
                    "Clustering Coefficient",
                    "Node Degree")


  return(list(dissimilarity.results = network.dissimilarities,
              global.results = global.results,
              degree.results = degree.differences,
              perms = perms))

}
