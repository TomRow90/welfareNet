#' This contains various functions to compute different topological descriptors of networks.
#' They are called in the compare.networks function.
#'
#' Global network density
#'
#' @param net The estimated network structure
#' @return Network density estimate
#' @export

unweighted.density <- function(net) sum(abs(sign(net)[lower.tri(sign(net))])) / length(net[lower.tri(net)])

#' Global weighted density
#'
#' @param net The estimated network structure
#' @return Network weighted density estimate
#' @export

weighted.density <- function(net) sum(abs(net[lower.tri(net)]))

#' Proportion positive edges. Only relevant for signed networks.
#'
#' @param net The estimated network structure
#' @return Proportion of positive edges
#' @export

proportion.pos <- function(net) length(which(net > 0)) / length(net[lower.tri(net)])

#' Proportion negative edges. Only relevant for signed networks.
#'
#' @param net The estimated network structure
#' @return Proportion of negative edges
#' @export

proportion.neg <- function(net) length(which(net < 0)) / length(net[lower.tri(net)])

#' Average clustering coefficient
#'
#' @param net The estimated network structure
#' @export

global.cc <- function(net) NetworkToolbox::clustcoeff(net, weighted = FALSE)$CC

#' Average weighted clustering coefficient
#'
#' @param net The estimated network structure
#' @return Average weighted clustering coefficient
#' @export

global.weighted.cc <- function(net) mean(qgraph::clustcoef_auto(net)$clustZhang)

#' Node degree
#'
#' @param net The estimated network structure
#' @return Degree of each node. Only relevant in thresholded networks.
#' @export

degree <- function(net) NetworkToolbox::degree(net)

#' Node strength
#'
#' @param net The estimated network structure
#' @return Strength of each node. Only relevant in weighted networks.
#' @export

strength <- function(net) NetworkToolbox::strength(net)


#' Frobenius norm dissimilarity
#'
#' @param net1 The estimated network structure
#' @param net2 The estimated network structure for group 2
#' @param p The number of nodes
#' @return Estimated Frobenius norm dissimilarity. Function adapted from:
#' @export

sF <- function(net1, net2,p){
  # n1: weighted adjacency matrix for network 1
  # n2: weighted adjacency matrix for network 2
  # p: number of nodes
  sF <- 1 - (1/(1+(norm(net1 - net2, type="F")/sqrt(p/2))))

  return(sF)

}

#' Hamming distance
#'
#' @param net1 The estimated network structure
#' @param net2 The estimated network structure for group 2
#' @return Estimated Hamming dissimilarity. Only relevant for thresholded unweighted networks
#' @export

hamming.distance <- function(net1, net2) {

  nodes <- ncol(net1)
  edges <- (nodes * (nodes-1)) / 2

  net1 <- sign(net1)
  net1 <- net1[lower.tri(net1)]

  net2 <- sign(net2)
  net2 <- net2[lower.tri(net2)]

  d <- sum(net1 != net2) / edges

  return(d)

}

#' Jaccard dissimilarity
#'
#' @param net1 The estimated network structure for group 1
#' @param net2 The estimated network structure for group 2
#' @return Estimated Jaccard dissimilarity. Only relevant for thresholded unweighted networks
#' @export

jaccard.similarity <- function(net1, net2) {

  nodes <- ncol(net1)
  edges <- (nodes * (nodes-1)) / 2

  net1 <- sign(net1)
  net1 <- net1[lower.tri(net1)]

  net2 <- sign(net2)
  net2 <- net2[lower.tri(net2)]

  jaccard <- 1 - (sum(net1 & net2) / sum(net1 | net2))

  return(jaccard)

}



#' Extract perms
#'
#' @param perms the permuted data object
#' @param property the network statistic to extract
#' @return the permuted null distribution for the request network statistic
#' @export

extract.property <- function(test, property) {

  sapply(test, function(x) x[[property]])

}



#' One sided P value from permuted null distribution
#'
#' @param perms the permutation distribution
#' @param emp.diff the emperical estimate of the statistic
#' @return the permutatin test one sided P value
#' @export

one.sided.p.function <- function(perms, emp.diff) {

  ind <- which(perms >= emp.diff)

  p <- length(ind) / length(perms)

  return(p)

}

#' Two sided P value from permuted null distribution
#'
#' @param perms the permutation distribution
#' @param emp.diff the emperical estimate of the statistic
#' @return the permutatin test two sided P value
#' @export

two.sided.p.function <- function(perms, emp.diff) {

  ind <- which(abs(perms) >= abs(emp.diff))

  p <- length(ind) / length(perms)

  return(p)

}
