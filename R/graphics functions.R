#' This contains various functions for graphical displays
#'
#' Change the node sizes given a node level descriptors
#'
#' @param data The data used in the network analyis (extracts columns names from this)
#' @param centrality A vector with a centrality score for each node
#' @param p1 value of parameter 1
#' @param p2 value of parameter 2
#'
#' @return vector of node sizes to be supplied to qgraph
#' @export


node.size <- function(data, centrality,p1,p2) {

  node.size <- rep(0, ncol(data))

  for (i in 1:ncol(data)) {
    if(centrality[i] <= 0) {
      node.size[i] <- (p1*exp(-ncol(data)/80) + p2*centrality[i])
    }
    else {
      node.size[i] <- (p1*exp(-ncol(data)/80) + p2*centrality[i])
    }
  }

  return(node.size)
}
