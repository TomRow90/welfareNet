#' Estimate conditional network structure with dcor based ARACNE algorithim
#'
#'
#' @param data data as matrix or data frame
#' @param alpha alpha level for selecting the initial zero order structure
#' @param adjust the P value adjustment for selecting the initial zero order structure. Set alpha to NULL if using an adjustment.
#' @param fdr the false discovery rate level edge selection will be controlled at
#' @param epsilon the error in the ARACNE algorithim
#' @param mod additive or multiplicative form of ARACNE algorithim
#'
#' @return structure = the selected network structure
#' @return dc.net = the zero order dcor network structure
#'
#' @export

aracne <- function(data, alpha = 0.05, adjust="none", fdr=NULL, eps=0.15, mod="multiplicative", perms=1000) {

  # Define the helper function outside the loop
  test_edge <- function(sub.net, i, j, k, eps, mod) {
    if (mod == "additive") {
      if (sub.net[i, j] <= min(sub.net[i, k] - eps, sub.net[j, k] - eps)) {
        sub.net[i, j] <- 0
        sub.net[j, i] <- 0
      }
    } else if (mod == "multiplicative") {
      if (sub.net[i, j] <= min(sub.net[i, k] * (1 - eps), sub.net[j, k] * (1 - eps))) {
        sub.net[i, j] <- 0
        sub.net[j, i] <- 0
      }
    }
    return(sub.net)
  }

  nodes <- ncol(data)
  edges <- nodes * (nodes - 1) / 2

  structure <- matrix(NA, nrow=nodes, ncol=nodes)
  dimnames(structure) <- list(colnames(data), colnames(data))

  #estimate marginal dcor network
  dc.net <- welfareNet::estimate.dcor.network(data, select="sig", dcor.permutations=perms, alpha=alpha, adjust=adjust, adjust.threshold=fdr)
  net <- dc.net$selected.net

  #extract all connected triplets
  g <- igraph::graph_from_adjacency_matrix(net, mode="undirected", weighted=TRUE, diag=FALSE)
  trips <- igraph::cliques(g, min = 3, max = 3)


  if(length(trips) > 0) {


    trip.m <- as.data.frame(matrix(names(unlist(trips)), nrow=length(trips), ncol=3, byrow=TRUE))
    colnames(trip.m) <- c("v1","v2","v3")

    #for each triplet perform the DPI pruning steps
    for(i in 1:nrow(trip.m)) {

      v1 <- trip.m$v1[i]
      v2 <- trip.m$v2[i]
      v3 <- trip.m$v3[i]

      sub.net <- net[c(v1,v2,v3), c(v1,v2,v3)]

      # Apply the helper function to all three edges in the triangle
      sub.net <- test_edge(sub.net, 1, 2, 3, eps, mod)
      sub.net <- test_edge(sub.net, 1, 3, 2, eps, mod)
      sub.net <- test_edge(sub.net, 2, 3, 1, eps, mod)

      ind <- which(structure[c(v1,v2,v3), c(v1,v2,v3)] == 0)
      sub.net[ind] <- 0

      # Update the global structure with the pruned submatrix for these nodes
      structure[c(v1,v2,v3), c(v1,v2,v3)] <- sub.net

    }

    ind <- which(is.na(structure) == TRUE)
    structure[ind] <- 0

    #return the final pruned network
    return(list(structure = structure,
                dcor.net = dc.net))

  } else {

    if(length(trips) == 0) {

      return(list(structure = net,
                  dcor.net = dc.net))
    }

  }

}
