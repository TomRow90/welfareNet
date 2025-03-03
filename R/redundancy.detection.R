#' Identify possibly redundant nodes based on zero-order correlation
#'
#' This function estimates pairwise zero-order distance correlations.
#' It returns the correlation matrix, a matrix of all pairwise results, and a matrix of possibly redundant
#' variables/nodes based on a user supplied correlation threshold. Users can also supply a character vector
#' containing the names of variables/nodes which should not be suggested to be redundant with any other node.
#'
#' @param data the data to estimate correlations from
#' @return redundancy.results = dataframe with all dcors and suggested redundancies
#' @return sig.redundancy.results = dataframe with only pairs over threshold shown
#' @return dcor.net = pairwise zero-order distance correlation network/matrix
#' @return threshold = the threshold used to decide about possible redundancies
#' @export



redundancy.detection <- function(data, retain.nodes=NULL,retain.edges=NULL, threshold=0.4, permutations=1000) {

  nodes <- ncol(data)
  edges <- (nodes * (nodes - 1)) / 2
  node.names <- colnames(data)

  edge.set <- as.data.frame(matrix(0, nrow = edges, ncol = 4))
  colnames(edge.set) <- c("edge", "dcor","p.value", "action")
  pairs <- combn(colnames(data), 2, FUN = function(x) paste(x, collapse = "-"), simplify = FALSE)
  pairs <- as.character(pairs)
  edge.set$edge <- pairs

  for (i in 1:nrow(edge.set)) {

    #identify individual nodes in pair
    node.pair <- strsplit(edge.set[i,1], "-")
    v1 <- node.pair[[1]][1]
    v2 <- node.pair[[1]][2]

    v1.ind <- which(colnames(data) == v1)
    v2.ind <- which(colnames(data) == v2)

    d <- data[,c(v1,v2)]
    d <- na.omit(d)

    dc.test <- energy::dcor.test(d[,v1], d[,v2], R = permutations)

    edge.set$dcor[i] <- dc.test$statistic
    edge.set$p.value[i] <- dc.test$p.value

  }

  #create pairwise dcor network
  saturated.net <- matrix(0, nrow=nodes, ncol=nodes)
  saturated.net[lower.tri(saturated.net)] <- edge.set$dcor
  saturated.net[upper.tri(saturated.net)] <- t(saturated.net)[upper.tri(saturated.net)]
  dimnames(saturated.net) <- list(colnames(data), colnames(data))

  #create p.value matrix
  p.vals <- matrix(0, nrow=nodes, ncol=nodes)
  p.vals[lower.tri(p.vals)] <- edge.set$p.value
  p.vals[upper.tri(p.vals)] <- t(p.vals)[upper.tri(p.vals)]

  #finds "significant" values from results object based on threshold
  ind <- which(edge.set$dcor >= threshold)
  sig.red <- edge.set[ind,]
  sig.red$action <- "choose"

  #create structure of zero order network
  zero.order.adj <- matrix(0, nrow=nodes, ncol=nodes)
  ind1 <- which(p.vals < 0.05)
  zero.order.adj[ind1] <- 1

  #weighted zero order network
  zero.order.wadj <- zero.order.adj * saturated.net

  #select redundancies from list

  for(i in 1:nrow(sig.red)) {

    split.string <- strsplit(sig.red$edge[i], "-")
    node1 <- split.string[[1]][1]
    node2 <- split.string[[1]][2]

    if(!is.null(retain.edges) && length(retain.edges) > 0 && sig.red$pair[i] %in% retain.edges) {

      sig.red$action[i] <- "no.action"

    } else if(node1 %in% retain.nodes & node2 %in% retain.nodes) {

      sig.red$action[i] <- "no.action"

    } else {

      #overrides
      if(node1 %in% retain.nodes & !node2 %in% retain.nodes) sig.red$action[i] <- node2
      if(node2 %in% retain.nodes & !node1 %in% retain.nodes) sig.red$action[i] <- node1

    }
  }


  return(list(redundancy.results=edge.set,
              sig.redundancy.results=sig.red[order(sig.red$dcor, decreasing=TRUE), ],
              dcor.net = saturated.net,
              zero.order.structure = zero.order.adj,
              zero.order.wadj = zero.order.wadj,
              threshold = threshold))

}
