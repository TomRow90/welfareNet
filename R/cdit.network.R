#' Estimate a dependency network using conditional distance independence tests with multivariate control
#'
#' This function takes data and estimates the dependence structure using conditional
#' distance independence tests. Missing data is removed.
#'
#' @param data data as matrix or data frame
#' @param boots number of bootstraps for local bootstrapping P value computation

#' @return adjacency matrix
#'
#' @export



cdit.network <- function(data, boots, chores) {

  nodes <- ncol(data)
  edges <- nodes * (nodes - 1) / 2

  #p value matrix
  p.values <- matrix(NA, nrow=nodes, ncol=nodes)

  #edge list - saves calculating each edge weight twice in symmetric matrix form
  edge.set <- as.data.frame(matrix(0, nrow = edges, ncol = 3))
  colnames(edge.set) <- c("edge", "value","P value")
  pairs <- combn(colnames(data), 2, FUN = function(x) paste(x, collapse = "-"), simplify = FALSE)
  pairs <- as.character(pairs)
  edge.set$edge <- pairs

  #save d as new dataframe because if qp is selected then can use original data there without removing all rows
  d <- na.omit(data)

  # estimate pdcor and P value for each edge weight in edge list
  for (i in 1:nrow(edge.set)) {

    #get first node pair
    node.pair <- strsplit(edge.set[i,1], "-")
    v1 <- node.pair[[1]][1]
    v2 <- node.pair[[1]][2]

    v1.ind <- which(colnames(data) == v1)
    v2.ind <- which(colnames(data) == v2)

    #data for node pair and high dimensional control variable
    x <- d[,v1.ind]
    y <- d[,v2.ind]
    z <- d[,-c(v1.ind,v2.ind)]


    cdc <- cdcsis::cdcov.test(x,y,z, num.bootstrap=boots, num.threads = chores)

    edge.set$`P value`[i] <- cdc$p.value


  }

  #fill in saturated pdcor network matrix
  selected.net <- matrix(0, nrow=nodes, ncol=nodes)
  ind <- which(edge.set$`P value` < 0.05)
  selected.net[lower.tri(selected.net)][ind] <- 1
  selected.net[upper.tri(selected.net)] <- t(selected.net)[upper.tri(selected.net)]
  dimnames(selected.net) <- list(colnames(data), colnames(data))

  #p vals
  p.values[lower.tri(p.values)] <- edge.set$`P value`
  p.values[upper.tri(p.values)] <- t(p.values)[upper.tri(p.values)]

  return(list(selected.net = selected.net,
              p.values = p.values))

}


