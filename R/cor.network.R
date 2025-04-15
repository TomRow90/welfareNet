#' Estimate a zero order correlation network (pearson, spearman, or kendall)
#'
#'
#' @param data data as matrix or data frame
#' @param select method for present/absent edge selection. Can either be "saturated" for no thresholding, "sig" for significance testing (default), of "qp" for q1-graph structure estimation
#' @param alpha if "sig" for select, the alpha level for thresholding the edges (default = 0.05)
#' @param threshold if qp graph procedure selected, threshold for NRR
#'
#' @return saturated.net = The saturated pdcor network
#' @return selected.net = The thresholded network if select = "sig" (or saturated network if select = "saturated")
#' @return p.vals = p value matrix if thresholding is selected
#' @return edge.set or edge.signs = edge list with estimates and P values
#'
#' @export


cor.network <- function(data, method = "spearman", use = "pairwise.complete.obs", select="sig", alpha=0.05) {


  #calculate node and edge sizes
  nodes <- ncol(data)
  edges <- nodes *(nodes - 1) / 2

  #create matrices/dataframes to store results in
  saturated.net <- matrix(0, nrow=ncol(data), ncol=ncol(data))
  dimnames(saturated.net) <- list(colnames(data), colnames(data))

  selected.net <- matrix(0, nrow=ncol(data), ncol=ncol(data))
  dimnames(selected.net) <- list(colnames(data), colnames(data))

  net.p <- matrix(0, nrow=ncol(data), ncol=ncol(data))
  dimnames(net.p) <- list(colnames(data), colnames(data))

  #edge list - saves calculating each edge weight twice in symmetric matrix form
  edge.set <- as.data.frame(matrix(0, nrow = edges, ncol = 3))
  colnames(edge.set) <- c("edge", "value","P value")
  pairs <- combn(colnames(data), 2, FUN = function(x) paste(x, collapse = "-"), simplify = FALSE)
  pairs <- as.character(pairs)
  edge.set$edge <- pairs

  #check for missing data, remove and notify user
  if(any(is.na(data))) print("correlation network estimated using pairwise complete observations")


  # estimate pcor and P value for each edge weight in edge list
  for (i in 1:nrow(edge.set)) {

    #get first node pair
    node.pair <- strsplit(edge.set[i,1], "-")
    v1 <- node.pair[[1]][1]
    v2 <- node.pair[[1]][2]

    v1.ind <- which(colnames(data) == v1)
    v2.ind <- which(colnames(data) == v2)

    x <- d[,v1.ind]
    y <- d[,v2.ind]

    obj <- cor(x,y, method = method, use = use)

    edge.set$value[i] <- obj$estimate
    edge.set$`P value`[i] <- obj$p.value

  }

  #fill in saturated cor network matrix
  saturated.net[lower.tri(saturated.net)] <- edge.set$value
  saturated.net[upper.tri(saturated.net)] <- t(saturated.net)[upper.tri(saturated.net)]
  dimnames(saturated.net) <- list(colnames(data), colnames(data))


  #fill in P value matrix if "sig" used
  if(select == "sig") {

    net.p <- matrix(0, nrow=nodes, ncol=nodes)

    net.p[lower.tri(net.p)] <- edge.set$`P value`

    net.p[upper.tri(net.p)] <- t(net.p)[upper.tri(net.p)]

    #select edges based on alpha threshold
    non.sig.ind <- which(net.p >= alpha)
    selected.net <- saturated.net
    selected.net[non.sig.ind] <- 0

  }

  #if select = "saturated" return saturated matrix as selected network as well
  if(select == "saturated") {

    selected.net <- saturated.net

  }


  #returns

  if(select == "sig") {

    return(list(saturated.net = saturated.net,
                selected.net = selected.net,
                p.values = net.p,
                structure = abs(sign(selected.net)),
                edge.set = edge.set,
                alpha = alpha,
                select = select))
  }

  if(select == "saturated") {

    return(list(saturated.net = saturated.net,
                selected.net = selected.net,
                edge.set = edge.set,
                select = select))
  }


}
