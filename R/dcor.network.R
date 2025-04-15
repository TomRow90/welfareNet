#' Estimate a zero order dependency network from distance correlations
#'
#' This function takes data and estimates the zero order distance correlation network
#' It returns the saturated dcor network (used if network moderation analysis performed)
#' as well as the 'selected' network by thresholding at FDR < 0.1 (user can change this setting).
#' Missing data is handled using pairwise complete observations.
#'
#' @param data data as matrix or data frame
#' @param select method for present/absent edge selection. Can either be "saturated" for no thresholding, "sig" for significance testing (default), of "qp" for q1-graph structure estimation
#' @param dcor.permutations if "sig" for select, how many permutations should be performed for the permutation test for each edge weight (default = 1000).
#' @param alpha if "sig" for select, the alpha level for thresholding the edges (default = 0.05)
#' @param adjust which method should be used for P value adjustment? default FDR using BH procedure
#' @param adjust.threshold what threshold should be used for the adjusted P value? default is 10% FDR
#'
#' @return saturated.net = The saturated dcor network
#' @return selected.net = The thresholded network if select = "sig" (or saturated network if select = "saturated")
#' @return p.vals = p value matrix if thresholding is selected
#' @return edge.set = edge list with estimates and P values
#'
#' @export


dcor.network <- function(data, select="sig", dcor.permutations=NULL, alpha=0.05, adjust = "fdr", adjust.threshold = 0.1) {

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
  edge.set <- as.data.frame(matrix(0, nrow = edges, ncol = 4))
  colnames(edge.set) <- c("edge", "value","P value", adjust)
  pairs <- combn(colnames(data), 2, FUN = function(x) paste(x, collapse = "-"), simplify = FALSE)
  pairs <- as.character(pairs)
  edge.set$edge <- pairs
  edge.set[,4] <- NA

  # estimate dcor and P value for each edge weight in edge list
  for (i in 1:nrow(edge.set)) {

    #get first node pair
    node.pair <- strsplit(edge.set[i,1], "-")
    v1 <- node.pair[[1]][1]
    v2 <- node.pair[[1]][2]

    v1.ind <- which(colnames(data) == v1)
    v2.ind <- which(colnames(data) == v2)

    #data for node pair and high dimensional control variable
    x <- data[,v1.ind]
    y <- data[,v2.ind]

    d <- na.omit(cbind(x,y))
    colnames(d) <- c("x","y")

    x <- d[,1]
    y <- d[,2]

    if(select == "sig") {

      if(is.null(dcor.permutations)) dcor.permutations <- 1000

      dc <- energy::dcor.test(x,y, R = dcor.permutations)

      edge.set$value[i] <- dc$statistic
      edge.set$`P value`[i] <- dc$p.value

    } else {

      dc <- energy::dcor(x,y)

      edge.set$value[i] <- dc
      edge.set$`P value`[i] <- NA

    }

  }

  #fill in saturated dcor network matrix
  saturated.net[lower.tri(saturated.net)] <- edge.set$value
  saturated.net[upper.tri(saturated.net)] <- t(saturated.net)[upper.tri(saturated.net)]
  dimnames(saturated.net) <- list(colnames(data), colnames(data))


  if(select == "sig" & adjust != "none") {

    edge.set[,4] <- p.adjust(edge.set$`P value`, method=adjust)

  }


  if(select == "sig" & all(is.na(edge.set[,4]))) {

    #matrix to store P values
    net.p <- matrix(0, nrow=nodes, ncol=nodes)

    #fill in
    net.p[lower.tri(net.p)] <- edge.set$`P value`
    net.p[upper.tri(net.p)] <- t(net.p)[upper.tri(net.p)]

    #select edges based on alpha threshold
    non.sig.ind <- which(net.p >= alpha)
    selected.net <- saturated.net
    selected.net[non.sig.ind] <- 0

  } else {

    if(select == "sig" & !all(is.na(edge.set[,4]))) {

      net.p <- matrix(0, nrow=nodes, ncol=nodes)

      net.p[lower.tri(net.p)] <- edge.set$`P value`

      net.p[upper.tri(net.p)] <- t(net.p)[upper.tri(net.p)]

      adjust.p <- p.adjust(net.p[lower.tri(net.p)], method=adjust)

      adjust.values <- matrix(0, nrow=nodes, ncol=nodes)
      adjust.values[lower.tri(adjust.values)] <- adjust.p
      adjust.values[upper.tri(adjust.values)] <- t(adjust.values)[upper.tri(adjust.values)]

      #select edges based on adjust threshold
      non.sig.ind <- which(adjust.values >= adjust.threshold)
      selected.net <- saturated.net
      selected.net[non.sig.ind] <- 0

    }
  }



  #if select = "saturated" return saturated matrix as selected network as well
  if(select == "saturated") {

    selected.net <- saturated.net
  }


  return(list(saturated.net = saturated.net,
              selected.net = selected.net,
              p.values = net.p,
              structure = abs(sign(selected.net)),
              edge.set = edge.set,
              dcor.perms = dcor.permutations,
              alpha = alpha,
              select = select))

}


