#' Estimate a dependency network with distance precision network
#'
#' This function takes data and estimates the distance precision network. Missing data is removed.
#'
#' @param data data as matrix or data frame
#' @param select method for present/absent edge selection. Can either be "saturated" for no thresholding, "k-means" for clustering, of "qp" for q1-graph structure estimation
#' @param threshold if qp graph is method for edge selection, what is the threshold of the non-rejection rate to select the graph structure (default = 0.95)? The alpha level argument is used for individual q tests
#'
#' @return saturated.net = The saturated pdcor network
#' @return selected.net = The thresholded network (or saturated network if select = "saturated")
#' @return edge.set = edge list with estimates
#'
#' @export

estimate.dpm <- function(data, select="k-means") {


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
  edge.set <- as.data.frame(matrix(0, nrow = edges, ncol=2))
  colnames(edge.set) <- c("edge", "value")
  pairs <- combn(colnames(data), 2, FUN = function(x) paste(x, collapse = "-"), simplify = FALSE)
  pairs <- as.character(pairs)
  edge.set$edge <- pairs

  #check for missing data, remove and notify user
  if(any(is.na(data))) print("full order pdcors estimated by removing missing data row wise")

  data <- na.omit(data)

  dpm <- DPM::reg.dpm(data)

  #fill in saturated pdcor network matrix
  saturated.net[lower.tri(saturated.net)] <- dpm[lower.tri(dpm)]
  saturated.net[upper.tri(saturated.net)] <- t(saturated.net)[upper.tri(saturated.net)]
  dimnames(saturated.net) <- list(colnames(data), colnames(data))


  #fill in P value matrix if "sig" used
  if(select == "k-means") {

    k.means <- DPM::kmeans_links(saturated.net)
    threshold <- k.means$threshold

    ind <- which(saturated.net < threshold)

    selected.net <- saturated.net

    selected.net[ind] <-0

  }

  #if select = "saturated" return saturated matrix as selected network as well
  if(select == "saturated") {

    selected.net <- saturated.net
  }


    return(list(saturated.net = saturated.net,
                selected.net = selected.net,
                structure = abs(sign(selected.net)),
                select = select))


}
