#' Estimate a kendall based dependency network
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


estimate.kendall <- function(data, select="sig", alpha=0.05, threshold=0.95) {


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
  if(any(is.na(data))) print("full order pdcors estimated by removing missing data row wise")

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


    if(select == "sig") {

      pdc <- ppcor::pcor.test(x,y,z, method="kendall")

      edge.set$value[i] <- pdc$estimate
      edge.set$`P value`[i] <- pdc$p.value

    } else {

      pdc <- ppcor::pcor.test(x,y,z, method="kendall")

      edge.set$value[i] <- pdc$estimate
      edge.set$`P value`[i] <- NA

    }

  }

  #fill in saturated pdcor network matrix
  saturated.net[lower.tri(saturated.net)] <- edge.set$value
  saturated.net[upper.tri(saturated.net)] <- t(saturated.net)[upper.tri(saturated.net)]
  dimnames(saturated.net) <- list(colnames(data), colnames(data))


  #fill in P value matrix if "sig" used
  if(select == "sig") {

    #replace NA with 1 if no variability in binary node
    ind <- which(is.na(edge.set$`P value`) == TRUE)
    edge.set$`P value`[ind] <- 1

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

  #if select = qp return selected network as q1 graph structure
  if(select == "qp") {

    nrr <- matrix(0, ncol = nodes, nrow = nodes)

    # Loop through each pair of nodes (i, j)
    for (i in 1:(nodes - 1)) {

      v1 <- colnames(data)[i]

      for (j in (i + 1):nodes) {

        v2 <- colnames(data)[j]

        # Calculate partial associations conditioned on each other node k individually
        p.values <- rep(0, ncol(data))

        for (k in 1:nodes) {

          if (k != i && k != j) {

            #create dataset for (p-2)*2 q1 models and remove missing data rows at this stage
            d <- data[,c(i,j,k)]
            d <- na.omit(d)

            y <- d[,1]
            x <- d[,2]
            z <- d[,3]

            pdc.test <- ppcor::pcor.test(y,x,z, method="kendall")

            p.values[k] <- pdc.test$p.value

          }

        }

        p.values <- p.values[-c(i,j)]

        nrr[i,j] <- length(which(p.values > alpha)) / length(p.values)
        nrr[j,i] <- nrr[i,j]

      }

    }

    #graph structure based on nrr threshold
    ind <- which(nrr < threshold)
    selected.net[ind] <- 1

  }

  #if no sign - ignore sign section and return following
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



  #if no sign and qp select also return nrr
  if(select == "qp") {

    return(list(saturated.net = saturated.net,
                selected.net = selected.net * saturated.net,
                structure = selected.net,
                edge.set = edge.set,
                nrr = nrr,
                alpha = alpha,
                select = select,
                sign = sign))
  }


}
