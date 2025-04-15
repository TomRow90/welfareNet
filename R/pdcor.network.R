#' Estimate a partial distance correlation network
#'
#' This function takes data and estimates a network based on partial distance correlation

#'
#' @param data data as matrix or data frame
#' @param select method for present/absent edge selection. Can either be "saturated" for no thresholding, "sig" for significance testing (default), of "qp" for q1-graph structure estimation
#' @param dcor.permutations if "sig" for select, how many permutations should be performed for the permutation test for each edge weight (default = 1000).
#' @param alpha if "sig" for select, the alpha level for thresholding the edges (default = 0.05)
#' @param adjust which method should be used for P value adjustment? default FDR using BH procedure
#' @param adjust.threshold what threshold should be used for the adjusted P value? default is 10% FDR
#' @param qp.threshold if qp graph is method for edge selection, what is the threshold of the non-rejection rate to select the graph structure (default = 0.95)? The alpha level argument is used for individual q tests
#'
#'
#' @return saturated.net = The saturated pdcor network
#' @return selected.net = The thresholded network if select = "sig" (or saturated network if select = "saturated")
#' @return p.vals = p value matrix if thresholding is selected
#' @return edge.set or edge.signs = edge list with estimates and P values
#' @return nrr = non-rejection rate matrix if "qp" is method of edge selection
#'
#' @export


pdcor.network <- function(data, select="sig", dcor.permutations=NULL, alpha=0.05, adjust = "none", adjust.threshold = NULL, qp.threshold = 0.95, zero.order=FALSE) {


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

  #check for missing data, remove and notify user
  if(any(is.na(data))) print("full order pdcors estimated by removing missing data row wise")

  #save d as new dataframe because if qp is selected then can use original data there without removing all rows with missing
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

      if(is.null(dcor.permutations)) dcor.permutations <- 1000

      pdc <- energy::pdcor.test(x,y,z, R = dcor.permutations)

      edge.set$value[i] <- pdc$estimate
      edge.set$`P value`[i] <- pdc$p.value

    } else {

      pdc <- energy::pdcor(x,y,z)

      edge.set$value[i] <- pdc
      edge.set$`P value`[i] <- NA

    }

  }

  #fill in saturated pdcor network matrix
  saturated.net[lower.tri(saturated.net)] <- edge.set$value
  saturated.net[upper.tri(saturated.net)] <- t(saturated.net)[upper.tri(saturated.net)]
  dimnames(saturated.net) <- list(colnames(data), colnames(data))


  if(select == "sig" & adjust != "none") {

    edge.set[,4] <- p.adjust(edge.set$`P value`, method=adjust)

  }


  if(select == "sig" & all(is.na(edge.set[,4]))) {

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

  } else {
    if(select == "sig" & !all(is.na(edge.set[,4]))) {

      #replace NA with 1 if no variability in binary node
      ind <- which(is.na(edge.set$`P value`) == TRUE)
      edge.set$`P value`[ind] <- 1

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

            pdc.test <- energy::pdcor.test(y,x,z, R=1000)

            p.values[k] <- pdc.test$p.value

          }

        }

        p.values <- p.values[-c(i,j)]

        nrr[i,j] <- length(which(p.values > alpha)) / length(p.values)
        nrr[j,i] <- nrr[i,j]

      }

    }

    #graph structure based on nrr threshold
    ind <- which(nrr < qp.threshold)
    selected.net[ind] <- 1

  }



    #if select = "saturated" return saturated matrix as selected network as well
    if(select == "saturated") {

      zero.order.selected.net <- zero.order.saturated.net
    }



  if(select == "sig") {

    return(list(saturated.net = saturated.net,
                selected.net = selected.net,
                p.values = net.p,
                structure = abs(sign(selected.net)),
                edge.set = edge.set,
                dcor.perms = dcor.permutations,
                alpha = alpha,
                select = select))
  }




  if(select == "qp") {

    return(list(saturated.net = saturated.net,
                selected.net = selected.net,
                structure = abs(sign(selected.net)),
                edge.set = edge.set,
                nrr = nrr,
                dcor.perms = dcor.permutations,
                alpha = alpha,
                select = select))
  }


}
