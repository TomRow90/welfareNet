#' Estimate a q1 graph network for network intervention analysis
#'
#' This function takes data which includes a binary intervention variable and estimates a q1-graph network
#'
#'
#' @param data data as matrix or data frame
#' @param alpha the alpha level for each individual test (default = 0.1)
#' @param nrr if this is equal to true then the non rejection rate procedure is used to select the network structure. Defaults to FALSE.
#' @param nrr.threshold if nrr is equal to TRUE this is the non rejection rate threshold for edge selection
#' @param num.boots number of bootstraps to calculate P value (default = 1000)
#' @param chores how many computer chores should be used when performing cdcor significance test
#'
#' @return structure = the undirected unweighted network structure. This is the network intervention analysis q1-graph network.
#' @return nrr = non-rejection rate matrix if nrr is equal to TRUE
#' @return alpha = the alpha level used in all significance tests
#'
#' @export


network.intervention.analysis <- function(data, alpha=0.1, nrr=FALSE, nrr.threshold=NULL, num.boots = 1000, chores=1) {


  nodes <- ncol(data)
  edges <- nodes * (nodes - 1) / 2

  if(nrr == TRUE) {

    nrr.m <- matrix(0, ncol = nodes, nrow = nodes)
    dimnames(nrr.m) <- list(colnames(data), colnames(data))

  }


  pairs <- unlist(pairs <- combn(colnames(data), 2, FUN = function(x) paste(x, collapse = "-"), simplify = FALSE))


  edge.set <- as.data.frame(matrix(0, nrow = edges, ncol = 2))
  colnames(edge.set) <- c("edge", "value")
  edge.set$edge <- pairs
  edge.set$value <- NA

  results <- as.data.frame(matrix(0, nrow = length(pairs) * (nodes-2), ncol=3))
  colnames(results) <- c("pair","test number","p value")
  results$pair <- rep(pairs, each = nodes -2)
  results$`test number` <- rep(1:(nodes-2), length(pairs))

  # Initialise progress bar
  pb <- txtProgressBar(min = 0, max = (nodes * (nodes - 1)) / 2, style = 3)
  progress <- 0  # Counter for progress tracking


  # Q1 graph procedure - Loop through each pair of nodes (i, j) to calculate their NRR
  for (i in 1:nrow(edge.set)) {

    v <- unlist(strsplit(edge.set$edge[i], "-"))
    v1 <- v[1]
    v2 <- v[2]

    p.vals <- vector()

    for (j in 1:(nodes-2)) {

      ind <- which(colnames(data) == v1 | colnames(data) == v2)
      names <- colnames(data)[-ind]

      v3 <- names[j]

      #create dataset for (p-2)*2 q1 models and remove missing data rows at this stage
      d <- data[,c(v1,v2,v3)]
      d <- na.omit(d)

      x <- d[,1]
      y <- d[,2]
      z <- d[,3]

      #conduct conditional distance independence test
      cdc.test <- cdcsis::cdcov.test(x,y,z, num.bootstrap=num.boots, num.threads = chores)

      #add to p value vector
      p.vals <- append(p.vals, cdc.test$p.value)

    }

    ind <- which(results$pair == edge.set$edge[i])

    results$`p value`[ind] <- p.vals

    # Update progress bar
    progress <- progress + 1
    setTxtProgressBar(pb, progress)

  }


  adj <- matrix(0,nrow=ncol(data), ncol=ncol(data))
  dimnames(adj) <- list(colnames(data), colnames(data))


  for(i in 1:nrow(edge.set)) {

    ind <- which(results$pair == edge.set$edge[i])

    res <- results[ind,]

    if(nrr == TRUE) {

      #calculate and store the NRR statistic in the corresponding element of the NRR matrix
      nrr.value <- length(which(res$`p value` > alpha)) / length(res$`p value`)
      edge.set$value[i] <- nrr.value

    } else {

      ifelse(any(res$`p value` > alpha), edge.set$value[i] <- 0, edge.set$value[i] <- 1)

    }

  }


  if(nrr == TRUE) {

    nrr.m[lower.tri(nrr.m)] <- edge.set$value
    nrr.m[upper.tri(nrr.m)] <- t(nrr.m)[upper.tri(nrr.m)]

    ind <- which(nrr.m < nrr.threshold)
    adj[ind] <- 1

  } else {

    adj[lower.tri(adj)] <- edge.set$value
    adj[upper.tri(adj)] <- t(adj)[upper.tri(adj)]

  }


  close(pb)

  if(nrr == TRUE) {

    return(list(structure = adj,
                nrr = nrr.m,
                alpha = alpha))
  } else {

    return(list(structure = adj,
                alpha = alpha))
  }


}
