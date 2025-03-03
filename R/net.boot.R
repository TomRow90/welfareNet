#' Perform bootstrap resampling of network estimate
#'
#' This function performs bootstrapping analysis (using the boot package) to assess the stability and accuracy of edge weights.
#'
#' @param data the dataset used to estimate the emperical network (as matrix or data frame)
#' @param net the estimate.network object which contains the saturated and selected network estimates
#' @param nBoots the number of bootstrap resamples to perform. Default is 1000.
#' @param stability Should a network structure be selected on each resample such that edge weight stability can be calculated?
#'
#' @return stability results = a dataframe with parameters for each edge weight
#' @return boot.object = the object from the boot package
#'
#' @export


net.boot <- function(data, net, nBoots, stability = TRUE) {

  pb <- txtProgressBar(min = 0, max = nBoots, style = 3)

  nodes <- ncol(data)
  edges <- nodes * (nodes - 1) / 2

  #initialise some arguements
  if(stability == TRUE) {

    select <- net$select
    dcor.permutations <- net$dcor.perms
    alpha <- net$alpha
    sign <- net$sign

  } else {

    select <- "saturated"
    dcor.permutations <- NULL
    alpha <- NULL
    sign <- net$sign

  }


  #define estimated network
  saturated.net <- net$saturated.net
  selected.net <- net$selected.net


  boot.fun <- function(data, indices, select, dcor.permutations, alpha, sign) {

    #pb bar
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)

    # Subset the data based on the bootstrap indices
    data <- data[indices,]

    # Fit the model using qp.net and the type argument
    mod <- estimate.network(data, select=select, dcor.permutations = dcor.permutations, alpha = alpha, sign = sign)

    #return saturated estimates
    sat <- mod$saturated.net[lower.tri(mod$saturated.net)]

    #return selected estimates
    sel <- mod$selected.net[lower.tri(mod$selected.net)]

    # Return the lower triangular part of the adjacency matrix
    return(c(sat,sel))

  }

  boot.obj <- boot::boot(data, statistic=boot.fun, R=nBoots, select=select, dcor.permutations = dcor.permutations, alpha = alpha, sign = sign)

  #create summary dataframe from obtained boots
  stability.results <- as.data.frame(matrix(0, ncol=6, nrow=edges))
  colnames(stability.results) <- c("Edge", "Estimate", "Median", "Lower bound", "Upper bound", "Stability")
  pairs <- combn(colnames(data), 2, FUN = function(x) paste(x, collapse = "-"), simplify = FALSE)
  pairs <- as.character(pairs)
  stability.results$Edge <- pairs
  stability.results$Estimate <- saturated.net[lower.tri(saturated.net)]
  stability.results$Stability <- NA

  for(i in 1:edges) {

    stability.results$Median[i] <- median(boot.obj$t[,i], na.rm = TRUE)

    stability.results$`Lower bound`[i] <- boot::boot.ci(boot.obj, index=i, type="perc")$percent[4]

    stability.results$`Upper bound`[i] <- boot::boot.ci(boot.obj, index=i, type="perc")$percent[5]

  }

  if(stability == "TRUE") {

    for(i in 1:edges) {

      stability.results$Stability[i] <- length(which(sign(boot.obj$t[,i+edges]) == sign(stability.results$Estimate[i]))) / nBoots

    }
  }

  close(pb)

  return(list(stability.results = stability.results,
              boot.object = boot.obj))

}
