#' Estimate a q1 graph network using node wise generalized linear models
#'#'
#'
#' @param data data as matrix or data frame
#' @param alpha the alpha level of the maximum P value to select edges
#' @param nrr if this is equal to true then the non rejection rate procedure is used to select the network structure. Defaults to FALSE and the maximum P value method is used.
#' @param type character string of variable types e.g. "g" for gaussian, "b" for binary, "pr" for proportion/beta distributions, and "p" for Poisson
#' @param nrr.threshold if nrr is equal to TRUE this is the non rejection rate threshold for edge selection
#'
#' @return structure = the undirected unweighted network structure
#' @return nrr = non-rejection rate matrix if nrr is equal to TRUE
#' @return alpha = the alpha level used for the maximum P value
#'
#' @export


qp.mgm <- function(data, alpha=0.1, type, nrr=FALSE, nrr.threshold=NULL) {

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

  print(i)
  v <- unlist(strsplit(edge.set$edge[i], "-"))
  v1 <- v[1]
  v2 <- v[2]

  p.vals <- vector()

  for (j in 1:(nodes-2)) {

    print(j)
    ind <- which(colnames(data) == v1 | colnames(data) == v2)
    names <- colnames(data)[-ind]

    v3 <- names[j]

    ind.y <- which(colnames(data) == v1)
    ind.x <- which(colnames(data) == v2)

    node.type.y <- type[ind.y]
    node.type.x <- type[ind.x]

    #create dataset for (p-2)*2 q1 models and remove missing data rows at this stage
    d <- data[,c(v1,v2,v3)]
    d <- na.omit(d)

    y <- d[,1]
    x <- d[,2]
    z <- d[,3]


    if(node.type.y == "p" & node.type.x == "g") {

      mod <- lm(x ~ y + z)
      p.vals <- append(p.vals, summary(mod)$coefficients["y","Pr(>|t|)"])

    }

    if(node.type.y == "p" & node.type.x == "b") {

      mod <- glm(y ~ x + z, family="poisson")
      p.vals <-append(p.vals, summary(mod)$coefficients["x","Pr(>|z|)"])

    }

    if(node.type.y == "p" & node.type.x == "pr") {

      mod <- betareg::betareg(x ~ y + z, link="logit")
      p.vals <- append(p.vals, summary(mod)$coefficients$mean["y", "Pr(>|z|)"])


    }

    if(node.type.y == "g" & node.type.x == "p") {

      mod <- lm(y ~ x + z)
      p.vals <-append(p.vals, summary(mod)$coefficients["x","Pr(>|t|)"])

    }

    if(node.type.y == "pr" & node.type.x == "p") {

      mod <- betareg::betareg(y ~ x + z, link="logit")
      p.vals <- append(p.vals, summary(mod)$coefficients$mean["x", "Pr(>|z|)"])

    }

    if(node.type.y == "b" & node.type.x == "p") {

      mod <- glm(x ~ y + z, family="poisson")
      p.vals <-append(p.vals, summary(mod)$coefficients["y","Pr(>|z|)"])

    }

    if(node.type.y == "p" & node.type.x == "p") {

      mod <- glm(y ~ x + z, family="poisson")
      p.vals <-append(p.vals, summary(mod)$coefficients["x","Pr(>|z|)"])

    }


    if(node.type.y == "g" & node.type.x == "g") {

      mod <- lm(y ~ x + z)
      p.vals <- append(p.vals, summary(mod)$coefficients["x","Pr(>|t|)"])

    }

    if(node.type.y == "g" & node.type.x == "pr") {

      mod <- lm(y ~ x + z)
      p.vals <-append(p.vals, summary(mod)$coefficients["x","Pr(>|t|)"])

    }


    if(node.type.y == "g" & node.type.x == "b") {

      mod <- lm(y ~ x + z)
      p.vals <- append(p.vals, summary(mod)$coefficients["x","Pr(>|t|)"])

    }


    if(node.type.y == "b" & node.type.x == "g") {

      mod <- lm(x ~ y + z)
      p.vals <- append(p.vals, summary(mod)$coefficients["y","Pr(>|t|)"])

    }

    if(node.type.y == "pr" & node.type.x == "g") {

      mod <- lm(x ~ y + z)
      p.vals <- append(p.vals, summary(mod)$coefficients["y","Pr(>|t|)"])

    }

    if(node.type.y == "pr" & node.type.x == "pr") {

      mod <- betareg::betareg(y ~ x + z)
      p.vals <- append(p.vals, summary(mod)$coefficients$mean["x", "Pr(>|z|)"])

    }

    if(node.type.y == "pr" & node.type.x == "b") {

      mod <- betareg::betareg(y ~ x + z)
      p.vals <- append(p.vals, summary(mod)$coefficients$mean["x", "Pr(>|z|)"])
    }


    if(node.type.y == "b" & node.type.x == "pr") {

      mod <- betareg::betareg(x ~ y + z)
      p.vals <- append(p.vals, summary(mod)$coefficients$mean["y", "Pr(>|z|)"])

    }


    if(node.type.y == "b" & node.type.x == "b") {


      mod <- tryCatch({
        # Try fitting the model
        fit <- glm(y ~ x + z, family = "binomial")

      }

      , warning = function(w) {
        # Handle warnings gracefully (e.g., log them and proceed)
        message("Warning encountered: ", w$message)

        # Decide whether to return NA, or the fit model, depending on the warning
        fit <- glm(y ~ x + z, family = "binomial", method = "brglmFit", type = "AS_mixed")

      }, error = function(e) {
        # Handle errors by returning NA
        message("Error encountered: ", e$message)
        NA
      })

      p.vals <- append(p.vals, summary(mod)$coefficients["x","Pr(>|t|)"])


    }

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
