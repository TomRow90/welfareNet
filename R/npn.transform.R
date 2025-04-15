#' Transform non-gaussian continuous data to gaussian via the non-paranormal transformation
#'
#' This function takes the raw non-gaussian data and applies the non-paranormal transformation.
#' Missing data is maintained in the same entries as the raw data. Only supply the variables
#' you want to transform.
#'
#' @param data Path to the input file
#' @return a data frame of the same dataset after npn transformation
#' @export

npn.transform <- function(data) {

  for(i in 1:ncol(data)) {

    na.ind <- which(is.na(data[,i] == TRUE))

    if(length(na.ind) > 0) {

      transform <- huge::huge.npn(data[-na.ind,i])

      data[-na.ind,i] <- transform

    }

    if(length(na.ind) == 0) {

      transform <- huge::huge.npn(data[,i])

      data[,i] <- transform

    }

  }

  data <- as.data.frame(data)

  return(data=data)

}
