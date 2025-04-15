#' Estimate a network structure
#'
#' This function is a wrapper around the various network estimation procedures implemented in this package
#' and other network analysis R packages. See details for estimation algorithms implemented.
#'
#' @param estimator = the network estimation algorithm to use
#' @param ... = the arguments specific for the estimation function selected
#'
#' @return structure = the selected network structure
#'
#' @export


estimate.network <- function(estimator, ...) {

  estimator.map <- list(
    dcor = welfareNet::dcor.network,
    cor = welfareNet::cor.network,
    qp.mgm = welfareNet::qp.mgm,
    dcor.aracne = welfareNet::dcor.aracne,
    qp.cdit = welfareNet::qp.cdit.network,
    qp.mgm.nrr = welfareNet::qp.mgm,
    qp.cdit.nrr = welfareNet::qp.cdit,
    pcor.pear = GGMnonreg::ggm_inference,
    pcor.spear = GGMnonreg::ggm_inference,
    pcor.kend = GGMnonreg::ggm_inference,
    npn.GGM = GGMnonreg::ggm_inference,
    npn.regGGM = corpcor::pcor.shrink,
    mgm.OR = mgm::mgm,
    mgm.AND = mgm::mgm,
    npn.logo = NetworkToolbox::LoGo,
    npn.tmfg = NetworkToolbox::TMFG
  )

  if (!estimator %in% names(estimator.map)) {
    stop(paste("Unknown estimator:", estimator))
  }

  if(estimator == "npn.GGM" | estimator == "npn.regGGM" | estimator == "npn.logo" | estimator == "npn.tmfg") {

    data <- npn.transform(data)

  }


  if(estimator == "mgm.OR" | estimator == "mgm.AND") {

    data <- as.matrix(data)

  }

  do.call(estimator.map[[estimator]], data = data, list(...))

}
