% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dcor.aracne.R
\name{dcor.aracne}
\alias{dcor.aracne}
\title{Estimate conditional network structure with dcor based ARACNE algorithim}
\usage{
dcor.aracne(
  data,
  alpha = 0.05,
  adjust = "none",
  fdr = NULL,
  eps = 0.15,
  mod = "multiplicative",
  perms = 1000
)
}
\arguments{
\item{data}{data as matrix or data frame}

\item{alpha}{alpha level for selecting the initial zero order structure}

\item{adjust}{the P value adjustment for selecting the initial zero order structure. Set alpha to NULL if using an adjustment.}

\item{fdr}{the false discovery rate level edge selection will be controlled at}

\item{mod}{additive or multiplicative form of ARACNE algorithim}

\item{epsilon}{the error in the ARACNE algorithim}
}
\value{
structure = the selected network structure

dc.net = the zero order dcor network structure
}
\description{
Estimate conditional network structure with dcor based ARACNE algorithim
}
