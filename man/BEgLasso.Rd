% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BEgLasso.R
\name{BEgLasso}
\alias{BEgLasso}
\title{Batch Effect Correction using Group Graphical Lasso + glassoFast}
\usage{
BEgLasso(
  X0.glist,
  X1.glist,
  penal.rho,
  penal.ksi,
  penal.gamma,
  eps = 1e-04,
  print.detail = TRUE
)
}
\arguments{
\item{X0.glist}{A list of reference-batch matrices (each is samples × features).}

\item{X1.glist}{A list of matrices to correct (each is samples × features).}

\item{penal.rho}{Regularization parameter for glasso.}

\item{penal.ksi}{Regularization parameter for coefficient update.}

\item{penal.gamma}{Additional penalty parameter for coefficient update.}

\item{eps}{Convergence threshold (default 1e-4).}

\item{print.detail}{Logical flag to print iteration details.}
}
\value{
A list with elements:
  \item{Theta}{List of precision matrices for each group.}
  \item{X1.cor}{List of corrected matrices for group 1.}
  \item{coef.a}{Scaling coefficients.}
  \item{coef.b}{Offset coefficients.}
}
\description{
This function applies a batch effect correction method using the 
Group Graphical Lasso (GGM) approach, but replaces the inner coordinate-
descent with a fast Fortran-accelerated glasso implementation (glassoFast).
}
