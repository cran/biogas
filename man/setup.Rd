\name{setup}
\alias{setup}
\docType{data}
\title{
  Setup Details for Batch Reactors
}
\description{
  Information on reactor substrates and substrate and inoculum masses for nine anaerobic reactors.
}
\usage{data("setup")}
\format{
  A data frame with 9 observations on the following 4 variables:
  \describe{
    \item{\code{id}}{identification code, a unique value for each reactor in the dataset. A factor with levels \code{2_1}, \code{2_2}, \code{2_3}, \code{2_4}, \code{2_5}, \code{2_6}, \code{2_10}, \code{2_11}, and \code{2_12} }
    \item{\code{descrip}}{description of substrate, a factor with levels \code{A} \code{B} \code{inoc}.}
    \item{\code{msub}}{mass of substrate volatile solids (VS) added to the reactor at the start in g, a numeric vector}
    \item{\code{minoc}}{mass of inoculum added to the reactor at the start in g, a numeric vector}
  }
}
\details{
  These data are meant to be example data for \code{\link{summBg}}.
  Reactors were 500 mL glass serum bottles with butyl rubber septa and screw caps. 
  Data in \code{vol}, \code{mass}, \code{massw}, \code{comp}, and \code{setup} are from the same reactors.
}
\source{
  Measurements made by Charlotte Rennuit and Ali Heidarzadeh Vazifehkhoran
}
\examples{
data(setup)
}
\keyword{datasets}
\concept{biogas}
