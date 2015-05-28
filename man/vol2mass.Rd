\name{vol2mass}
\alias{vol2mass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate Mass of Biogas
}
\description{
  \code{vol2mass} calculates the mass of biogas removed, based on its composition, temperature, and pressure.
  This function is the inverse of \code{mass2vol}.
}
\usage{
vol2mass(volBg, xCH4, temp, pres = 1, unit.pres = "atm", unit.temp = "C")
}
\arguments{
  \item{volBg}{
  measured (not standardized) biogas volume in mL.
  Numeric vector
}
  \item{xCH4}{
  mole fraction of methane within biogas (dry, methane and carbon dioxide only).
  Numeric vector
}
  \item{temp}{
    temperature of biogas just prior to removal, in the units specified in \code{unit.temp} (default of degrees C).
    Numeric vector
}
  \item{pres}{
    pressure of biogas just prior to removal, in the units specified in \code{unit.pres} (default of atm).
    Numeric vector
}
  \item{unit.pres}{
    pressure units.
    Options are \code{"atm"} (the default), \code{"Pa"}, \code{"kPa"}, \code{"hPa"}, and \code{"bar"}.
    Length one character vector
}
  \item{unit.temp}{
    temperature units.
    Options are \code{"C"} (degrees Celcius, the default), \code{"F"}, and \code{"K"}.
    Length one character vector
}
}
\details{
  This function is vectorized. 
  Argument elements will be recycled as needed.
  Note that this function is conceptually but not numerically the inverse of \code{mass2vol}, because the \code{volBg} argument here is not standardized, and is assumed to saturated with water vapor just prior to removal.
}
\value{
  biogas mass in g as a numeric vector.
}
\references{
  Hafner, S.D., Rennuit, C., Triolo, J.M., Richards, B.K. In review. A gravimetric method for measuring biogas production. Biomass and Bioenergy.
}
\author{
  Sasha D. Hafner and Charlotte Rennuit
}

\seealso{
  \code{\link{mass2vol}}
}
\examples{
# mass of 100 mL biogas at 35 degrees C
vol2mass(100,xCH4=0.65,temp=35)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\concept{biogas}
