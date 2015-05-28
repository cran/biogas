\name{mass2vol}
\alias{mass2vol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate Biogas Production Gravimetrically
}
\description{
  \code{mass2vol} calculates biogas volume based on measured reactor mass loss.
}
\usage{
mass2vol(mass, xCH4, temp, pres = 1, unit.pres = "atm", unit.temp = "C", pres.std = NULL, 
         temp.std = NULL, value = "CH4", headspace = NULL, headcomp = "N2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mass}{
    reactor mass loss in g. A numeric vector.
}
  \item{xCH4}{
    mole fraction of methane in biogas (dry, methane + carbon dioxide only). A numeric vector.
}
  \item{temp}{
    the temperature of biogas as it exited the reactor. A length-one numeric vector.
    Degrees Celcius by default (see \code{unit.temp} argument).
}
  \item{pres}{
    the pressure of biogas as it exited the reactor. A length-one numeric vector.
    Atmospheres (atm) by default (see \code{unit.pres} argument).
}
  \item{unit.pres}{
    pressure units for \code{pres} and \code{pres.std} arguments.
    Options are \code{"atm"} (the default), \code{"Pa"}, \code{"kPa"}, \code{"hPa"}, and \code{"bar"}.
}
  \item{unit.temp}{
    temperature units for \code{temp} and \code{temp.std} arguments.
    Options are \code{"C"} (degrees Celcius, the default), \code{"F"}, and \code{"K"}.
}
  \item{pres.std}{
    standard pressure for presentation of biogas and methane results.
    Default value is 1.0 atm.
    Length one numeric vector.
}
  \item{temp.std}{
    standard temperature for presentation of biogas and methane results.
    Default value is 0 degrees C.
}
  \item{value}{
    what should be returned?
    Default is \code{"CH4"}, which returns methane volume in mL.
    \code{"all"} returns biogas, methane, and carbon dioxide volumes.
}
  \item{headspace}{
    (optional) reactor headspace volume in mL. 
    Used to correct for effect of initial reactor headspace on mass loss.
    A numeric vector.
}
  \item{headcomp}{
    (optional) composition of the initial reactor headspace.
    Only required if \code{headspace} argument is specified.
    Currently, the only option is \code{"N2"} for dinitrogen gas.
    Use of any other value will be ignored with a warning.
}
}
\details{
  This function uses the method described in Hafner et al. (in review) to calculate biogas production from reactor mass loss.
  It is essential that the only change in reactor mass is due to biogas removal!
  Users are advised to read the original reference before applying the method.
  This function is vectorized.
  Instead of using this function directly when working with multiple measurements on multiple reactors, use the \code{\link{cumBg}} function (which can call up \code{mass2vol}).
}
\value{
volume of methane in mL as a numeric vector (if \code{value="CH4"}) or a matrix with columns for biogas, methane, and carbon dioxide volumes in mL.
  \item{vBg}{standardized volume of biogas in mL}
  \item{vCH4}{standardized volume of methane in mL} 
  \item{vCH4}{standardized volume of carbon dioxyde in mL} 
  All volumes are standardized using \code{pres.std} and \code{temp.std} arguments (see \code{\link{stdVol}}).
}
\references{
  Hafner, S.D., Rennuit, C., Triolo, J.M., Richards, B.K. In review. A gravimetric method for measuring biogas production. Biomass and Bioenergy.
}
\author{
  Sasha D. Hafner and Charlotte Rennuit
}

\seealso{
  \code{\link{vol2mass}},
  \code{\link{stdVol}},
  \code{\link{cumBg}}
}
\examples{
  # Volume of methane if measured mass loss was 3.1 g
  mass2vol(3.1,xCH4=0.65,temp=35)
  # More details
  mass2vol(3.1,xCH4=0.65,temp=20,value="all")
  
  # Vectorized
  data("massw")
  massw$massloss <- massw$start - massw$end
  massw$vCH4 <- mass2vol(massw$massloss,xCH4=0.65,temp=37)
  massw
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\concept{biogas}
