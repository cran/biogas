\name{stdVol}
\alias{stdVol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Correct Gas Volume to 'Standard' Conditions
}
\description{
  \code{stdVol} corrects gas volumes to dry conditions at a specified temperature and pressure.
}
\usage{
stdVol(vol, temp = 0, pres = 1, rh = 1, temp.std = NULL, pres.std = NULL, 
       unit.pres = "atm", unit.temp = "C")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vol}{
    measured gas volume at \code{temp} temperature, \code{pres} pressure, and \code{rh} relative humidity.
    A numeric vector.
}
  \item{temp}{
    temperature of gas in degrees C by default (see \code{unit.temp}).
    A numeric vector.
}
  \item{pres}{
    pressure of gas in atm by default (see \code{unit.pres}).
    A numeric vector.
}
  \item{rh}{
    relative humidity of the gas (\code{0<=rh<=1}).
    Typically should be 1.0 (the default).
    A numeric vector.
}
  \item{temp.std}{
    "standard" temperature.
    Default value is 0 degrees C.
}
  \item{pres.std}{
    "standard" pressure.
    Default value is 1.0 atm.
}
  \item{unit.pres}{
    pressure units for \code{pres} and \code{pres.std} arguments.
    Options are \code{"atm"} (the default), \code{"Pa"}, \code{"kPa"}, \code{"hPa"}, and \code{"bar"}.
}
  \item{unit.temp}{
    temperature units for \code{temp} and \code{temp.std} arguments.
    Options are \code{"C"} (degrees Celcius), \code{"F"} (degrees Fahrenheit), and \code{"K"} (Kelvin). 
}
}
\details{
  Standard temperature and pressure can be defined by the user.
  Default values are 0 degrees C and 1.0 atm.
  \code{stdVol} is vectorized.
}
\value{
  corrected gas volume in the same units as \code{vol}.
  A numeric vector.
}
\references{
  Hafner, S.D., Rennuit, C., Triolo, J.M., Richards, B.K. In review. A gravimetric method for measuring biogas production. \emph{Biomass and Bioenergy}.

  Richards, B.K., Cummings, R.J., White, T.E., Jewell, W.J. 1991. Methods for kinetic analysis of methane fermentation in high solids biomass digesters. \emph{Biomass and Bioenergy} 1: 65-73.
}
\author{
  Sasha D. Hafner and Charlotte Rennuit
}

\seealso{
  \code{\link{cumBg}},
  \code{\link{summBg}}
}
\examples{
  # 100 mL, measured at 35 C
  stdVol(100,temp=35)

  # And 1.02 atm
  stdVol(100,temp=35,pres=1.02)

  # Or, with different pressure units
  stdVol(100,temp=35,pres=103,unit.pres="kPa")

  # With different units and standard conditions
  stdVol(100,temp=100,pres=1,unit.pres="bar",pres.std=1,unit.temp="F",temp.std=68)

  # Vectorized
  data(vol)
  head(vol)
  vol$nvol<-stdVol(vol$vol,temp=20,pres=1.02)
  head(vol)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\concept{biogas}
