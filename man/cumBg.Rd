\name{cumBg}
\alias{cumBg}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
  Calculate Cumulative Biogas Production 
}
\description{
  \code{cumBg} (for \emph{cum}ulative \emph{b}io\emph{g}as) calculates cumulative biogas, methane production and production rates from individual volume and composition measurements for any number of reactors.
}
\usage{
cumBg(dat, dat.type="vol", comp="NULL", temp="NULL", pres = 1, id.name="id", 
      time.name="time",dat.name=dat.type, comp.name="xCH4", headspace=NULL, 
      vol.hs.name="vol.hs", temp.std = NULL, pres.std = NULL, unit.pres = "atm", 
      unit.temp = "C", cmethod = "removed", imethod = "linear", extrap = FALSE, 
      addt0= TRUE, showt0 = if(any(dat[,time.name]==0)) TRUE else FALSE) 
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{dat}{
    a data frame with the columns reactor identification code; time of measurement (as \code{numeric}, or \code{POSIX}); and measured biogas volume or total reactor mass, in g (see \code{dat.type} argument).
    The names of these columns can be specified with \code{id.name}, \code{time.name}, and \code{dat.name}.
    Or, the default names can be used and if the column are not identified matching can be based on order (\code{id.name}, \code{time.name}, followed by \code{dat.name}).
    Additional columns can be present--these will be returned in the output data frame.
}
  \item{dat.type}{
  the type of data contained in \code{dat}. 
  Value can be \code{"vol"} if the data is measured biogas volume or \code{"mass"} if the reactor masses are used. 
  Default is \code{"vol"}.
}
  \item{comp}{
    (optional) a data frame with the columns reactor identification code; time of measurement, (as \code{numeric}, or \code{POSIX}); and methane concentration within dry biogas as a mole fraction, considering only methane and carbon dioxide (unless \code{cmethod="total"}).
    If omitted, cumulative biogas volume will still be calculated and returned (but no methane data will be returned).
    The names of these columns can be specified with \code{id.name}, \code{time.name}, and \code{comp.name}.
    Or, the default names can be used and if the column are not identified matching can be based on order (\code{id.name}, \code{time.name}, followed by \code{comp.name}).
}
  \item{temp}{
    the temperature at which biogas volume was measured (when \code{dat.type="vol"}), or of biogas as it exited the reactor (when \code{dat.type="mass"}).
    A length-one numeric vector.
    Degrees Celcius by default (see \code{unit.temp} argument). 
    Default is \code{NULL}, which suppresses correction for temperature and pressure. 
}
  \item{pres}{
    the pressure at which biogas volume was measured (when \code{dat.type="vol"}), or of biogas as it exited the reactor (when \code{dat.type="mass"}).
    A length-one numeric vector.
    Presure is set to 1 atmosphere (atm) by default (see \code{unit.pres} argument). 
}
  \item{id.name}{
    name of the reactor identification code column in \code{dat}. Should be the same in all data frames used in the function. Default is \code{"id"}.
}
  \item{time.name}{
    name of column containing time data in \code{dat} and \code{comp} data frames. Default is \code{"time"}.
}
  \item{dat.name}{
    name of column containing the data in \code{dat} data frame. Default is \code{dat.type} value. See \code{dat.type} argument.
}
  \item{comp.name}{
    name of column containing gas composition data in \code{comp} data frame. Default is \code{"xCH4"}.
}
  \item{headspace}{
    an optional data frame or numeric vector with reactor headspace volumes.
    If a data frame is used, it should at least contain a \code{"id"} (reactor identification code) column (see \code{"id.name"})  and headspace volume column (see \code{vol.hs.name} argument).
    If both of these columns namesare not correctly specified, matching will be done based on order.
}
  \item{vol.hs.name}{
    name of column containing headspace volume data in optional \code{headspace} dataframe. Default is \code{"vol.hs"}.
}
  \item{temp.std}{
    standard temperature for presentation of biogas and methane results.
    Default value is 0 degrees C.
    See \code{\link{stdVol}}.
    Length one numeric vector.
}
  \item{pres.std}{
    standard pressure for presentation of biogas and methane results.
    Default value is 1.0 atm.
    See \code{\link{stdVol}}.
    Length one numeric vector.
}
  \item{unit.pres}{
    pressure units for \code{pres} and \code{pres.std} arguments.
    Options are \code{"atm"}, \code{"Pa"}, \code{"kPa"}, \code{"hPa"}, and \code{"bar"}. Default is \code{"atm"}.
}
  \item{unit.temp}{
    temperature units for \code{temp} and \code{temp.std} arguments.
    Options are \code{"C"} (degrees Celcius), \code{"F"} (degrees Fahrenheit), and \code{"K"} (Kelvin). 
    Default is \code{"C"}. 
    See \code{\link{stdVol}}.
}
  \item{cmethod}{
    the method for calculating cumulative methane production.
    Use \code{"removed"} to base production on \code{xCH4} and gas volumes removed.
    Use \code{"total"} to base it on the sum of methane removed and methane remaining in the reactor headspace.
    \emph{In this case, \code{xCH4} should be calculated including all biogas components except water.}
    Length one character vector. Default is \code{"removed"}.
}
  \item{imethod}{
    method used for interpolation of \code{xCH4}.
    This is passed as the \code{method} argument to \code{\link{interp}}.
    Length one character vector. Default is \code{linear} for linear interpolation.
}
  \item{extrap}{
    should \code{comp.name} be extrapolated? 
    Set to \code{TRUE} if extrapolation wanted. 
    Length one logical vector.
    This is passed as the \code{extrap} argument to \code{\link{interp}}. 
    Default is \code{FALSE}.
}
 \item{addt0}{
    is the earliest time in \code{dat} data frame ``time zero'' (start time)? 
    If not, this argument adds a row with \code{time.name=0} for each reactor in order to calculate production rates for the first observation.
    This addition is only made when \code{dat.type="vol"} and \code{time.name} is \code{numeric} (or \code{integer}).
    Length-one logical vector.
    Default is \code{TRUE}.
    To return these additional rows in the output, see \code{showt0}.
}
  \item{showt0}{
    should ``time zero'' rows be returned in the output?
    Only applies if \code{time.name} is \code{numeric} (or \code{integer}).
    These rows may have been present in the original data (\code{dat}) or added by the function (see   \code{addt0}).
    Default value depends on \code{dat} \code{time.name} column content. 
    If \code{time.name} column is numeric and contains 0 then the default value is \code{TRUE} and otherwise \code{FALSE}. 
}
}
\details{
  Using volume or mass data from \code{dat} and, optionally, gas composition from \code{comp}, this function will calculate normalized biogas and methane production (if \code{comp} is provided) for each observation,  interpolating \code{comp.name} (from \code{comp} argument) to each \code{time.name} in \code{dat} if needed, and summing these for cumulative values.
  Use of volumetric measurements is the default case. 
  Alternatively, if reactor mass is included in \code{dat}, a gravimetric approach (Hafner et al. in review) can be applied by setting \code{dat.type} to \code{"mass"}.
  In this case, to determine total cumulative biogas production for batch reactors, an initial mass should be included for each reactor.
  All volumes (input and output) have the same units, e.g., mL, L, SCF.
  To surpress volume correction to a 'standard' temperature and pressure, leave \code{temp} as \code{NULL} (the default).
}
\value{
  \code{cumBg} returns a data frame with all the columns originally present in \code{dat}, plus these others:
  \item{vBg}{Normalized volume of biogas production for individual event.}
  \item{xCH4}{Interpolated mole fraction of methane in biogas.}
  \item{vCH4}{Normalized volume of methane production for individual event. Only if \code{comp} is provided.}
  \item{cvBg}{Normalized cumulative volume of biogas production.}
  \item{cvCH4}{Normalized cumulative volume of methane production. Only if \code{comp} is provided.}
  \item{rvBg}{Production rate of biogas.}
  \item{rvCH4}{Production rate of methane. Only if \code{comp} is provided.}
}
%NTS: are there any references for these methods?
\references{
  Hafner, S.D., Rennuit, C., Triolo, J.M., Richards, B.K. In review. A gravimetric method for measuring biogas production. \emph{Biomass and Bioenergy}.

  Richards, B.K., Cummings, R.J., White, T.E., Jewell, W.J. 1991. Methods for kinetic analysis of methane fermentation in high solids biomass digesters. \emph{Biomass and Bioenergy} 1: 65-73.
}
\author{
  Sasha D. Hafner and Charlotte Rennuit
}

\seealso{
  \code{\link{summBg}},
  \code{\link{interp}},
  \code{\link{stdVol}}
}
\examples{
data("vol")
data("comp")
head(vol)
head(comp)
# Example using extrapolation
cumBg(vol,temp=20,comp=comp,time.name="days",extrap=TRUE)
# Compare to extrap=FALSE (default)
cumBg(vol,temp=20,comp=comp,time.name="days")
# Time 0 was assumed to be at days=0 (default addt0 argument), can show this row
cumBg(vol,temp=20,comp=comp,time.name="days",extrap=TRUE,showt0=TRUE)
# Leave out composition data
cumBg(vol,temp=20,time.name="days",extrap=TRUE)

# Example with different column names
data(vol2)
data(comp2)
head(vol2)
head(comp2)
cumBg(vol2,temp=20,comp=comp2,id.name="bottle",time.name="days",dat.name="meas.vol",
      comp.name="CH4.conc",extrap=TRUE)

# Mass example
data("mass")
mass
cumBg(mass,dat.type="mass",temp=38,comp=comp,time.name="days")
cumBg(mass,dat.type="mass",temp=38,comp=comp,time.name="days",showt0=FALSE)
}
%% Add one or more standard keywords, see file "KEYWORDS" in the
%% R documentation directory.
\keyword{chron}
\keyword{manip}
\concept{biogas}
