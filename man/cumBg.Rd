\name{cumBg}
\alias{cumBg}
%- Also NEED an "\alias" for EACH other topic documented here.
\title{
  Calculate Cumulative Biogas Production
}
\description{
  \code{cumBg} (for \emph{cum}ulative \emph{b}io\emph{g}as) calculates cumulative biogas, methane production and production rates from individual volume (or mass) and composition measurements for any number of reactors.
}
\usage{
cumBg(dat, dat.type = 'vol', comp = NULL, temp = NULL, pres = NULL,
  id.name = 'id', time.name = 'time', dat.name = dat.type, comp.name = 'xCH4', 
  headspace = NULL, vol.hs.name = 'vol.hs', headcomp = 'N2', temp.init = NULL,
  temp.std = getOption('temp.std', 0.0), pres.std = getOption('pres.std', 1.0),
  unit.temp = getOption('unit.temp', 'C'), unit.pres = getOption('unit.pres', 'atm'),
  cmethod = 'removed', imethod = 'linear', extrap = FALSE,
  addt0 = TRUE, showt0 = TRUE, std.message = TRUE, check = TRUE) 
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{dat}{
    a data frame with the columns reactor identification code; time of measurement (as \code{numeric}, or \code{POSIX}); and measured biogas volume or total reactor mass, in g (see \code{dat.type} argument).
    The names of these columns are specified with \code{id.name}, \code{time.name}, and \code{dat.name}.
    Additional columns can be present--these will be returned in the output data frame.
}
  \item{dat.type}{
  the type of data contained in \code{dat}. 
  Use \code{"vol"} for biogas volume (volumetric method will be used) or \code{"mass"} for reactor masses (gravimetric method will be used).
  Default is \code{"vol"}.
}
  \item{comp}{
    (optional) a data frame with the columns reactor identification code; time of measurement, (as \code{numeric}, or \code{POSIX}); and methane concentration within dry biogas as a mole fraction, considering only methane and carbon dioxide (unless \code{cmethod = "total"}).
    If omitted, cumulative biogas volume will still be calculated and returned (but no methane data will be returned).
    The names of these columns are specified with \code{id.name}, \code{time.name}, and \code{comp.name}.
    Default is \code{NULL}.
}
  \item{temp}{
    the temperature at which biogas volume was measured (when \code{dat.type = "vol"}), or of biogas just prior to exiting the reactor (when \code{dat.type = "mass"}).
    A length-one numeric vector.
    Degrees Celcius by default (see \code{unit.temp} argument). 
    Default is \code{NULL}, which suppresses correction for temperature and pressure. 
}
  \item{pres}{
    the pressure at which biogas volume was measured (when \code{dat.type = "vol"}), or of biogas just prior to exiting the reactor (when \code{dat.type = "mass"}).
    A length-one numeric vector.
    Atmosphere by default (see \code{unit.pres} argument). 
    Default is \code{NULL}, which suppresses correction for temperature and pressure. 
}
  \item{id.name}{
    name of the reactor identification code column in \code{dat}. Should be the same in all data frames used in the function. Default is \code{"id"}.
}
  \item{time.name}{
    name of column containing time data in \code{dat} and \code{comp} data frames. Default is \code{"time"}.
}
  \item{dat.name}{
    name of column containing the primary response variable (volume or mass) in \code{dat} data frame. Default is \code{dat.type} value. See \code{dat.type} argument.
}
  \item{comp.name}{
    name of column containing biogas mole fraction of methane in \code{comp} data frame. Default is \code{"xCH4"}. Must be normalised so xCH4 + xCO2 = 1.0 unless \code{cmethod = "total"}.
}
  \item{headspace}{
    (optional) a data frame or numeric vector with reactor headspace volumes.
    If a data frame is used, it should at least contain a \code{"id"} (reactor identification code) column (see \code{"id.name"}) and headspace volume column (see \code{vol.hs.name} argument).
    Required if \code{method = "total"} for the volumetric method, and for initial headspace correction for the gravimetric method (see \code{headcomp} and \code{temp.init}).
    Default is \code{NULL}.
}
  \item{vol.hs.name}{
    name of column containing headspace volume data in optional \code{headspace} dataframe. 
    Default is \code{"vol.hs"}.
}
  \item{headcomp}{
    Composition of headspace for correction of initial headspace mass for gravimetric method only.
    If provided, correction will be applied to first mass loss observation for each individual reactor.
    See argument of the same name in \code{\link{mass2vol}}.
    Currently, the only option is \code{"N2"} for dinitrogen gas.
    Use of any other value will be ignored with a warning.
}
  \item{temp.init}{
    (optional) initial headspace temperature in \code{unit.temp} units.
    Used to correct for effect of initial reactor headspace on mass loss.
}
  \item{temp.std}{
    standard temperature for presentation of biogas and methane results.
    Length one numeric vector.
    Default value is 0 degrees C.
    Argument is passed to \code{\link{stdVol}}.
}
  \item{pres.std}{
    standard pressure for presentation of biogas and methane results.
    Length one numeric vector.
    Default value is 1.0 atm (101325 Pa).
    Argument is passed to \code{\link{stdVol}}.
}
   \item{unit.temp}{
    temperature units for \code{temp} and \code{temp.std} arguments.
    %Options are \code{"C"} (degrees Celcius), \code{"F"} (degrees Fahrenheit), and \code{"K"} (Kelvin). 
    Default is \code{"C"}. 
    Argument is passed to \code{\link{stdVol}}.
}
    \item{unit.pres}{
    pressure units for \code{pres} and \code{pres.std} arguments.
    %Options are \code{"atm"}, \code{"Pa"}, \code{"kPa"}, \code{"hPa"}, and \code{"bar"}. 
    Default is \code{"atm"}.
    Argument is passed to \code{\link{stdVol}}.
}
  \item{cmethod}{
    method for calculating cumulative methane production.
    Use \code{"removed"} to base production on \code{xCH4} and gas volumes removed (default).
    Use \code{"total"} to base it on the sum of methane removed and methane remaining in the reactor headspace.
    For \code{"removed"}, \emph{\code{xCH4} should be calculated based on methane and CO2 only (xCH4 + xCO2 = 1.0).}
    For \code{"total"}, \emph{\code{xCH4} should be calculated including all biogas components (CH4, CO2, N2, etc.) except water.}
    Length one character vector. 
    %Default is \code{"removed"}.
}
  \item{imethod}{
    method used for interpolation of \code{xCH4}.
    This is passed as the \code{method} argument to \code{\link{interp}}.
    Length one character vector. 
    Default is \code{"linear"} for linear interpolation.
}
  \item{extrap}{
    should \code{comp.name} be extrapolated? 
    %Set to \code{TRUE} if extrapolation wanted. 
    Length one logical vector.
    This is passed as the \code{extrap} argument to \code{\link{interp}}. 
    Default is \code{FALSE}.
}
 \item{addt0}{
    is the earliest time in \code{dat} data frame \dQuote{time zero} (start time)? 
    If not, this argument adds a row with \code{time.name = 0} for each reactor in order to calculate production rates for the first observation.
    This addition is only made when \code{dat.type = "vol"} and \code{time.name} is \code{numeric} (or \code{integer}).
    Length-one logical vector.
    Default is \code{TRUE}.
    To return these additional rows in the output, see \code{showt0}.
}
  \item{showt0}{
    should \dQuote{time zero} rows be returned in the output?
    Only applies if \code{time.name} is \code{numeric} (or \code{integer}).
    These rows may have been present in the original data (\code{dat}) or added by the function (see \code{addt0}).
    Default value depends on \code{dat} \code{time.name} column content. 
    If \code{time.name} column is numeric and contains 0 then the default value is \code{TRUE} and otherwise \code{FALSE}. 
}

  \item{std.message}{
    should a message with the standard conditions be displayed? 
    Default is \code{TRUE}. 
  }
  
  \item{check}{
    should input data be checked for unreasonable values (with warnings)? 
    Default is \code{TRUE}.
  }
    
}
\details{
  Using volume or mass data from \code{dat} and gas composition from \code{comp}, this function will calculate standardised biogas and methane production (if \code{comp} is provided) for each observation,  interpolating \code{comp.name} (from \code{comp} argument) to each \code{time.name} in \code{dat} if needed, and summing these for cumulative values.
  Use of volumetric measurements is the default case. 
  Alternatively, if reactor mass is given in \code{dat}, a gravimetric approach (Hafner et al. in review) can be applied by setting \code{dat.type} to \code{"mass"}.
  In this case, to determine total cumulative biogas production for batch reactors, an initial mass should be included for each reactor.

  For the volumetric method, all volumes (input and output) have the same units, e.g., mL, L, SCF.
  For the gravimetric method, output volumes will be in mL if input masses are in g.
  Biogas composition (\code{comp.name} column in \code{comp}) is specified as the mole fraction of methane in dry biogas, normalised so mole fractions of methane and carbon dioxide sum to unity (Richards et al. 1991).
 
  Standard values and units for temperature and pressure can be globally set using the function \code{\link{options}}.
  See \code{\link{stdVol}}.
  To surpress volume correction to a \dQuote{standard} temperature and pressure, leave \code{temp} as \code{NULL} (the default).

  If \code{check = TRUE}, the input values of mole fraction of methane in biogas are checked, and a warning is returned if the are outside 0, 1.
  
}


\value{
  a data frame with all the columns originally present in \code{dat}, plus these others:
  \item{vBg}{Standardised volume of biogas production for individual event.}
  \item{xCH4}{Interpolated mole fraction of methane in biogas.}
  \item{vCH4}{Standardised volume of methane production for individual event. Only if \code{comp} is provided.}
  \item{cvBg}{Standardised cumulative volume of biogas production.}
  \item{cvCH4}{Standardised cumulative volume of methane production. Only if \code{comp} is provided.}
  \item{rvBg}{Production rate of biogas.}
  \item{rvCH4}{Production rate of methane. Only if \code{comp} is provided.}
}
\references{
  Hafner, S.D., Rennuit, C., Triolo, J.M., Richards, B.K. In review. Validation of a simple gravimetric method for measuring biogas production in laboratory experiments. \emph{Biomass and Bioenergy}.

  Richards, B.K., Cummings, R.J., White, T.E., Jewell, W.J. 1991. Methods for kinetic analysis of methane fermentation in high solids biomass digesters. \emph{Biomass and Bioenergy} 1: 65-73.
}
\author{
  Sasha D. Hafner and Charlotte Rennuit
}

\seealso{
  \code{\link{summBg}},
  \code{\link{interp}},
  \code{\link{stdVol}},
  \code{\link{options}}
}

\examples{

data("vol")
data("comp")

head(vol)
head(comp)

# Calculate cumulative production and rates from vol and comp
cum.prod <- cumBg(vol, comp = comp, temp = 20, pres = 1, id.name = "id", time.name = "days", 
		  comp.name = "xCH4", dat.name = "vol")
head(cum.prod)

# Note warnings and related NAs in results
# Set extrap = TRUE to extrapolate xCH4 to earliest times

# Calculate cumulative production and rates from vol and comp
cum.prod <- cumBg(vol, comp = comp, temp = 20, pres = 1, id.name = "id", time.name = "days", 
		  comp.name = "xCH4", dat.name = "vol", extrap = TRUE)
head(cum.prod)

# In this case, we can use default values for some column names, so this call is identical
cum.prod <- cumBg(vol, comp = comp, temp = 20, pres = 1, time.name = "days", extrap = TRUE)

# Plot results
library(ggplot2)
qplot(x = days, y = cvCH4, data = cum.prod, xlab = "Time (d)", 
      ylab = "Cumulative methane production (mL)",color = id, geom = "line")

# Omit added time zero rows
cum.prod <- cumBg(vol, comp = comp, temp = 20, pres = 1, time.name = "days", extrap = TRUE, 
		  showt0 = FALSE)
head(cum.prod)

qplot(x = days, y = cvCH4, data = cum.prod, xlab = "Time (d)", 
      ylab = "Cumulative methane production (mL)", color = id, geom = "line")

# Previous is different from never adding them in the first place (rates not calculated for first 
# observations here)
cum.prod <- cumBg(vol, comp = comp, temp = 20, pres = 1, time.name = "days", extrap = TRUE, 
		  addt0 = FALSE)
head(cum.prod)

qplot(x = days, y = cvCH4, data = cum.prod, xlab = "Time (d)", 
      ylab = "Cumulative methane production (mL)", color = id, geom = "line")

# Can use POSIX objects for time (but cumBg cannot add t0 rows here)
class(vol$date.time)
class(comp$date.time)
cum.prod <- cumBg(vol, comp = comp, temp = 20, pres = 1, time.name = "date.time", extrap = TRUE)
head(cum.prod)

qplot(x = date.time, y = cvCH4, data = cum.prod, xlab = "Time (d)", 
      ylab = "Cumulative methane production (mL)", color = id, geom = "line")

# Can leave out composition data, and then CH4 is not included in results
cum.prod <- cumBg(vol, temp = 20, pres = 1, time.name = "days")
head(cum.prod)

# Leave out pres or temp, and results are not standardised
cum.prod <- cumBg(vol, time.name = "days")
head(cum.prod)

# Example with input data frames with different column names
data("vol2")
data("comp2")

head(vol2)
head(comp2)

cum.prod <- cumBg(vol2, comp = comp2, temp = 20, pres = 1, id.name = "bottle", time.name = "days", 
		  dat.name = "meas.vol", comp.name = "CH4.conc")
head(cum.prod)
tail(cum.prod)

# Note warnings and related NAs in results
warnings()

# Set extrap = TRUE to avoid
cum.prod <- cumBg(vol2, comp = comp2, temp = 20, pres = 1, id.name = "bottle", time.name = "days", 
		  dat.name = "meas.vol", comp.name = "CH4.conc", extrap = TRUE)
head(cum.prod)

qplot(x = days, y = cvCH4, data = cum.prod, xlab = "Time (d)", 
      ylab = "Cumulative methane production (mL)", color = bottle, geom = "line")

# Mass example
data("mass")

mass

# Need to specify data type with dat.type argument (using default values for id.name, dat.name, 
# and comp.name)
cum.prod <- cumBg(mass, dat.type = "mass", comp = comp, temp = 35, pres = 1, time.name = "days")
cum.prod

# Drop time 0 rows
cum.prod <- cumBg(mass, dat.type = "mass", comp = comp, temp = 35, pres = 1, time.name = "days", 
		  showt0 = FALSE)
cum.prod

# Add initial headspace correction (alternatively, headspace could be a data frame with a different 
# volume for each reactor)
cum.prod <- cumBg(mass, dat.type = "mass", comp = comp, temp = 35, pres = 1, time.name = "days", 
		  headspace = 300, headcomp = "N2", temp.init = 20, showt0 = FALSE)
cum.prod

}
%% Add one or more standard keywords, see file "KEYWORDS" in the
%% R documentation directory.
\keyword{chron}
\keyword{manip}
\concept{biogas}
