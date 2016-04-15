\name{summBg}
\alias{summBg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summarise and Normalise Cumulative Methane Production
}
\description{
    From cumulative gas production, use \code{summBg} to standardise, interpolate, subtract innoculum contribution, normalise by substrate mass, and summarise the output calculating mean and standard devations for each type of sample (groups of replicates). 
    The function is flexible: some, all, or none of these operations can be carried out in a call.
    Typically \code{summBg} is used to calculate biochemical methane potential (BMP) from cumulative methane production.
  }
\usage{
summBg(vol, setup, id.name = "id", time.name = "time", 
       descrip.name = "descrip", inoc.name = NULL, 
       norm.name = NULL, inoc.m.name = "minoc", vol.name = "cvCH4", 
       imethod = "linear", extrap = FALSE, when = 30, 
       show.obs = FALSE, sort = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vol}{
    a data frame with the columns reactor identification code; time of measurement (as \code{numeric}, or \code{POSIX}); and the response variable biogas volume or cumulative volume.
    The names of these columns can be specified with \code{id.name}, \code{time.name}, and \code{vol.name}. Or, the default names can be used. 
    Additional columns can be present-- these will be returned in the output data frame.
}
  \item{setup}{ %maybe could be shorned by just giving the column names (descriptions) of the columns are afterwards
    a data frame containing information to summarise, substract inoculum effect or normalise the data. Should contain at least the column \code{id.name} (reactor identification code), and the additional columns depending on the operations desired. See \sQuote{Details} for more information.
}
  
  \item{id.name}{
    name of the reactor identification code column in \code{vol} and \code{setup}. Should be the same in all data frames used in the function. Default is \code{"id"}.
}
  \item{time.name}{
    name of column containing time data in \code{vol}. Default is \code{"time"}.
}
  \item{descrip.name}{
    (optional) name of column containing a description of reactor substrate (or a code for this) in \code{setup} data frame. The summary of the data will be made following this column. Default is \code{"descrip"}.
}
  \item{inoc.name}{
    (optional) the value in the \code{descrip.name} (\code{setup} data frame) column that is used to describe the inoculum-only reactors. Length one character vector. Argument not used by default.  
}
  \item{norm.name}{
    (optional) the name of the column in \code{setup} that has the mass to perform the normalisation (typically substrate mass). Length one character vector. Argument not used by default.
}
  \item{inoc.m.name}{
    (optional) the name of the column in \code{setup} that has inoculum mass.
    Length one character vector. Default is \code{"minoc"}.
}
  \item{vol.name}{
    the name of the column in \code{vol} that contains the response variable --typically cumulative methane volume.
    The default is \code{"cvCH4"} for \emph{c}umulative \emph{v}olume of methane (CH4).
}
  \item{imethod}{
    the interpolation method to be used. 
    This is passed as the \code{method} argument to \code{\link{interp}}. Length one character vector. Default is \code{"linear"} for linear interpolation.
}
  \item{extrap}{
    should extrapolation be carried out? Set to \code{TRUE} if extrapolation wanted.
    Length one logical vector.
    This is passed as the \code{extrap} argument to \code{\link{interp}}. Default is \code{FALSE}
}
  \item{when}{
    value(s) of \code{time.name} in \code{vol} data frame at which the results should be evaluated.
    Alternatively, set to \code{"end"} (to return the latest values for each reactor) or \code{"meas"} (to return a value for each measured time for each reactor).
    Numeric or character vector.
}
  \item{show.obs}{
    set to \code{TRUE} to return all observations, otherwise means and standard deviations are returned. Default is \code{FALSE}. 
}
  \item{sort}{
    controls sorting of results, which is by \code{descrip.name} and then \code{time.name} by default. Set to \code{FALSE} to use order of \code{descrip.name} from the \code{setup} data frame. Default is \code{TRUE}.
}
}
\details{
  \code{summBg} was primarily designed to calculate the biochemical methane potential (BMP) from cumulative methane production of a set of batch reactors through these steps: interpolation of cumulative production to a specified time (if needed), subtratction of apparent innoculum contribution, normalisation of the results by substrate mass (typically volatile solids (VS) mass, but could be the mass of anything within the reactor) and calculation of mean and standard deviation for each sample type (set of replicates, identified by \code{descrip.name}). 

If needed \code{summBg} can return values for all observations and be used for simpler operations e.g., determining cumulative biogas production at some specified time or normalising gas volume by different substrate characteristics.

To summarise data, the \code{setup} data frame should have a column with a description of reactor substrate (or a code for this). 
Its name is set by \code{descrip.name}.  
If the inoculum effect is to be subtracted out, a column named \code{descrip.name} (for identifying replicates) and a column with the mass of inoculum present (any units), with a name set by \code{inoc.m.name} are both required in \code{setup}. 
To normalise by substrate mass (or any mass to be used for normalisation), an additional column with the mass of substrate is needed in \code{setup} --its name is set by \code{norm.name}.
  
This function is probably easier to understand by example.
See \sQuote{Examples}.
}
\value{
  a data frame, with the colums:
  \item{descrip}{from the input data frames}
  \item{mean}{mean of the response variable \code{vol.name}}
  \item{sd}{standard deviation of the response variable \code{vol.name}}
  \item{n}{number of reactors}
}
\author{
  Sasha D. Hafner and Charlotte Rennuit
}
\note{
  Reported standard deviation includes an estimate of variability from subtracting the inoculum contribution when this is done.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{cumBg}}
}
\examples{
data("vol")
data("comp") 
data("setup") 

# First need to calculate cumulative methane production data
cum.prod <- cumBg(vol, comp = comp, temp = 20, pres = 1, 
                  time.name = "days", extrap = TRUE)

head(cum.prod)

# Cumulative methane production at 30 d (default)
summBg(vol = cum.prod, setup = setup, time.name = "days", when = 30)

# Or total cumulative biogas
summBg(vol = cum.prod, setup = setup, time.name = "days", 
       vol.name = "cvBg", when = 30)

# Cumulative CH4 only, subtract inoculum contribution
summBg(vol = cum.prod, setup = setup, time.name = "days", 
       inoc.name = "inoc", when = 30)

# And normalise by mvs.sub column (mass of substrate VS here) (so the result is BMP)
summBg(vol = cum.prod, setup = setup, time.name = "days", 
       inoc.name = "inoc", when = 30, norm.name = "mvs.sub")

# Same example, but return results for three times
summBg(vol = cum.prod, setup = setup, time.name = "days", 
       inoc.name = "inoc", when = c(10, 30, 60), 
       norm.name = "mvs.sub")

# Back to earlier example, but return all individual observations 
# (and total production and individual contributions of substrate 
# and inoculum) 
summBg(vol = cum.prod, setup = setup, time.name = "days", 
       inoc.name = "inoc", when = 30, 
       norm.name = "mvs.sub", show.obs = TRUE)

# Something different: interpolated biogas production rates
summBg(vol = cum.prod, setup = setup, time.name = "days", 
       vol.name = "rvBg", when = 30, show.obs = TRUE)

# The when argument could also be 'meas' or 'end' for any of these examples 
summBg(vol = cum.prod, setup = setup, time.name = "days", 
       inoc.name = "inoc", when = "end", 
       norm.name = "mvs.sub", show.obs = TRUE)

# Example with dataset with different column names
data("vol2")
data("comp2")
data("setup2")

# First need to calculate cumulative methane production data
cum.prod <- cumBg(vol2, comp = comp2, temp = 20, pres = 1, 
                  id.name = "bottle", time.name = "days", 
		  dat.name = "meas.vol", comp.name = "CH4.conc", 
                  extrap = TRUE)

head(cum.prod)

# Cumulative CH4 production at 30 d, subtract inoculum contribution 
# and normalise by sub.vs column (mass of substrate VS here) (look at setup2).
summBg(vol = cum.prod, setup = setup2, id.name = "bottle", 
       time.name = "days", descrip.name = "description", 
       inoc.name = "Inoculum", inoc.m.name = "inoc.mass", 
       norm.name = "sub.vs", when = 30)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{manip}
\keyword{chron}
\concept{biogas}
