\name{biogas-package}
\alias{biogas-package}
\alias{biogas}
\docType{package}
\title{
Analyze Biogas Data and Predict Biogas Production
}
\description{
Functions for working with biogas data. Both low- and high-level functions are included for carrying out common tasks for analysis of biogas and related data. Molar mass and calculated oxygen demand (COD') can be determined from a chemical formula. Measured gas volume can be corrected for water vapor and to (possibly user-defined) standard temperature and pressure. Gas composition, cumulative production, or other variables can be interpolated to a specified time. Cumulative biogas and methane production (and rates) can be calculated from individual biogas volume (or reactor mass) and biogas composition measurements for any number of reactors. With cumulative methane production data and data on reactor contents, biochemical methane potential (BMP) can be calculated and summarized, including subtraction of the inoculum contribution and normalization by substrate mass. Cumulative production and production rates can be summarized in several different ways (e.g., omitting normalization) using the same function. Lastly, biogas and methane production can be predicted from substrate composition and additional, optional data.
}
\details{
\tabular{ll}{
Package: \tab biogas\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-05\cr
License: \tab GPL-2\cr
}
    The \code{biogas} package provides several low-level functions for working with biogas measurements and related data.
    Use the function \code{\link{stdVol}} to correct biogas volumes to 'standard' temperature and pressure, and dry conditions.
    For interpolating biogas composition or cumulative production, use \code{\link{interp}}.
    The function \code{\link{molMass}} and \code{\link{calcCOD}} can be used to calculate molar mass and chemical oxygen demand based on a chemical formula.
    And for gravimetric measurement of biogas production, there is \code{\link{mass2vol}} and \code{\link{vol2mass}}.

    Additionally, there are two high-level functions in the \code{biogas} package that simplify analysis of biogas measurements. 
    \code{\link{cumBg}} calculates cumulative biogas and methane production from separate, indvidual observations of biogas volume (or reactor mass) and biogas composition for any number of batch reactors.
    Interpolation is automatically applied as needed. 
    And \code{\link{summBg}} summarizes cumulative methane production by providing cumulative values for a specified time.
    \code{\link{summBg}} can calculate biochemical methane potential (BMP) from measurements by subtracting the contribution of an inoculum, and normalizing results by substrate mass.

    Lastly, \code{\link{predBg}} can be used to predict biogas production and composition based on substrate composition, along with degradability and other optional parameters.
    \code{\link{predBg}} will predict biogas and methane volume (and mass) based on an empirical formula or 'macromolecule'-based composition.
    If other parameters are provided (e.g., degradability, microbial substrate partitioning, solution pH, or substrate concentration), additional predictions can be made, e.g., biogas composition and nitrogen requirement.
}
\author{
Sasha D. Hafner, Charlotte Rennuit, Jin Mi Triolo, and Ali Heidarzadeh Vazifehkhoran.

Maintainer: Sasha D. Hafner <saha@kbm.sdu.dk>
}
\references{
  Hafner, S.D., Rennuit, C., Triolo, J.M., Richards, B.K. In review. A gravimetric method for measuring biogas production. \emph{Biomass and Bioenergy}.

  Richards, B.K., Cummings, R.J., White, T.E., Jewell, W.J. 1991. Methods for kinetic analysis of methane fermentation in high solids biomass digesters. \emph{Biomass and Bioenergy} 1: 65-73.

  Rittmann, B., McCarty, P. 2001 \emph{Environmental Biotechnology}. McGraw-Hill.
}
\keyword{package}
\examples{
# Calculate BMP with data from a laboratory BMP test, and calculate and compare results for 
# cellulose to the expected value
data("vol2")
data("comp2")
data("setup2")

# vol2 has measured biogas volumes
head(vol2)

# comp2 has measured biogas composition
head(comp2)

# setup2 has information on batch reactor content
head(setup2)

# Calculate cumulative biogas and methane production volumes at dry conditions and 
# 0 degrees C and 1.0 atm
cprod <- cumBg(vol2,temp=20,comp=comp2,id.name="bottle",time.name="days",dat.name="meas.vol",
               comp.name="CH4.conc",extrap=TRUE)
head(cprod)

# Now we can get methane production from a common time for all reactors (60 d here) 
# (interpolation is used if needed), calculate and subtract out the inoculum contribution, 
# and normalize by substrate VS mass
bmps <- summBg(vol=cprod,setup=setup2,when=60,id.name="bottle",time.name="days",
               descrip.name="description",inoc.name="Inoculum",inoc.m.name="inoc.mass",
               norm.name="sub.vs")
bmps

# Or, to see all cumulative values at 60 d
bmpobs <- summBg(vol=cprod,setup=setup2,when=60,id.name="bottle",time.name="days",
                 descrip.name="description",inoc.name="Inoculum",inoc.m.name="inoc.mass",
                 norm.name="sub.vs",show.obs=TRUE)
bmpobs
# When individual observations are shown like this, the separate contributions of substrate
# and inoculum are also returned, along with total production (including inoculum
# contribution).

# It is possible to get results for more than one time
bmps2 <- summBg(vol=cprod,setup=setup2,when=c(10,30,60),id.name="bottle",time.name="days",
                descrip.name="description",inoc.name="Inoculum",inoc.m.name="inoc.mass",
                norm.name="sub.vs")
bmps2

# bmps, bmps2, and bmpobs results are corrected for inoculum and normalized to 1 g substrate VS
# To just get all cumulative methane values, don't provide values for inoc.name, inoc.m.name, 
# and norm.name
allresults <- summBg(vol=cprod,setup=setup2,when=60,id.name="bottle",time.name="days",
                     descrip.name="description", show.obs=TRUE)
allresults

# How much methane should we expect to get from pure cellulose?
predBg("C6H10O5") 

# Let's assume 10% is used to produce microbial biomass
predBg("C6H10O5",fs=0.05) 
predBg("C6H10O5",fs=0.05,value="all") 

# Compare to measurement
bmps
bmps[bmps$description=="Cellulose","mean"]/predBg("C6H10O5",fs=0.05)
# About 15% higher than expected
}
\keyword{univar}
\concept{biogas}
