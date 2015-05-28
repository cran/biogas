\name{molMass}
\alias{molMass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate Molar Mass of Compound
}
\description{
  \code{molMass} calculates the molar mass of any organic compound (and some salts) based on its chemical formula.
}
\usage{
molMass(form)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{form}{a chemical formula, as a character vector of length one, e.g., \code{"C6H12O6"} for glucose or \code{"CH3COOH"} for acetic acid. 
    Not case-sensitive for single letter elements (see 'Details').
}
}
\details{
  Atomic weights are from CIAAW, and were rounded to four digits or more if possible. 
  \code{molMass} can handle two letter elements (see 'Examples') but in this case \code{form} is case-sensitive.
}
\value{
  Length one numeric vector with molar mass (g/mol)
}
\references{
  CIAAW <http://www.ciaaw.org/atomic-weights.htm>
}
\author{
  Charlotte Rennuit and Sasha D. Hafner
}

\seealso{
  \code{\link{calcCOD}}
}
\examples{
  molMass("C6H12O6")
  molMass("CH3COOH")
  molMass("CH3CH2OH")
  molMass("CH4")
  molMass("ch4")

  # Case-sentitive for two letter elements
  molMass("NaHCO3")
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\concept{biogas}
\concept{chemistry}
