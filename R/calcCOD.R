# Modified: 22 JULY 2015

calcCOD <-
function(form) {

  # Check argument
  checkArgClassValue(form, 'character')

  # If and only if first letter of form is lowercase, entire string is capitalized
  if(grepl('[a-z]', substr(form, 1,1))) form <- toupper(form)

  # Read chemical formula
  fc <- readFormula(form, elements = c('C', 'H', 'O', 'N'))
  mmass <- molMass(form)

  # Calculate COD based on Rittmann and McCarty
  COD <- as.vector((2*fc['C'] + 0.5*fc['H'] - 1.5*fc['N'] - fc['O'])*molMass('O')/mmass)

  return(COD)
}
