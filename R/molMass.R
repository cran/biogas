molMass <-
function(form) {

  # If and only if first letter of form is lowercase, entire string is capitalized
  if(grepl('[a-z]',substr(form,1,1))) form<-toupper(form)

  # Define standard atomic weights
  # These are from CIAAW (http://www.ciaaw.org/atomic-weights.htm), rounded to four digits if range was large 
  amass<-c(C=12.01,H=1.008,O=16.00,N=14.007,S=32.06,P=30.974,Na=22.990,K=39.098,Cl=35.45)

  # Get coefficients of formula
  fc<-readFormula(form)

  # Check for unidentified element
  if(any(!names(fc) %in% names(amass))) stop('One or more elements in \"form\" is not in the database. (You can add it to the amass vector if you want to modify the function.)')

  # Calculate molar mass, using names of fc for indexing
  mmass<-signif(sum(amass[names(fc)]*fc),5)
  return(mmass)
}
