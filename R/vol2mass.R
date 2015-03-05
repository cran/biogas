vol2mass <-
function(
	volBg,
	xCH4,
	temp,pres=1,
	unit.pres="atm",
	unit.temp="C"
	) {

  # Hardwire rh for now
  rh<-1

  # Convert pressure to atm and temp to C
  pres<-unitConvert(x=pres,unit=unit.pres,to='atm')
  temp<-unitConvert(x=temp,unit=unit.temp,to='C')
  
  mmb<-xCH4*16.0425 + (1 - xCH4)*44.0095
  mvb<-xCH4*22361 + (1 - xCH4)*22263
  db<-mmb/mvb

  # Calculate saturated water vapor pressure in atm (based on NIST)
  pH2O<-rh*10^(6.203913 - 2354.731/(temp + 280.709))

  mH2O<-18.0152*pH2O/((pres - pH2O)*mvb)

  mass<-volBg*(db + mH2O)

  return(mass)
}
