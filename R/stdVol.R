# Modified: 5 MAR 2015

stdVol <-
function(vol,temp=0,pres=1,rh=1,temp.std=NULL,pres.std=NULL,unit.pres="atm",unit.temp="C") {
  # Convert pressure to atm and temp to C
  pres<-unitConvert(x=pres,unit=unit.pres,to='atm')
  temp<-unitConvert(x=temp,unit=unit.temp,to='C')

  # Next two lines avoid problems when unit.__ is changed but __.norm isn't
  if(is.null(pres.std)) pres.std<-1.0 else pres.std<-unitConvert(x=pres.std,unit=unit.pres,to='atm')
  if(is.null(temp.std)) temp.std<-0.0 else temp.std<-unitConvert(x=temp.std,unit=unit.temp,to='C')
  
  # Check for RH problem
  if(rh>1 | rh<0) stop('rh argument must be between 0 and 1, but is ',rh,'.')

  # Calculate water vapor pressure in atm (based on NIST)
  pH2O<-rh*10^(6.203913 - 2354.731/(temp + 280.709))

  # Correct volume for water and to 1.0 atm
  vol.dry<-vol*(pres - pH2O)/pres.std

  # Correct dry volume for temperature
  vol.norm<-vol.dry*(temp.std + 273.15)/(temp + 273.15)

  return(vol.norm)
}
