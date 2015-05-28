mass2vol <-
function(
	mass,
	xCH4,
	temp,pres=1,
	unit.pres="atm",
	unit.temp="C",
	pres.std=NULL,
	temp.std=NULL,
	value="CH4",
	headspace=NULL,
	headcomp='N2'
	) {

  # Hardwire rh, maybe will add as argument again later
  rh<-1

  # Convert pressure to atm and temp to C
  pres<-unitConvert(x=pres,unit=unit.pres,to='atm')
  temp<-unitConvert(x=temp,unit=unit.temp,to='C')

  if(is.null(pres.std)) pres.std<-1.0 else pres.std<-unitConvert(x=pres.std,unit=unit.pres,to='atm')
  if(is.null(temp.std)) temp.std<-0.0 else temp.std<-unitConvert(x=temp.std,unit=unit.temp,to='C')
 
  mmb<-xCH4*16.0425 + (1 - xCH4)*44.0095
  mvBg<-xCH4*22361 + (1 - xCH4)*22263
  db<-mmb/mvBg

  # Calculate water vapor pressure in atm (based on NIST)
  pH2O<-rh*10^(6.203913 - 2354.731/(temp + 280.709))

  mH2O<-18.0152*pH2O/((pres - pH2O)*mvBg)

  # Biogas volume
  vBg<-mass/(db + mH2O)

  # If headspace data provided
  if(!is.null(headspace)) {
    if(headcomp=='N2') {
      # Reduce mass loss for N2 loss
      mass<-mass - stdVol(headspace,temp=20,pres=1,rh=0,temp.std=0,pres.std=1,unit.pres='atm',unit.temp='C')*0.0012504 # NTS: fixed T and P now, need to improve. Assumes dry headspace initially.
      vBg<-mass/(db + mH2O) + stdVol(headspace,temp=temp,pres=1.0,rh=1,temp.std=0,pres.std=1,unit.pres='atm',unit.temp='C') # NTS: assumes final residual pressure of 1.0 atm
    } else {
      warning('headcomp argument was given as ',headcomp,', but \"N2\" is the only option available, so no correction was applied for initial headspace composition.')
    }
  }

  # Normalize (is 1 atm 0C by default, based on molar volumes used above--so in this case stdVol does nothing.)
  vBg<-stdVol(vBg,temp=0,pres=1,rh=0,temp.std=temp.std,pres.std=pres.std,unit.pres='atm',unit.temp='C')
  vCH4<-xCH4*vBg*22361/mvBg
  vCO2<-(1 - xCH4)*vBg*22263/mvBg

  if(tolower(value)=="ch4") {
    return(vCH4) 
  } else if(tolower(value)=="bg") {
    return(vBg) 
  } else if(tolower(value)=="co2") {
    return(vCO2) 
  } else {
    return(cbind(vBg=vBg,vCH4=vCH4,vCO2=vCO2))
  }
}
