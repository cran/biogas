# Modified: 4 MAR 2015

predBg <-
function(
  form,     # Character chemical formula of substrate
  mass=1,   # Mass of substrate (g)
  mol,      # Moles of substrate 
  fs=0,     # Fraction substrate going to cell synthesis, fs in Rittmann and McCarty
  fd=1,     # Biodegradable fraction of substrate (g/g)
  mcomp,    # Macromolecular composition, vector with named elements NTS add details
  COD,      # Substrate total COD (g) (alternative to form or mcomp)
  conc.sub, # Substrate concentration in g per kg water (only needed for CO2 partitioning)
  pH,       # Solution pH (only needed for CO2 partitioning) 
  temp,     # Temperature in degrees C (only needed for CO2 partitioning) 
  mu=0.1,   # Solution ionic strength (mol/kgw) (only needed for CO2 partitioning)
  value='CH4' # Output format, default is CH4 volume (mL)
  ){
  # Warnings
  # Check for errors in actual arguments
  
  # fe from Rittmann and McCarty (2001)
  fe<-1 - fs
  
  # If COD is given, just return CH4 volume in mL
  if(!missing(COD)) {
    vCH4<-fd*fe*COD*22361*0.5/(15.9994*2)
    if(tolower(value)=='ch4') {
      return(vCH4)
    } else {
      return(data.frame(COD=COD,fd=fd,fs=fs,fe=fe,vCH4))
    }
  } 

  # Read chemical formula
  if(!missing(form)) {
    # Capitalize form
    form<-toupper(form)
    fc<-t(mapply(readFormula,form=form,MoreArgs=list(elements=c('C','H','N','O'),min.elements=c('C','H')),USE.NAMES=FALSE))
  } else if(!missing(mcomp)) {
    # Make sure it is a vector. Cannot provide multiple values for mcomp.
    if(!is.vector(mcomp)) stop('mcomp can only be a vector, and can only be used to specify a single composition.')
    # Check names, add if missing
    mm<-c(vfa=0,protein=0,carbohydrate=0,lipid=0,lignin=0)
    mm[names(mcomp)]<-mcomp
    mcomp<-mm
    if(abs(sum(mcomp)-1)>1E-5) {
      warning('Sum of mcomp != 1.0 so dividing all elements by the sum.')
      mcomp<-mcomp/sum(mcomp)
    }
    names(mcomp)<-tolower(names(mcomp))
    fc<-mcomp['vfa']*c(C=2,H=4,O=2,N=0) + mcomp['protein']*c(C=5,H=7,O=2,N=1) + mcomp['carbohydrate']*c(C=6,H=10,O=5,N=0) + mcomp['lipid']*c(C=57,H=104,O=6,N=0) + mcomp['lignin']*c(C=10,H=13,O=3,N=0)
    fc<-round(fc/min(fc[fc!=0]))
    form<-paste0(names(fc)[fc!=0],fc[fc!=0],collapse='') # So we have a text-only empirical formula for writing out
    # Turn into a matrix, to be compatible with code below for vectorized formulas (mcomp is not vectorize)
    fc<-t(as.matrix(fc))
  } else stop('Must provide one of these arguments: form, mcomp, or COD.')

  # Molar mass of substrate
  mmass<-mapply(molMass,form,USE.NAMES=FALSE)

  # Moles and mass of substrate, if mol is provided, it will override default mass
  if(missing(mol)) mol<-mass/mmass else mass<-mol/mmass

  # Calculate COD'
  COD<-mapply(calcCOD,form,USE.NAMES=FALSE)*mass

  # Coefficients in overall reaction in Eq. 13.5 in Rittmann and McCarty (2001)
  d<-as.vector(4*fc[,'C'] + fc[,'H'] - 2*fc[,'O'] - 3*fc[,'N'])
  cH2O<-as.vector(2*fc[,'C'] + fc[,'N'] - fc[,'O'] - 9*d*fs/20 - d*fe/4)
  cCH4<-as.vector(d*fe/8)
  cCO2<-as.vector(fc[,'C'] - fc[,'N'] - d*fs/5 - d*fe/8)
  cbio<-as.vector(d*fs/20)
  cNH4<-as.vector(fc[,'N'] - d*fs/20)
  cHCO3<-as.vector(fc[,'N'] - d*fs/20)

  # Results
  # CH4 in mL at 0C and 1 atm
  nCH4<-cCH4*fd*mol
  vCH4<-22361*nCH4

  if(toupper(value)=='CH4') {
    return(vCH4)
  }

  # Hydrolytic water consumption, g H2O
  h<-cH2O*18.0152*fd*mol # NTS double check
  # Ammonia requirement, g N substrate (if < 0 then is produced not required)
  mNH4<- -cNH4*14.0067*fd*mol # NTS double check
  # Microbial biomass production
  mBio<-cbio*113*fd*mol
  # CH4, CO2, and biogas in g
  mCH4<-16.0423*nCH4
  nCO2<-cCO2*fd*mol
  mCO2<-44.0095*nCO2
  # Biogas composition (per production, but does not consider TIC in solution)
  fCH4<-nCH4/(nCH4 + nCO2)

  # Partition CO2 between headspace and solution if H2O mass, pH, and temperature are available
  if(!missing(conc.sub) & !missing(pH) & !missing(temp)) {
    # Note that method is hard-wired to 'cont' for now. 
    out<-t(mapply(partCO2,nCO2=nCO2,nCH4=nCH4,mass.H2O=mass/conc.sub,temp.c=temp,pH=pH,mu=mu,pres=1,method='cont',value='all',USE.NAMES=FALSE))
    nCO2Bg<-out[,'nCO2Bg']
    nCO2.sol<-out[,'nCO2.sol']
    cTIC<-out[,'cTIC']
    xCH4<-out[,'xCH4']

    # Volume of CO2 left in biogas
    vCO2<-22236*nCO2Bg
    # Dry biogas volume
    vBg<-vCH4 + vCO2

    # Masses of CO2 in solution and biogas
    mCO2Bg<-44.0095*nCO2Bg
    mCO2.sol<-44.0095*nCO2.sol

    if(all(fs==0)) {
      results<-data.frame(form=form,mass=mass,mol.mass=mmass,moles=mol,COD=COD,fd=fd,conc.sub=conc.sub,temp=temp,pH=pH,hydro=h,fCH4=fCH4,xCH4=xCH4,vCH4=vCH4,vCO2=vCO2,vBg=vBg,mCH4=mCH4,mCO2=mCO2,mCO2Bg=mCO2Bg,mCO2.sol=mCO2.sol,cTIC=cTIC)
    } else {
      results<-data.frame(form=form,mass=mass,mol.mass=mmass,moles=mol,COD=COD,fs=fs,fe=fe,fd=fd,conc.sub=conc.sub,temp=temp,pH=pH,hydro=h,fCH4=fCH4,xCH4=xCH4,vCH4=vCH4,vCO2=vCO2,vBg=vBg,mCH4=mCH4,mCO2=mCO2,mCO2Bg=mCO2Bg,mCO2.sol=mCO2.sol,cTIC=cTIC,m.bio=mBio,N.req=mNH4)
    }
    rownames(results)<-1:nrow(results)
    return(results)
  }

  # Not enough info for CO2 partitioning
  if(all(fs==0)) {
    results<-data.frame(form=form,mass=mass,mol.mass=mmass,moles=mol,COD=COD,hydro=h,fCH4=fCH4,vCH4=vCH4,mCH4=mCH4,mCO2=mCO2)
  } else {
    results<-data.frame(form=form,mass=mass,mol.mass=mmass,moles=mol,COD=COD,fs=fs,fe=fe,fd=fd,hydro=h,fCH4=fCH4,vCH4=vCH4,mCH4=mCH4,mCO2=mCO2,m.bio=mBio,N.req=mNH4)
  }
  rownames(results)<-1:nrow(results)
  return(results)
}
