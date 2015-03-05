# Modified: 5 MAR 2015

cumBg <-
function(
  dat,
  dat.type='vol',
  comp=NULL,
  temp=NULL,
  pres=1,
  id.name='id',
  time.name='time',
  dat.name=dat.type,
  comp.name='xCH4', 
  headspace=NULL,
  vol.hs.name='vol.hs',
  temp.std=NULL, 
  pres.std=NULL,
  unit.pres='atm',
  unit.temp='C',
  cmethod='removed',
  imethod='linear',
  extrap=FALSE,
  addt0=TRUE,
  showt0=if(any(dat[,time.name]==0)) TRUE else FALSE
  ){

  if(!dat.type %in% c('vol','mass')) stop('Argument \"dat.type\" must be either \"vol\" or \"mass\" but you gave ',dat.type,'.')

  # Hard-wire rh for now at least
  rh<-1
  ## Check for RH problem
  #if(rh>1 | rh<0) stop('rh argument must be between 0 and 1, but is ',rh,'.')

  # Next two lines avoid problems when unit.__ is changed but __.norm isn't
  if(is.null(pres.std)) pres.std<-1.0 else pres.std<-unitConvert(x=pres.std,unit=unit.pres,to='atm')
  if(is.null(temp.std)) temp.std<-0.0 else temp.std<-unitConvert(x=temp.std,unit=unit.temp,to='C')

  # Check column names in argument data frames
  # comp needs id time xCH4 
  if(!is.null(comp) && class(comp)=='data.frame') {
    if(any(!c(id.name,time.name,comp.name) %in% names(comp))){
      message('Column identities inferred from order for comp: id = comp$',names(comp)[1],', time = comp$',names(comp)[2],', and xCH4 = comp$',names(comp)[3],'.')
      names(comp)<-c(id.name,time.name,comp.name)
    }
  }

  # Check for headspace argument if it is needed
  if(is.null(headspace) & cmethod=='total') stop('cmethod is set to \"total\" but headspace argument is not provided.')

  # Sort out composition data
  # First sort so can ignore first observation for mass data
  dat<-dat[order(dat[,id.name],dat[,time.name]),]
  dat[,comp.name]<-NA
  if(!is.null(comp) && class(comp)=='data.frame'){
    # Interpolate gas composition to times of datume measurements
    for(i in unique(dat[,id.name])) {
      if(dat.type=='mass' & nrow(dat[dat[,id.name]==i,])<2) stop('There are < 2 observations for reactor ',i,' but dat.type=\"mass\". You need at least 2 observations to apply the gravimetric method.')
      dc<-comp[comp[,id.name]==i,]
      if(nrow(dc)==0) stop('No biogas composition data for reactor ',i,' so can\'t interpolate!') 
      if(nrow(dc)>1) {
	if(dat.type=='vol') {
          dat[dat[,id.name]==i,comp.name]<-interp(dc[,time.name],dc[,comp.name],time.out=dat[dat[,id.name]==i,time.name],method=imethod,extrap=extrap)
	} else if (dat.type=='mass') {
	  # Then ignore first point, since it isn't used anyway--this is just to avoid warning with interp if extrap=FALSE
          dat[dat[,id.name]==i,comp.name][-1]<-interp(dc[,time.name],dc[,comp.name],time.out=dat[dat[,id.name]==i,time.name][-1],method=imethod,extrap=extrap)
	}
      } else { # If only one xCH4 value is available, use it for all dat obs if extrap=TRUE or times match, but warn if times don't match
        for(j in 1:nrow(dat[dat[,id.name]==i,])) {
	  if(j > 1 | dat.type=='vol') { # This just to avoid warning for first observation for mass data
            if(dc[,time.name]==dat[dat[,id.name]==i,time.name][j]) {
              dat[dat[,id.name]==i,comp.name][j]<-dc[,comp.name]
            } else {
              if(extrap) {
                dat[dat[,id.name]==i,comp.name][j]<-dc[,comp.name]
              } else {
                dat[dat[,id.name]==i,comp.name][j]<-NA
      	        warning('Not enough xCH4 data to interpolate for reactor ',i,' so results will be missing.\n If you prefer, you can use extrapolation by setting extrap=TRUE.')
              }
            }
	  }
        }
      }
    }
  } else if(!is.null(comp) && class(comp) %in% c('numeric','integer') && length(comp)==1) {
    # Or if a single value is given, use it
    dat[,comp.name]<-comp
  } else {
    # If no composition data is given, just use NA
    dat[,comp.name]<-NA 
  }

  # Volumetric
  # Function will work with vol and add columns
  if(dat.type=='vol') {
    message('Working with volume data, applying volumetric method.')
    # Column names
    # vol needs id time vol
    vol<-dat
    vol.name<-dat.name
    if(any(!c(id.name,time.name,vol.name) %in% names(vol))){
      message('Column identities inferred from order for vol: id = vol$',names(vol)[1],', time = vol$',names(vol)[2],', and vol = vol$',names(vol)[3],'.')
      names(vol)<-c(id.name,time.name,vol.name,names(vol)[-1:-3])
    } 

    # Add t0 row if requested
    # Not added if there are already zeroes present!
    if(addt0 & !class(vol[,time.name])[1] %in% c('numeric','integer')) addt0<-FALSE
    if(addt0 & !any(vol[,time.name]==0)) {
      t0<-data.frame(id=unique(vol[,id.name]),time=0)
      names(t0)<-c(id.name,time.name)
      t0[,vol.name]<-0
      vol<-rbindf(vol,t0)
    }

    # Convert pressure to atm and temp to C
    pres<-unitConvert(x=pres,unit=unit.pres,to='atm')
    temp<-unitConvert(x=temp,unit=unit.temp,to='C')
 
    # Add headspace if needed
    if(cmethod=='total') {
      if(missing(headspace)) stop('cmethod = \'total\' but headspace argument is not provided.')

      if(is.numeric(headspace)) {
       vol[,vol.hs.name]<-headspace
      } else if(is.data.frame(headspace)) {       
        # headspace needs id vol
        if(any(!c(id.name,vol.name) %in% names(headspace))){
          message('Column identities inferred from order for headspace: id = headspace$',names(headspace)[1],', vol = headspace$',names(headspace)[2],'.')
          names(headspace)<-c(id.name,vol.name,names(headspace)[-1:-2])
        }
        vol<-merge(vol,headspace,by=id.name,suffixes=c('','.hs'))
      } else stop('headspace actual argument not recognized. What is it?')
    }

    # Normalize total gas volumes
    if(!is.null(temp)) {
      vol$vBg<-stdVol(vol[,vol.name],temp=temp,pres=pres,rh=rh,pres.std=pres.std,temp.std=temp.std)
    } else {
      vol$vBg<-vol[,vol.name]
      message('No temperature provided (temp argument) so volumes are NOT normalized.')
    }

    # Calculate interval gas production
    vol$vCH4<-vol$vBg*vol[,comp.name] 
    if(addt0) vol[vol[,time.name]==0,'vCH4']<-0

    # Cumulative gas production
    vol<-vol[order(vol[,id.name],vol[,time.name]),]
    # Calculate delta t for rates
    if(class(vol[,time.name])[1] %in% c('numeric','integer')) {
      dt<-c(NA,diff(vol[,time.name]))
    } else if(class(vol[,time.name])[1] %in% c('POSIXct','POSIXlt')) {
      dt<-c(NA,as.numeric(diff(vol[,time.name]),units='days'))
    } else {
      dt<-NA
      warning('time column in vol data frame not recognized, so rates will not be calculated.')
    }
    # Set dt to NA for first observations for each reactor
    dt[c(TRUE,vol[,id.name][-1] != vol[,id.name][-nrow(vol)])]<-NA 
    
    for(i in unique(vol[,id.name])) {
      vol[vol[,id.name]==i,'cvBg'] <-cumsum(vol[vol[,id.name]==i,'vBg' ])
      vol[vol[,id.name]==i,'cvCH4']<-cumsum(vol[vol[,id.name]==i,'vCH4'])
      vol[vol[,id.name]==i,'rvBg'] <- vol[vol[,id.name]==i,'vBg' ]/dt[vol[,id.name]==i]
      vol[vol[,id.name]==i,'rvCH4']<- vol[vol[,id.name]==i,'vCH4']/dt[vol[,id.name]==i]
    }

    # Drop t0 if not requested (whether originally present or added)
    if(!showt0) {
      vol<-vol[vol[,time.name]!=0,]
    }

    if(cmethod=='total') {
      vol$vhsCH4<-vol[,comp.name]*vol[,vol.hs.name]
      vol$vtCH4<-vol$vCH4 + vol$vhsCH4
    } 

    # Sort and return results
    vol<-vol[order(vol[,id.name],vol[,time.name]),]
    if(is.null(comp)) {
      vol<-vol[,!names(vol) %in% c(comp.name,'vCH4','cvCH4','rvCH4')]
    }
    rownames(vol)<-1:nrow(vol)
    return(vol)

  } else if(dat.type=='mass') {
    # Gravimetric
    # Work with mass
    message('Working with mass data (applying gravimetric approach).')
    mass<-dat
    mass.name<-dat.name
    # Column names
    if(any(!c(id.name,time.name,mass.name) %in% names(mass))){
      message('Column identities inferred from order for mass: id = mass$',names(mass)[1],', time = mass$',names(mass)[2],', and mass = mass$',names(mass)[3],'.')
      names(mass)<-c(id.name,time.name,mass.name)
    } 

    # Calculate mass loss
    mass<-mass[order(mass[,id.name],mass[,time.name]),]
    # starts data frame is binary, used to track first observation for each reactor, considered the start
    starts<-mass[,c(id.name,time.name)]
    starts$start<-FALSE
    for(i in unique(mass[,id.name])) {
      mass[mass[,id.name]==i,'massloss']<-c(0,-diff(mass[mass[,id.name]==i,mass.name]))
      starts[starts[,id.name]==i,'start'][1]<-TRUE
    }

    # Calculate biogas production
    mass[,c('vBg','vCH4')]<-mass2vol(mass=mass$massloss,xCH4=mass[,comp.name],temp=temp,pres=pres,value='all')[,c('vBg','vCH4')]
    # Set time zero volumes to zero--necessary because xCH4 is always missing
    mass[mass$massloss==0,c('vBg','vCH4')]<-0

    # Cumulative gas production and rates
    mass<-mass[order(mass[,id.name],mass[,time.name]),]
    # Calculate delta t for rates
    if(class(mass[,time.name])[1] %in% c('numeric','integer')) {
      dt<-c(NA,diff(mass[,time.name]))
    } else if(class(mass[,time.name])[1] %in% c('POSIXct','POSIXlt')) {
      dt<-c(NA,as.numeric(diff(mass[,time.name]),units='days'))
    } else {
      dt<-NA
      warning('time column in mass data frame not recognized, so rates will not be calculated.')
    }
    # Set dt to NA for the first observation for each reactor
    dt[c(TRUE,mass[,id.name][-1] != mass[,id.name][-nrow(mass)])]<-NA
    for(i in unique(mass[,id.name])) {
      mass[mass[,id.name]==i,'cvBg']<- cumsum(mass[mass[,id.name]==i,'vBg' ])
      mass[mass[,id.name]==i,'cvCH4']<-cumsum(mass[mass[,id.name]==i,'vCH4'])
      mass[mass[,id.name]==i,'rvBg']<- mass[mass[,id.name]==i,'vBg' ]/dt[mass[,id.name]==i]
      mass[mass[,id.name]==i,'rvCH4']<-mass[mass[,id.name]==i,'vCH4']/dt[mass[,id.name]==i]
    }

    # Drop time 0 or initial times, works even if time column not recognized
    if(!showt0) {
      mass<-mass[!starts$start,]
    }

    # Sort and return results
    mass<-mass[order(mass[,id.name],mass[,time.name]),]
    if(is.null(comp)) {
      vol<-vol[,!names(vol) %in% c(comp.name,'vCH4','cvCH4','rvCH4')]
    }
    rownames(mass)<-1:nrow(mass)
    return(mass)

  } else stop('Either \"vol\" or \"mass\" must be provided.')
}
