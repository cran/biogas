summBg <-
function(vol,
         setup,
         id.name='id',
	 time.name='time',
         descrip.name='descrip',
         inoc.name=NULL,
         norm.name=NULL,
         inoc.m.name='minoc',
         vol.name='cvCH4',
         imethod='linear',
         extrap=FALSE,
         when=30,
         show.obs=FALSE) 
{
  message('Response variable (volume) is vol$',vol.name,'.')

  if(class(when) %in% c('numeric','integer')) {
    if(any(!c(id.name,time.name,vol.name) %in% names(vol))){
      message('Column identities inferred from order for vol: id = vol$',names(vol)[1],', time = vol$',names(vol)[2],', and ',vol.name,' = vol$',names(vol)[3],'.')
      names(vol)[1:3]<-c(id.name,time.name,vol.name)
    } 
  } else {
    if(any(!c(id.name,vol.name) %in% names(vol))){
      message('Column identities inferred from order for vol: id = vol$',names(vol)[1],', and ',vol.name,' = vol$',names(vol)[2],'.')
      names(vol)[1:2]<-c(id.name,vol.name)
    } 
  }

  if(any(!c(id.name,descrip.name) %in% names(setup))){
    message('Column identities inferred from order for setup: id = setup$',names(setup)[1],', descrip = setup$',names(setup)[2],'.')
    names(setup)[1:2]<-c(id.name,descrip.name)
  } 

  # Interpolate cvCH4 to common time for each reactor
  s1<-data.frame(id=ids<-unique(setup[,id.name]))
  names(s1)[1]<-id.name
  if(class(when) %in% c('numeric','integer')) {
    # Then interpolate
    for(i in ids) {
      dc<-vol[vol[,id.name]==i,]
      s1[s1[,id.name]==i,vol.name]<-interp(dc[,time.name],dc[,vol.name],time.out=when,method=imethod,extrap=extrap)
    }
  } else if(tolower(when)=='end') { # User just wants to use latest values of volume
    # Sort, in order to find latest values
    vol<-vol[order(vol[,id.name],vol[,vol.name]),]
    for(i in ids) {
      dc<-vol[vol[,id.name]==i,]
      s1[s1[,id.name]==i,vol.name]<-dc[nrow(dc),vol.name]
    }
  }
  # Merge to add mass inoculum and VS in substrate
  s1<-merge(setup,s1,by=id.name)

  # Pull out inoculum data
  if(!is.null(inoc.name) && inoc.name %in% setup[,descrip.name]) {
    message('Inoculumn contribution subtracted based on setup$',inoc.m.name,'.') 
    is1<-s1[s1[,descrip.name]==inoc.name,]
    is1$vol.mi<-is1[,vol.name]/is1[,inoc.m.name]
    inocvol<-c(mn=mean(is1$vol.mi),s=sd(is1$vol.mi))

    # Remaining data
    s1<-s1[s1[,descrip.name]!=inoc.name,]
    if(is.factor(s1[,id.name])) s1[,id.name]<-droplevels(s1[,id.name])

    # Calculate and substract inoc contribution
    s1[,vol.name]<-s1[,vol.name] - inocvol['mn']*s1[,inoc.m.name]
    s1[,'sd.inoc']<-inocvol['s']*s1[,inoc.m.name]

  } else {
    message('Inoculumn contribution not subtracted.') 
    s1[,'sd.inoc']<-0
  }

  if(!is.null(norm.name)) {
    message('CH4 production normalized by setup$',norm.name,'.')
    # Normalize remaining vol by norm.name (typically by substrate VS)
    s1[,vol.name]<-s1[,vol.name]/s1[,norm.name]
    s1[,'sd.inoc']<-s1[,'sd.inoc']/s1[,norm.name]
  } else {
    message('No normalization by substrate mass.')
  }

  if(show.obs) {
    s1$sd.inoc<-NULL
    s1<-s1[order(s1[,descrip.name],s1[,id.name]),]
    return(s1)
  }

  # Summarize by description
  s2<-data.frame(descrip=unique(s1[,descrip.name]))
  names(s2)[1]<-descrip.name
  for(i in unique(s1[,descrip.name])){
    dd<-s1[s1[,descrip.name]==i,]
    s2[s2[,descrip.name]==i,'mean']<-signif(mean(na.omit(dd[,vol.name])),5)
    #s2[s2[,descrip.name]==i,'sd']<-signif(sd(na.omit(dd[,vol.name])),5) # Original code with no contribution of inoc to sd
    s2[s2[,descrip.name]==i,'sd']<-signif(sqrt(sd(na.omit(dd[,vol.name]))^2 + mean(dd[,'sd.inoc'])^2),5) # NTS: Probably there is a better way to do this
    s2[s2[,descrip.name]==i,'n']<-sum(!is.na(dd[,vol.name]))  
  }

  return(s2)
}
