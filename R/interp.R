# Modified: 4 MAR 2015

interp <-
function(times,y,time.out,method='linear',extrap=FALSE) {

  # Because spline will return values without y, make sure it is a valid vector
  if(!is.numeric(y) & !is.integer(y)) stop('Problem with y argument, class is: ',class(y))

  # Warn about extrapolation
  #if(extrap & (max(time.out)>max(times) | min(time.out)<min(times))) warning('Some results are extrapolated. Set extrap = FALSE to prevent this.')

  if(method=='linear') {
    y.out<-approx(times,y,xout=time.out,rule=ifelse(extrap,2,1))$y
  } else { 
    # Or for natural, hyman, or fmm methods (others will return an error)
    y.out<-spline(times,y,xout=time.out,method=method)$y
    # If no extrapolation, delete extrapolated values
    if(!extrap) y.out[time.out<min(times) | time.out>max(times)]<-NA
  }

  # Check for missing values
  if(any(is.na(y.out)) & !extrap) warning(sum(is.na(y.out)),' NA(s) in calculated values. If you don\'t mind extrapolation, set extrap=TRUE.')
  # Add names to output
  if(is.numeric(time.out) | is.integer(time.out)) {
    names(y.out)<-paste0('t',signif(time.out,3)) 
  } else {
    names(y.out)<-paste0('t',1:length(time.out)) 
  }

  return(y.out)
}
