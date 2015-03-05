readFormula <-
function(
  form,
  elements=NULL,        # Set of elements returned, all others ignored, e.g., c('C','H','N','O')
  min.elements=NULL     # Minimum set of elements, will return error if these at least are not included
  ) {
  #fc<-as.list(rep(0,length(elements)))

  form.orig<-form
  #form<-toupper(form)

  # Add implied coefficients of 1
  form<-gsub('([a-zA-Z])([A-Z])','\\11\\2',form)
  form<-gsub('([a-zA-Z])([A-Z])','\\11\\2',form) # Repeated for e.g., COOH
  form<-gsub('([a-zA-Z])$','\\11',form)

  # Extract integer coefficients 
  cc<-as.numeric(strsplit(form,'[A-Za-z]+')[[1]][-1])
  names(cc)<-strsplit(form,'[0-9]+')[[1]]

  # Sort out elements to return
  if(is.null(elements)) elements<-unique(names(cc))
  fc<-numeric(length(elements))
  names(fc)<-elements

  # Fill in fc, summing elements of cc if required (if elements are repeated)
  for(i in elements) {
    for(j in 1:length(cc)) {
      if(names(cc)[j]==i) fc[i]<-fc[i] + cc[j]
    }
  }
  
  # Check for minimum set of elements
  if(!is.null(min.elements)) if(any(!min.elements %in% names(fc)) | any(fc[min.elements] == 0)) stop('Minimum elements required are ',min.elements,' (from min.elements argument), but form is ',form.orig,', interpreted as ',form)

  return(fc)
}
