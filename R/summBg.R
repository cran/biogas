# Modified: 28 JULY 2015 SDH

summBg <- function(
  vol,
  setup,
  id.name = 'id',
  time.name = 'time',
  descrip.name = 'descrip',
  inoc.name = NULL,
  norm.name = NULL,
  inoc.m.name = 'minoc',
  vol.name = 'cvCH4',
  imethod = 'linear',
  extrap = FALSE,
  when = 30,
  show.obs = FALSE) 
{

  # Argument checks
  checkArgClassValue(vol, 'data.frame')
  checkArgClassValue(setup, 'data.frame')
  checkArgClassValue(id.name, 'character')
  checkArgClassValue(time.name, c('character', 'NULL'))
  checkArgClassValue(descrip.name, c('character', 'NULL'))
  checkArgClassValue(inoc.name, c('character', 'NULL'))
  checkArgClassValue(norm.name, c('character', 'NULL'))
  checkArgClassValue(inoc.m.name, c('character', 'NULL'))
  checkArgClassValue(vol.name, 'character')
  # Skip imethod, checked in interp()
  checkArgClassValue(extrap, 'logical')
  checkArgClassValue(when, c('numeric', 'integer', 'character', 'NULL'))
  checkArgClassValue(show.obs, 'logical')

  # Echo response variable
  message('Response variable (volume) is ', deparse(substitute(vol)), '$', vol.name, '.')

  # Check for missing columns in vol
  if(class(when) %in% c('numeric', 'integer')) {
    if(any(missing.col <- !c(id.name, time.name, vol.name) %in% names(vol))){
      stop('Specified columns in vol data frame (', deparse(substitute(vol)), ') not found: ', c(id.name, time.name, vol.name)[missing.col], '.')
    } 
  } else { # when is 'end' or NULL here
    if(any(missing.col <- !c(id.name, vol.name) %in% names(vol))){
      stop('Specified columns in vol data frame (', deparse(substitute(vol)), ') not found: ', c(id.name, vol.name)[missing.col], '.')
    } 
  }

  # Check for missing columns in setup
  if(any(missing.col <- !c(id.name, descrip.name) %in% names(setup))){
    stop('Specified columns in setup data frame (', deparse(substitute(setup)), ') not found: ', c(id.name, descrip.name)[missing.col], '.')
  } 

  # Check that inoc.name and norm.name can be found in setup data frame
  if(!is.null(inoc.name) && !inoc.name %in% setup[, descrip.name]) stop('inoc.name ', deparse(substitute(inoc.name)), ' not found in ', deparse(substitute(setup)), '$', descrip.name, '.')
  if(!is.null(norm.name) && !norm.name %in% names(setup)) stop('norm.name ', deparse(substitute(norm.name)), ' not found in the column names of ', deparse(substitute(setup)), '.')

  # NTS: Any other checks to add?

  # Find reactor/bottle IDs present in both vol and setup
  ids <- intersect(setup[, id.name], vol[, id.name])

  # Trim setup and check again for inoc.name and norm.name
  setup <- setup[setup[, id.name] %in% ids, ]
  if(!is.null(inoc.name) && !inoc.name %in% setup[, descrip.name]) stop('inoc.name ', deparse(substitute(inoc.name)), ' no longer in setup after trimming--are reactors present in setup missing in vol?')
  if(!is.null(norm.name) && !norm.name %in% names(setup)) stop('norm.name ', deparse(substitute(norm.name)), 'no longer in setup after trimming--are reactors present in setup missing in vol?')

  # Interpolate cvCH4 to common time for each reactor
  # Next line only for time-independent option with one observation per reactor
  if(is.null(when)) when <- 'one.obs'
  summ1 <- expand.grid(id = ids, time = when)
  names(summ1) <- c(id.name, time.name)

  if(class(when) %in% c('numeric', 'integer')) {
    # Then interpolate
    for(i in ids) {
      dc <- vol[vol[, id.name]==i, ]
      # Interpolate if more than one value is present
      if(nrow(dc)>1) {
        summ1[summ1[, id.name]==i, vol.name] <- interp(dc[, time.name], dc[, vol.name], time.out = when, method = imethod, extrap = extrap)
      } else {
	if(dc[, time.name]==when) { # `when` argument matches the single time present
          summ1[summ1[, id.name]==i, vol.name] <- dc[, vol.name]
	} else {
          summ1[summ1[, id.name]==i, vol.name] <- NA
      	  warning('There is only a single ', vol.name, ' value for reactor ', i,', and it does not match the specified when (', when, '). Interpolation is not possible.')
	}
      }
    }
  } else if(length(when) == 1 && tolower(when) == 'end') { # User just wants to use latest values of volume
    # Sort, in order to find latest values
    vol <- vol[order(vol[, id.name], vol[, vol.name]), ]
    for(i in ids) {
      dc <- vol[vol[, id.name]==i, ]
      # Select the last row from sorted data frame
      summ1[summ1[, id.name]==i, vol.name] <- dc[nrow(dc), vol.name]
    }
  } else if(length(when) == 1 && when == 'one.obs') { # Only one value provided per reactor, use it
    for(i in ids) {
      dc <- vol[vol[, id.name]==i, ]
      if(nrow(dc) > 1) stop('Argument when is NULL, which implies one observation per reactor, but for reactor ', i, ' at least, there are multiple observations.')
      summ1[summ1[, id.name]==i, vol.name] <- dc[, vol.name]
    }
  } else stop('when argument not recognized. Options are numeric or integer vector, \"end\".')

  # Merge to add mass inoculum and VS in substrate
  summ1 <- merge(setup, summ1, by = id.name)

  # Messages about inoculum and normalization
  if(!is.null(inoc.name) && inoc.name %in% setup[, descrip.name]) { # Inoculum contribution subtracted
    #message('Inoculum contribution subtracted based on ', deparse(substitute(setup.orig)), '$', inoc.m.name, '.') 
    message('Inoculum contribution subtracted based on setup$', inoc.m.name, '.') 
  } else {
      message('Inoculum contribution not subtracted.') 
  }

  # Message about normalization
  if(!is.null(norm.name)) { 
    #message('Response normalized by ', deparse(substitute(setup)), '$', norm.name, '.')
    message('Response normalized by setup$', norm.name, '.')
  } else {
    message('No normalization by substrate mass.')
  }

  # Loop through all times within when argument
  # NTS: there may be a more efficient way to do this
  summ2 <- data.frame()
  for (w in when) {
    s1 <- summ1[summ1[, time.name]==w, ]

    # Pull out inoculum data
    if(!is.null(inoc.name) && inoc.name %in% setup[, descrip.name]) { # Inoculum contribution subtracted
      # Subset with just inoculum data
      is1 <- s1[s1[, descrip.name]==inoc.name, ]

      # Volume contribution per unit inoculum mass
      is1$vol.mi <- is1[, vol.name]/is1[, inoc.m.name]

      # Mean and sd volume contribution per unit inoc mass
      inocvol <- c(mn = mean(is1$vol.mi), s = sd(is1$vol.mi))

      # Remaining data (drop out inoculum data)
      s1 <- s1[s1[, descrip.name]!=inoc.name, ]
      if(is.factor(s1[, id.name])) s1[, id.name] <- droplevels(s1[, id.name])

      # Calculate and substract inoc contribution
      # Next three lines only for returning additional info when show.obs = TRUE
      s1[, paste0(vol.name, '.tot')] <- s1[, vol.name]
      s1[, paste0(vol.name, '.inoc')] <- inocvol['mn']*s1[, inoc.m.name]
      s1[, 'fv.inoc'] <- s1[, paste0(vol.name, '.inoc')]/s1[, paste0(vol.name, '.tot')]

      # Correct vol for inoculum
      s1[, vol.name] <- s1[, vol.name] - inocvol['mn']*s1[, inoc.m.name]

      # Add sd in volume produced by inoculum for use below in error propagation
      s1[, 'sd.inoc'] <- inocvol['s']*s1[, inoc.m.name]

    } else {
      s1[, 'sd.inoc'] <- 0
    }

    if(!is.null(norm.name)) { 
      # Normalize remaining vol by norm.name (typically by substrate VS)
      s1[, vol.name] <- s1[, vol.name]/s1[, norm.name]

      # Normalize sd contribution from inoc by the same value
      s1[, 'sd.inoc'] <- s1[, 'sd.inoc']/s1[, norm.name]

      # Next two lines only for returning additional info when show.obs = TRUE
      # Only have the .tot and .inoc columns when inoc is subtracted out
      if(!is.null(inoc.name) && inoc.name %in% setup[, descrip.name]) { 
        s1[, paste0(vol.name, '.tot')] <- s1[, paste0(vol.name, '.tot')]/s1[, norm.name]
        s1[, paste0(vol.name, '.inoc')] <- s1[, paste0(vol.name, '.inoc')]/s1[, norm.name]
      }
    } 

    # Calculate means and sd for a summary
    if(!show.obs) {
      # Summarize by description
      s2 <- data.frame(descrip = unique(s1[, descrip.name]), time = w)
      names(s2) <- c(descrip.name, time.name)
      for(i in unique(s1[, descrip.name])){
        dd <- s1[s1[, descrip.name]==i, ]
        s2[s2[, descrip.name]==i, 'mean'] <- signif(mean(na.omit(dd[, vol.name])), 5)
        #s2[s2[, descrip.name]==i, 'sd'] <- signif(sd(na.omit(dd[, vol.name])), 5) # Original code with no contribution of inoc to sd
        s2[s2[, descrip.name]==i, 'sd'] <- signif(sqrt(sd(na.omit(dd[, vol.name]))^2 + mean(dd[, 'sd.inoc'])^2), 5) # NTS: Probably there is a better way to do this
        s2[s2[, descrip.name]==i, 'n'] <- sum(!is.na(dd[, vol.name]))  
      }
    } else { # If show.obs = TRUE, just return individual observations
      #s1$sd.inoc <- NULL
      s2 <- s1[order(s1[, descrip.name], s1[, id.name]), ]
    }

    summ2 <- rbind(summ2, s2)
  }

  # If there is no time, drop time column
  if(length(when) == 1 && when == 'one.obs') summ2[, time.name] <- NULL

  return(summ2)
}
