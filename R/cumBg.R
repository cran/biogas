# Modified: 23 May 2016 SDH

cumBg <- function(
  # Main arguments
  dat,
  dat.type = 'vol',
  comp = NULL,
  temp = NULL,
  pres = NULL,
  interval = TRUE,
  # Column names for volumetric method
  id.name = 'id',
  time.name = 'time',
  dat.name = dat.type,
  comp.name = 'xCH4', 
  # Additional arguments for manometric and gravimetric methods
  pres.resid = NULL,
  temp.init = NULL,
  pres.init = NULL,
  rh.resid.init = 1,
  headspace = NULL,
  vol.hs.name = 'vol.hs',
  headcomp = 'N2',
  absolute = TRUE,
  # Calculation method and other settings
  cmethod = 'removed',
  imethod = 'linear',
  extrap = FALSE,
  addt0 = TRUE,
  showt0 = TRUE,
  # Warnings and messages
  std.message = TRUE,
  check = TRUE,
  # Units and standard conditions
  temp.std = getOption('temp.std', 0.0), 
  pres.std = getOption('pres.std', 1.0),
  unit.temp = getOption('unit.temp', 'C'),
  unit.pres = getOption('unit.pres', 'atm')
){
  # Check arguments
  checkArgClassValue(dat, 'data.frame')
  checkArgClassValue(dat.type, 'character', expected.values = c('vol', 'mass', 'pres', 'volume', 'pressure'))
  checkArgClassValue(comp, c('data.frame', 'integer', 'numeric', 'NULL'))
  checkArgClassValue(temp, c('integer', 'numeric', 'NULL'))
  checkArgClassValue(pres, c('integer', 'numeric', 'NULL'))
  checkArgClassValue(id.name, 'character')
  checkArgClassValue(time.name, 'character')
  checkArgClassValue(dat.name, 'character')
  checkArgClassValue(comp.name, c('character', 'NULL'))
  checkArgClassValue(headspace, c('data.frame', 'integer', 'numeric', 'NULL'))
  checkArgClassValue(vol.hs.name, 'character')
  checkArgClassValue(headcomp, 'character')
  checkArgClassValue(temp.init, c('integer', 'numeric', 'NULL'))
  checkArgClassValue(temp.std, c('integer', 'numeric'))
  checkArgClassValue(pres.std, c('integer', 'numeric'))
  checkArgClassValue(pres.resid, c('integer', 'numeric', 'character', 'NULL'))
  checkArgClassValue(pres.init, c('integer', 'numeric', 'NULL'))
  checkArgClassValue(rh.resid.init, c('integer', 'numeric', 'NULL'), expected.range = c(0, 1))
  checkArgClassValue(unit.temp, 'character')
  checkArgClassValue(unit.pres, 'character')
  checkArgClassValue(cmethod, 'character', expected.values = c('removed', 'total'))
  # Skip imethod, checked in interp
  checkArgClassValue(extrap, 'logical')
  checkArgClassValue(addt0, 'logical')
  checkArgClassValue(showt0, 'logical')
  checkArgClassValue(std.message, 'logical')
  checkArgClassValue(check, 'logical')
  checkArgClassValue(absolute, 'logical')
  checkArgClassValue(interval, 'logical')

  # Hard-wire rh for now at least
  rh <- 1

  # Check column names in argument data frames
  # comp needs id (time) xCH4, time optional
  if(!is.null(comp) && class(comp)=='data.frame') {
    if(any(missing.col <- !c(id.name, comp.name) %in% names(comp))){
      stop('Specified column(s) in comp data frame (', deparse(substitute(comp)), ') not found: ', c(id.name, comp.name)[missing.col], '.')
    }
  }

  # dat (volume or mass)
  if(any(missing.col <- !c(id.name, time.name, dat.name) %in% names(dat))){
    stop('Specified columns in dat data frame (', deparse(substitute(dat)), ') not found: ', c(id.name, time.name, dat.name)[missing.col], '.')
  } 

  # Check for headspace argument if it is needed
  if(is.null(headspace) & cmethod=='total') stop('cmethod is set to \"total\" but headspace argument is not provided.')

  # Check for necessary arguments
  if(dat.type %in% c('pres', 'pressure')) {
    if(is.null(headspace)) stop('dat.type is set to \"pres\" but \"headspace\" is not provided.')
    if(is.null(temp.init)) stop('dat.type is set to \"pres\" but \"temp.init\" is not provided.')
    if(is.null(pres.resid)) stop('dat.type is set to \"pres\" but \"pres.resid\" is not provided.')
    if(is.null(pres.init)) stop('dat.type is set to \"pres\" but \"pres.init\" is not provided.')
  }

  # Check for input errors
  # NTS: Add other checks here
  if(check) {
    # Composition
    if(is.numeric(comp) | is.integer(comp)) {
      if(any(comp < 0 | comp > 1)) {
	warning('Check biogas composition in ', deparse(substitute(comp)), '. One or more values is outside of range 0.0-1.0.')
      }
    } else {
      if(any(comp[, comp.name] < 0 | comp[, comp.name] > 1)) {
	warning('Check biogas composition in ', deparse(substitute(comp)), '$', comp.name, '. One or more values is outside of range 0.0-1.0.')
      }
    }
  }

  # Convert date.type to lowercase so it is more flexible for users
  dat.type <- tolower(dat.type)

  # Sort out composition data
  mssg.no.time <- mssg.interp <- FALSE
  # First sort so can find first observation for mass data to ignore it
  dat <- dat[order(dat[, id.name], dat[, time.name]), ]
  dat[, comp.name] <- NA
  if(!is.null(comp) && class(comp)=='data.frame'){
    # Interpolate gas composition to times of volume measurements
    for(i in unique(dat[, id.name])) {
      if(dat.type=='mass' & nrow(dat[dat[, id.name]==i, ])<2) stop('There are < 2 observations for reactor ', i,' but dat.type = \"mass\". 
                                                                   You need at least 2 observations to apply the gravimetric method.')
      dc <- comp[comp[, id.name]==i, ]
      if(nrow(dc)==0) stop('No biogas composition data for reactor ', i,' so can\'t interpolate!') 
      if(nrow(dc)>1) {
	# If there is no time column
	if(!time.name %in% names(comp)) stop('Problem with comp  (', deparse(substitute(comp)), 
                                             '): a time column was not found but there is > 1 observation at least for reactor ',i, '.')
	if(dat.type %in% c('vol', 'volume', 'pres', 'pressure')) {
	  mssg.interp <- TRUE
          dat[dat[, id.name]==i, comp.name] <- interp(dc[, time.name], dc[, comp.name], time.out = dat[dat[, id.name]==i, time.name], method = imethod, extrap = extrap)
          # Set first value to zero if there is no biogas production (fixes problem with cmethod = total when there is a t0 observation included)
          dat[dat[, id.name]==i & dat[, dat.name]==0, comp.name][1] <- 0
	} else if (dat.type=='mass') {
	  # Then ignore first point, since it isn't used anyway--this is just to avoid warning with interp if extrap = FALSE
	  mssg.interp <- TRUE
          dat[dat[, id.name]==i, comp.name][-1] <- interp(dc[, time.name], dc[, comp.name], time.out = dat[dat[, id.name]==i, time.name][-1], method = imethod, extrap = extrap)
	}
      } else { # If only one xCH4 value is available, use it for all dat obs if extrap = TRUE or times match, but warn if times don't match
	if(!time.name %in% names(comp)) {
	  # If there is no time column in comp
	  mssg.no.time <- TRUE
          dat[dat[, id.name]==i, comp.name] <- dc[, comp.name]
	} else {
	  # There is a time column in dc/comp
          for(j in 1:nrow(dat[dat[, id.name]==i, ])) {
	    if(j > 1 | dat.type=='vol') { # This just to avoid warning for first observation for mass data
              if(dc[, time.name]==dat[dat[, id.name]==i, time.name][j]) { 
	        # If times match
                dat[dat[, id.name]==i, comp.name][j] <- dc[, comp.name]
              } else {
                if(extrap) {
                  dat[dat[, id.name]==i, comp.name][j] <- dc[, comp.name]
                } else {
                  dat[dat[, id.name]==i, comp.name][j] <- NA
      	          warning('Not enough xCH4 data (one observation) to interpolate for reactor ', i,' so results will be missing.\n If you prefer, you can use extrapolation by setting extrap = TRUE.')
                }
              }
	    }
          }
	}
      }
    }
  } else if(!is.null(comp) && class(comp) %in% c('numeric', 'integer') && length(comp)==1) {
    # Or if a single value is given, use it
    message('Only a single value was provided for biogas composition (', comp, '), so applying it to all observations.')
    dat[, comp.name] <- comp
  } else {
    # If no composition data is given, just use NA
    dat[, comp.name] <- NA 
  }
  if(mssg.no.time) message('A time column was not found in comp (', deparse(substitute(comp)), '), and a single value was used for each reactor.')
  if(mssg.interp) message('Biogas composition is interpolated.')

  # Add headspace if provided
  if(!is.null(headspace)) {
    if(is.numeric(headspace)) {
     dat[, vol.hs.name] <- headspace
    } else if(is.data.frame(headspace)) {       
      # headspace needs id vol
      if(any(missing.col <- !c(id.name, vol.hs.name) %in% names(headspace))){
        stop('Columns with names matching id.name or vol.hs.name are not present in headspace data frame: ', c(id.name, vol.hs.name)[missing.col], '.')
      }
      dat <- merge(dat, headspace[ , c(id.name, vol.hs.name)], by = id.name, suffixes = c('', '.hs'))
    } else stop('headspace actual argument not recognized. What is it?')
  }

  # Volumetric
  # Function will work with vol and add columns
  if(dat.type %in% c('vol', 'volume', 'pres', 'pressure')) {
    # vol dat needs id time vol

    # Standardize total gas volumes
    # Note that temperature and pressure units are not converted at all in cumBg (but are in stdVol of course)
    if(dat.type %in% c('vol', 'volume')) {
      message('Working with volume data, applying volumetric method.')
      if(!is.null(temp) & !is.null(pres)) {
        dat$vBg <- stdVol(dat[, dat.name], temp = temp, pres = pres, rh = rh, pres.std = pres.std, 
                          temp.std = temp.std, unit.temp = unit.temp, unit.pres = unit.pres, 
                          std.message = std.message)
      } else {
        dat$vBg <- dat[, dat.name]
        message('Either temperature or presure is missing (temp and pres arguments) so volumes are NOT standardized.')
      }
    } 

    if(dat.type %in% c('pres', 'pressure')) {
       message('Working with pressure data, pressure measurements are', if (absolute) ' ABSOLUTE' else ' GAUGE', 
              ' If this is incorrect, change \'absolute\' argument to ', !absolute, '.')

      # Add pres.resid to dat if it isn't already present
      if(is.numeric(pres.resid) | is.integer(pres.resid)) {
        dat$pres.resid <- pres.resid
        pres.resid <- 'pres.resid'
      }

      # If absolute != TRUE, calculate absolute pressure and add to dat
      if(!absolute) {
        dat[, aa <- paste0(dat.name, '.abs')] <- dat[, dat.name] + dat[, pres.resid]
        dat.name <- aa
      }

      # Add residual rh (after pressure measurement and venting)
      dat$rh.resid <- dat[, pres.resid]/dat[, dat.name]

      # Sort to add *previous* residual pressure, rh, and temperature columns
      dat <- dat[order(dat[, id.name], dat[, time.name]), ]

      # Add previous residual pressure, residual rh, and previous temperature columns
      for(i in unique(dat[, id.name])) {
        pr <- dat[dat[, id.name]==i, pres.resid]
        rhr <- dat[dat[, id.name]==i, 'rh.resid']
        dat[dat[, id.name]==i, 'pres.resid.prev'] <- c(pres.init, pr[-length(pr)])
        dat[dat[, id.name]==i, 'rh.resid.prev'] <- c(rh.resid.init, rhr[-length(rhr)])
        dat[dat[, id.name]==i, 'temp.prev'] <- c(temp.init, rep(temp, length(pr) - 1))
      }

      # Finally, calculate volume of gas in bottle headspace
      if(interval) {
        # Second call (subtracted bit) is standardized volume in bottle headspace at end of previous measurement (after venting)
        # Result then may not be exactly volume vented in current measurement, but total new gas volume since last measurement (only differ if pres.resid differs)
        dat$vBg <- stdVol(dat[, vol.hs.name], temp = temp, pres = dat[, dat.name], rh = rh, 
                          pres.std = pres.std, temp.std = temp.std, unit.temp = unit.temp, 
                          unit.pres = unit.pres, std.message = FALSE, warn = FALSE) - 
                   stdVol(dat[, vol.hs.name], temp = dat$temp.prev, pres = dat$pres.resid.prev, rh = dat$rh.resid.prev,  
                          pres.std = pres.std, temp.std = temp.std, unit.temp = unit.temp, 
                          unit.pres = unit.pres, std.message = std.message)
      } else {
        # Second call (subtracted bit) is original bottle headspace (standardized), assumed to start at first pres.resid (NTS: may have been sealed at RT--could be changed)
        dat$vBg <- stdVol(dat[, vol.hs.name], temp = temp, pres = dat[, dat.name], rh = rh, 
                          pres.std = pres.std, temp.std = temp.std, unit.temp = unit.temp, 
                          unit.pres = unit.pres, std.message = FALSE, warn = FALSE) - 
                   stdVol(dat[, vol.hs.name], temp = temp, pres = pres.init, rh = rh.resid.init,  
                          pres.std = pres.std, temp.std = temp.std, unit.temp = unit.temp, 
                          unit.pres = unit.pres, std.message = std.message, warn = FALSE)
      }
    }

    # Calculate interval (or cum if interval = FALSE) gas production
    dat$vCH4 <- dat$vBg*dat[, comp.name]*vol.mol['CH4'] / (dat[, comp.name]*vol.mol['CH4'] +(1 - dat[, comp.name])*vol.mol['CO2'])  # CH4 and CO2 molar volumes in ml/mol
    ##if(addt0) dat[dat[, time.name]==0, 'vCH4'] <- 0 # Not clear why this was ever needed

    # For cmethod = 'total', add headspace CH4 to vCH4
    if(cmethod=='total') {
      if(dat.type %in% c('vol', 'volume')) {
        warning('For cmethod = \"total\", headspace temperature is taken as temp (', temp, unit.temp, '), pressure as \"pres\" (', pres, unit.pres, '), and relative humidity as 1.0 (100%).')

        # Small problem with rh assumption here. Will actually be < 1 after gas removal
        # Also assume vol meas pressure pres = residual headspace pressure
        dat$vhsCH4 <- dat[, comp.name]*
                      stdVol(dat[, vol.hs.name], temp = temp, pres = pres, rh = 1, pres.std = pres.std, 
                            temp.std = temp.std, unit.temp = unit.temp, unit.pres = unit.pres, 
                            std.message = std.message)
      }

      if(dat.type %in% c('pres', 'pressure')) {
        dat$vhsCH4 <- dat[, comp.name]*
                      stdVol(dat[, vol.hs.name], temp = temp, pres = dat[ , pres.resid], rh = dat$rh.resid,  
                             pres.std = pres.std, temp.std = temp.std, unit.temp = unit.temp, 
                             unit.pres = unit.pres, std.message = std.message)
      }

      # vhsCH4 is added to cvCH4 below
      # Calculations are up here to avoid t0 issues
    } 

    # Add t0 row if requested
    # Not added if there are already zeroes present!
    if(addt0 & !class(dat[, time.name])[1] %in% c('numeric', 'integer')) addt0 <- FALSE
    if(addt0 & !any(dat[, time.name]==0)) {
      t0 <- data.frame(id = unique(dat[, id.name]), time = 0)
      names(t0) <- c(id.name, time.name)
      t0[, 'vBg'] <- t0[, 'vCH4'] <- 0

      # This messy, but needed for calculating vCH4 by diff when this method is used and interval = FALSE
      if(cmethod == 'total') {
        t0[, 'vhsCH4'] <- 0
      }
      dat <- rbindf(dat, t0)
    }

    # Calculate delta t for rates
    dat <- dat[order(dat[, id.name], dat[, time.name]), ]
    if(class(dat[, time.name])[1] %in% c('numeric', 'integer')) {
      dt <- c(NA, diff(dat[, time.name]))
    } else if(class(dat[, time.name])[1] %in% c('POSIXct', 'POSIXlt')) {
      dt <- c(NA, as.numeric(diff(dat[, time.name]), units = 'days'))
    } else {
      dt <- NA
      warning('time column in dat data frame not recognized, so rates will not be calculated.')
    }
    # Set dt to NA for first observations for each reactor
    dt[c(TRUE, dat[, id.name][-1] != dat[, id.name][-nrow(dat)])] <- NA 

    # May already have cumulative production, if so move it to cvCH4, and calculate vCH4 down below
    if(!interval) {
      dat$cvBg <- dat$vBg
      dat$cvCH4 <- dat$vCH4
    }
    
    # Calculate cumulative production or interval production (depending on interval argument)
    # And calculate rates
    if(interval) {
      for(i in unique(dat[, id.name])) {
        dat[dat[, id.name]==i, 'cvBg'] <- cumsum(dat[dat[, id.name]==i, 'vBg' ])
        dat[dat[, id.name]==i, 'cvCH4'] <- cumsum(dat[dat[, id.name]==i, 'vCH4'])
      } 
    }

    # For cmethod = 'total', add headspace CH4 to cvCH4
    if(cmethod == 'total') {
      dat$cvCH4 <- dat$cvCH4 + dat$vhsCH4
    }

    # For cumulative results, calculate interval production from cvCH4 (down here because it may have headspace CH4 added if cmethod = total) so cannot be combined with cvCH4 calcs above
    # vBg could be moved up, but that means more code
    if(!interval) {
      for(i in unique(dat[, id.name])) {
        dat[dat[, id.name]==i, 'vBg'] <- diff(c(0, dat[dat[, id.name]==i, 'cvBg' ]))
        dat[dat[, id.name]==i, 'vCH4'] <- diff(c(0, dat[dat[, id.name]==i, 'cvCH4']))
      }
    }

    # When cmethod = 'total', cvCH4 must be (re)calculated from cvCH4, because vhsCH4 is added to cvCH4 (correctly)
    # vBg is not affected by cmethod = 'total'
    if(cmethod == 'total') {
      for(i in unique(dat[, id.name])) {
        dat[dat[, id.name]==i, 'vCH4'] <- diff(c(0, dat[dat[, id.name]==i, 'cvCH4']))
      }
    }


    # Calculate rates for all cases 
    for(i in unique(dat[, id.name])) {
      dat[dat[, id.name]==i, 'rvBg'] <- dat[dat[, id.name]==i, 'vBg' ]/dt[dat[, id.name]==i]
      dat[dat[, id.name]==i, 'rvCH4']<- dat[dat[, id.name]==i, 'vCH4' ]/dt[dat[, id.name]==i]
    }

    # Drop t0 if not requested (whether originally present or added)
    if(!showt0) {
      dat <- dat[dat[, time.name] != 0, ]
    }

    # Sort and return results
    dat <- dat[order(dat[, id.name], dat[, time.name]), ]
    if(is.null(comp)) {
      dat <- dat[, ! names(dat) %in% c(comp.name, 'vCH4', 'cvCH4', 'rvCH4')]
    }
    if(all(is.na(dt))) {
      dat <- dat[, ! names(dat) %in% c('rvBg','rvCH4')]
    }
    rownames(dat) <- 1:nrow(dat)

    return(dat)

  } else if(dat.type=='mass') {
    # Gravimetric
    # Work with mass
    message('Working with mass data (applying gravimetric approach).')

    # Check for pressure and temperature--required, but default is NULL (for use of volumetric method without standardization) so must check here 
    if(is.null(temp)) stop('temp argument missing but is required for gravimetric method.')
    if(is.null(pres)) stop('pres argument missing but is required for gravimetric method.')
    if(is.null(comp)) stop('comp argument missing but is required for gravimetric method.')

    # In this section main data frame is saved to `dat`, and name of response (mass) to `mass.name`
    mass <- dat
    mass.name <- dat.name

    # Calculate mass loss
    mass <- mass[order(mass[, id.name], mass[, time.name]), ]
    # starts data frame is binary, used to track first observation for each reactor, considered the start
    starts <- mass[, c(id.name, time.name)]
    starts$start <- FALSE
    for(i in unique(mass[, id.name])) {
      mass[mass[, id.name]==i, 'massloss'] <- c(0, -diff(mass[mass[, id.name]==i, mass.name]))
      mass[mass[, id.name]==i, 'cmassloss'] <- cumsum(mass[mass[, id.name]==i, 'massloss'])
      starts[starts[, id.name]==i, 'start'][1] <- TRUE
    }

    # Calculate biogas production
    if(any(mass[, 'massloss'] < 0)) {
      mass[whichones <- which(mass$massloss < 0), 'massloss'] <- NA
      stop('Mass *gain* calculated for one or more observations. See ', paste('id.name column:', mass[whichones, id.name], ' and time.name column:', mass[whichones - 1, time.name], 'to', mass[whichones, time.name], sep = ' ', collapse = ', '), ' in dat data frame. ')
    }

    mass[, c('vBg', 'vCH4')] <- mass2vol(mass = mass[, 'massloss'], xCH4 = mass[, comp.name], temp = temp, pres = pres, temp.std = temp.std, pres.std = pres.std, unit.temp = unit.temp, unit.pres = unit.pres, value = 'all', std.message = std.message)[, c('vBg', 'vCH4')]
    if(!is.null(headspace)) {
      # Apply initial headspace correction only for times 1 and 2 (i.e., one mass loss measurement per reactor)
      which1and2 <- sort(c(which(starts$start), which(starts$start) + 1) )
      mass[which1and2, c('vBg', 'vCH4')] <- mass2vol(mass = mass$massloss[which1and2], xCH4 = mass[which1and2, comp.name], temp = temp, pres = pres, temp.std = temp.std, pres.std = pres.std, unit.temp = unit.temp, unit.pres = unit.pres, value = 'all', headspace = mass[which1and2, vol.hs.name], headcomp = 'N2', temp.init = temp.init, std.message = FALSE)[, c('vBg', 'vCH4')]
    }
    # Set time zero volumes to zero--necessary because xCH4 is always missing
    mass[mass$massloss==0, c('vBg', 'vCH4')] <- 0

    # Cumulative gas production and rates
    mass <- mass[order(mass[, id.name], mass[, time.name]), ]
    # Calculate delta t for rates
    if(class(mass[, time.name])[1] %in% c('numeric', 'integer')) {
      dt <- c(NA, diff(mass[, time.name]))
    } else if(class(mass[, time.name])[1] %in% c('POSIXct', 'POSIXlt')) {
      dt <- c(NA, as.numeric(diff(mass[, time.name]), units = 'days'))
    } else {
      dt <- NA
      warning('time column in mass data frame not recognized, so rates will not be calculated.')
    }
    # Set dt to NA for the first observation for each reactor
    dt[c(TRUE, mass[, id.name][-1] != mass[, id.name][-nrow(mass)])] <- NA
    for(i in unique(mass[, id.name])) {
      mass[mass[, id.name]==i, 'cvBg']<- cumsum(mass[mass[, id.name]==i, 'vBg' ])
      mass[mass[, id.name]==i, 'cvCH4'] <- cumsum(mass[mass[, id.name]==i, 'vCH4'])
      mass[mass[, id.name]==i, 'rvBg']<- mass[mass[, id.name]==i, 'vBg' ]/dt[mass[, id.name]==i]
      mass[mass[, id.name]==i, 'rvCH4'] <- mass[mass[, id.name]==i, 'vCH4']/dt[mass[, id.name]==i]
    }

    # Drop time 0 or initial times, works even if time column not recognized
    if(!showt0) {
      mass <- mass[!starts$start, ]
    }

    # Sort and return results
    mass <- mass[order(mass[, id.name], mass[, time.name]), ]
    if(is.null(comp)) {
      vol <- vol[, !names(vol) %in% c(comp.name, 'vCH4', 'cvCH4', 'rvCH4')]
    }
    rownames(mass) <- 1:nrow(mass)
    
    return(mass)

  } else stop('Either \"vol\" or \"mass\" must be provided.')
}
