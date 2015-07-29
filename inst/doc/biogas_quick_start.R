## ----setoptions, include = FALSE, cache = FALSE------------------
  library(knitr) # per http://stackoverflow.com/questions/16752765/knitr-r-package-check-error-object-opts-chunk-not-found 
  opts_chunk$set(cache = FALSE, tidy = FALSE, highlight = FALSE, options(width = 67))
  #library(biogas)

## ----------------------------------------------------------------
library(biogas)
  
data("vol")

dim(vol)

head(vol)

summary(vol)

## ----------------------------------------------------------------
data("comp")

dim(comp)

head(comp)

summary(comp)

## ----------------------------------------------------------------
data("setup")

setup

## ----------------------------------------------------------------
cum.prod <- cumBg(vol, comp = comp, time.name = "days", temp = 35, pres = 1, extrap = TRUE)

## ----------------------------------------------------------------
head(cum.prod)

dim(cum.prod)

## ----fig.width=6, fig.height=4, fig.align="center"---------------
library(ggplot2)
qplot(x = days, y = cvCH4, data = cum.prod, xlab = "Time (d)", 
      ylab = "Cumulative methane production (mL)",color = id, 
      geom = "line")

## ----cprod30-----------------------------------------------------
summBg(cum.prod, setup = setup, time.name = "days", 
       descrip.name = "descrip", when = 30)

## ----------------------------------------------------------------
BMP <- summBg(cum.prod, setup = setup, time.name = "days", 
	      inoc.name = "inoc", inoc.m.name = "minoc", 
	      norm.name = "mvs.sub", when = 60)
BMP

## ----------------------------------------------------------------
predBg("C6H10O5")

## ----cellbmpcomp-------------------------------------------------
BMP

## ----------------------------------------------------------------
predBg("C6H10O5", fs = 0.1)

## ----------------------------------------------------------------
predBg(COD = c(A = 1.439, B = 1.561, cellu = 1.184))

## ----------------------------------------------------------------
  BMP$mean/predBg(COD = c(A = 1.439, B = 1.561, cellu = 1.184))

