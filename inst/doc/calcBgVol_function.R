## ----include=FALSE, cache=FALSE--------------------------------
library(knitr)
#opts_chunk$set(cache=FALSE,tidy=FALSE,highlight=FALSE)
opts_chunk$set(cache = FALSE, tidy = FALSE, fig.align = "center")

library(biogas)
options(width=65)

## --------------------------------------------------------------
data("s3lcombo")

dim(s3lcombo)

s3lcombo

summary(s3lcombo)

## --------------------------------------------------------------
cum.prod.lc <- calcBgVol(s3lcombo, temp = 25, pres = 1, 
                         time.name = 'time.d', vol.name = 'vol.ml', 
                         comp.name = 'xCH4',
                         extrap = TRUE)

## --------------------------------------------------------------
head(cum.prod.lc)

dim(cum.prod.lc)

## ----fig.width=6, fig.height=4, fig.align="center"-------------
library(ggplot2)

ggplot(cum.prod.lc, aes(time.d, cvCH4, colour = factor(id))) + 
  geom_point() +
  geom_line(aes(group = id)) +
  labs(x = "Time [d]", y = "Cumulative methane production  [mL]", colour = "Bottle ID")  + 
  theme_bw() 

## --------------------------------------------------------------
data("feedVol")

dim(feedVol)

head(feedVol)

summary(feedVol)

## --------------------------------------------------------------
cum.prod.w <- calcBgVol(feedVol, comp = 1, temp = 0, pres = 1,
                       data.struct = "wide",
                       time.name = "time.d", vol.name = "1",
                       dry = TRUE,
                       interval = FALSE)

## --------------------------------------------------------------
head(cum.prod.w)

dim(cum.prod.w)

## ----fig.width=6, fig.height=4, fig.align="center"-------------
ggplot(cum.prod.w, aes(time.d, cvCH4, colour = factor(id))) + 
       geom_point() +
       geom_line(aes(group = id)) +
       labs(x = "Time [d]", y = "Cumulative methane production  [mL]", colour = "Bottle ID")  + 
       theme_bw() 

## --------------------------------------------------------------
data("vol")

dim(vol)

head(vol)

summary(vol)

## --------------------------------------------------------------
data("comp")

dim(comp)

head(comp)

summary(comp)

## --------------------------------------------------------------
cum.prod.l <- calcBgVol(vol, comp = comp, temp = 35, pres = 1,
                       data.struct = "long",
                       time.name = "days", vol.name = "vol", comp.name = "xCH4", 
                       extrap = TRUE)

## --------------------------------------------------------------
head(cum.prod.l)

dim(cum.prod.l)

## ----fig.width=6, fig.height=4, fig.align="center"-------------
ggplot(cum.prod.l, aes(days, cvCH4, colour = factor(id))) + 
      geom_point() +
      geom_line(aes(group = id)) +
      labs(x = "Time [d]", y = "Cumulative methane production  [mL]", 
           colour = "Bottle ID")  + 
      theme_bw() 

## --------------------------------------------------------------
data("vol")

dim(vol)

head(vol)

summary(vol)

## --------------------------------------------------------------
data("comp")

dim(comp)

head(comp)

summary(comp)

## --------------------------------------------------------------
vol$temp <- 35
vol$pres <- NA
vol$pres <- rnorm(vol$pres, mean = 1, sd = 0.001)

head(vol)

## --------------------------------------------------------------
cum.prod <- calcBgVol(vol, comp = comp, temp = "temp", pres = "pres",
                       data.struct = "long",
                       time.name = "days", vol.name = "vol", comp.name = "xCH4", 
                       extrap = TRUE)

## --------------------------------------------------------------
cum.prod <- calcBgVol(vol, comp = comp, temp = "temp", pres = "pres",
                       data.struct = "long",
                       time.name = "days", vol.name = "vol", comp.name = "xCH4", 
                       extrap = TRUE, showt0 = FALSE)

head(cum.prod)

## ----fig.width=6, fig.height=4, fig.align="center"-------------
ggplot(cum.prod, aes(days, cvCH4, colour = factor(id))) + 
      geom_point() +
      geom_line(aes(group = id)) +
      labs(x = "Time [d]", y = "Cumulative methane production  [mL]", 
           colour = "Bottle ID")  + 
      theme_bw() 

## --------------------------------------------------------------
cum.prod <- calcBgVol(vol, comp = comp, temp = "temp", pres = "pres",
                       data.struct = "long",
                       time.name = "days", vol.name = "vol", comp.name = "xCH4", 
                       extrap = TRUE, addt0 = FALSE)

head(cum.prod)

## ----fig.width=6, fig.height=4, fig.align="center"-------------
ggplot(cum.prod, aes(days, cvCH4, colour = factor(id))) + 
      geom_point() +
      geom_line(aes(group = id)) +
      labs(x = "Time [d]", y = "Cumulative methane production  [mL]", 
           colour = "Bottle ID")  + 
      theme_bw() 

## --------------------------------------------------------------
data("vol")
data("comp")

comp[10,"xCH4"] <- 1.5
cum.prod <- calcBgVol(vol, comp = comp, temp = 35, pres = 1,
                       data.struct = "long",
                       time.name = "days", vol.name = "vol", comp.name = "xCH4", 
                       extrap = TRUE)

## --------------------------------------------------------------
cum.prod <- calcBgVol(vol, comp = comp, temp = 35, pres = 1,
                       data.struct = "long",
                       time.name = "days", vol.name = "vol", comp.name = "xCH4", 
                       extrap = TRUE, check = FALSE)

## --------------------------------------------------------------
vol$temp <- 35 + 273.15
vol$pres <- NA
vol$pres <- rnorm(vol$pres, mean = 101.325, sd = 0.101325)

head(vol)

## --------------------------------------------------------------
cum.prod <- calcBgVol(vol, comp = comp, temp = "temp", pres = "pres",
                       data.struct = "long",
                       time.name = "days", vol.name = "vol", comp.name = "xCH4", 
                       extrap = TRUE)

## --------------------------------------------------------------
cum.prod <- calcBgVol(vol, comp = comp, temp = "temp", pres = "pres",
                       data.struct = "long",
                       time.name = "days", vol.name = "vol", comp.name = "xCH4", 
                       extrap = TRUE,
                       unit.pres = "kPa", unit.temp = "K")

