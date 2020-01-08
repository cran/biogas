## ----include=FALSE, cache=FALSE--------------------------------
library(knitr)
#opts_chunk$set(cache=FALSE,tidy=FALSE,highlight=FALSE)
opts_chunk$set(cache = FALSE, tidy = FALSE, fig.align = "center")
library(biogas)
options(width=65)

## --------------------------------------------------------------
data("sludgeTwoBiogas")

dim(sludgeTwoBiogas)

head(sludgeTwoBiogas)

summary(sludgeTwoBiogas)

## --------------------------------------------------------------
data("sludgeTwoSetup")

dim(sludgeTwoSetup)

head(sludgeTwoSetup)

summary(sludgeTwoSetup)

## --------------------------------------------------------------
cum.prod.lc <- calcBgMan(sludgeTwoBiogas, temp = 30,
                      time.name = "time.d", comp.name = "xCH4n",
                      temp.init = 30, pres.init = 0.0, 
                      pres.resid = 0,
                      headspace = sludgeTwoSetup,
                      pres.amb = 1013, absolute = FALSE,
                      unit.pres = "mbar")

## --------------------------------------------------------------
head(cum.prod.lc)

dim(cum.prod.lc)

## ----fig.width=6, fig.height=4, fig.align="center"-------------
library(ggplot2)

ggplot(cum.prod.lc, aes(time.d, cvCH4, colour = factor(id))) + 
  geom_point() +
  geom_line(aes(group = id)) +
  labs(x = "Time [d]", y = "cvCH4  [mL]", colour = "Bottle id")  + 
  theme_bw() 

## --------------------------------------------------------------
data("strawPressure")

dim(strawPressure)

head(strawPressure)

summary(strawPressure)

## --------------------------------------------------------------
data("strawComp")

dim(strawComp)

head(strawComp)

summary(strawComp)

## --------------------------------------------------------------
data("strawSetup")

dim(strawSetup)

head(strawSetup)

summary(strawSetup)

## --------------------------------------------------------------
cum.prod.l <- calcBgMan(strawPressure, comp = strawComp, temp = 31,
                        data.struct = "long",
                        time.name = "time", id.name = "bottle", comp.name = "xCH4",
                        temp.init = 21.55, pres.resid = "pres.resid", 
                        pres.init = 0.0,
                        headspace = strawSetup, vol.hs.name = "headspace",
                        pres.amb = 101.3, absolute = FALSE,
                        extrap = TRUE, 
                        unit.pres = "kPa")


head(cum.prod.l)

## ----fig.width=6, fig.height=4, fig.align="center"-------------
ggplot(cum.prod.l, aes(time, cvCH4, colour = factor(bottle))) + 
  geom_point() +
  geom_line(aes(group = bottle)) +
  labs(x = "Time [d]", y = "cvCH4  [mL]", colour = "Bottle id")  + 
  theme_bw() 

