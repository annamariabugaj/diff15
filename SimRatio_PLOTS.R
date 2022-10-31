# Comment: we add the negligible constant of 0.00001 to the denominator as a pre-caution to #ensure that the denominator is non-zero.
install.packages("corrplot")
library(corrplot)

## To calculate the topological overlap for all genes within each module for human and chimp: # Magenta
distTOMHumanmagenta <- TOMdist1(AdjMatHumanrestmagenta)
simTOMHumanmagenta = 1-distTOMHumanmagenta
diag(simTOMHumanmagenta)=0
distTOMChimpmagenta <- TOMdist1(AdjMatChimprestmagenta)
simTOMChimpmagenta = 1-distTOMChimpmagenta
diag(simTOMChimpmagenta)=0
## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# Magenta 
simRatiomagenta=simTOMHumanmagenta/mean(simTOMHumanmagenta)/(simTOMHumanmagenta/mean(simTOMHumanmagenta)+ 
                                                               simTOMChimpmagenta/mean(simTOMChimpmagenta)+0.00001)


## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# Brown
distTOMHumanbrown <- TOMdist1(AdjMatHumanrestbrown)
simTOMHumanbrown = 1-distTOMHumanbrown
diag(simTOMHumanbrown)=0
distTOMChimpbrown <- TOMdist1(AdjMatChimprestbrown)
simTOMChimpbrown = 1-distTOMChimpbrown
diag(simTOMChimpbrown)=0

simRatiobrown=simTOMHumanbrown/mean(simTOMHumanbrown)/(simTOMHumanbrown/mean(simTOMHumanbrown)+ 
                                                         simTOMChimpbrown/mean(simTOMChimpbrown)+0.00001)


###### BROWN ######################
corrplot(simRatiobrown, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_brown <- simRatiobrown > 0.65
### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat_brown <- bigratio_brown
dat_brown <- 1*dat_brown
head(dat_brown)
corrplot(dat_brown, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_brown <- simRatiobrown < 0.35
low_brown <- lowratio_brown
low_brown <- 1*low_brown
head(low_brown)
corrplot(low_brown, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))


## To calculate the topological overlap for all genes within each module for human and chimp: # Black
distTOMHumanblack <- TOMdist1(AdjMatHumanrestblack)
simTOMHumanblack = 1-distTOMHumanblack
diag(simTOMHumanblack)=0
distTOMChimpblack <- TOMdist1(AdjMatChimprestblack)
simTOMChimpblack = 1-distTOMChimpblack
diag(simTOMChimpblack)=0
## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# Black
simRatioblack=simTOMHumanblack/mean(simTOMHumanblack)/(simTOMHumanblack/mean(simTOMHumanblack)+ 
                                                         simTOMChimpblack/mean(simTOMChimpblack)+0.00001)


###### BLACK ######################
corrplot(simRatioblack, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_black <- simRatioblack > 0.65
write.csv(simRatioblack, file = "simRatioblack.csv")
write.csv(bigratio_black, file = "bigRatioblack.csv")
bigratio_black

### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat <- bigratio_black
dat <- 1*dat
head(dat)
corrplot(dat, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_black <- simRatioblack < 0.35
low_black <- lowratio_black
low_black <- 1*low_black
head(low_black)
corrplot(low_black, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))




## To calculate the topological overlap for all genes within each module for human and chimp: # Yellow
distTOMHumanyellow <- TOMdist1(AdjMatHumanrestyellow)
simTOMHumanyellow = 1-distTOMHumanyellow
diag(simTOMHumanyellow)=0
distTOMChimpyellow <- TOMdist1(AdjMatChimprestyellow)
simTOMChimpyellow = 1-distTOMChimpyellow
diag(simTOMChimpyellow)=0
## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# Yellow
simRatioyellow=simTOMHumanyellow/mean(simTOMHumanyellow)/(simTOMHumanyellow/mean(simTOMHumanyellow)+ 
                                                            simTOMChimpyellow/mean(simTOMChimpyellow)+0.00001)



## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# BLUE
distTOMHumanblue <- TOMdist1(AdjMatHumanrestblue)
simTOMHumanblue = 1-distTOMHumanblue
diag(simTOMHumanblue)=0
distTOMChimpblue <- TOMdist1(AdjMatChimprestblue)
simTOMChimpblue = 1-distTOMChimpblue
diag(simTOMChimpblue)=0

simRatioblue=simTOMHumanblue/mean(simTOMHumanblue)/(simTOMHumanblue/mean(simTOMHumanblue)+ 
                                                         simTOMChimpblue/mean(simTOMChimpblue)+0.00001)


###### BLUE ######################
corrplot(simRatioblue, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_blue <- simRatioblue > 0.65
### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat_blue <- bigratio_blue
dat_blue <- 1*dat_blue
head(dat_blue)
corrplot(dat_blue, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_blue <- simRatioblue < 0.35
low_blue <- lowratio_blue
low_blue <- 1*low_blue
head(low_blue)
corrplot(low_blue, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))



## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# GREEN
distTOMHumangreen <- TOMdist1(AdjMatHumanrestgreen)
simTOMHumangreen = 1-distTOMHumangreen
diag(simTOMHumangreen)=0
distTOMChimpgreen <- TOMdist1(AdjMatChimprestgreen)
simTOMChimpgreen = 1-distTOMChimpgreen
diag(simTOMChimpgreen)=0

simRatiogreen=simTOMHumangreen/mean(simTOMHumangreen)/(simTOMHumangreen/mean(simTOMHumangreen)+ 
                                                      simTOMChimpgreen/mean(simTOMChimpgreen)+0.00001)


###### GREEN ######################
corrplot(simRatiogreen, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_green <- simRatiogreen > 0.65
### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat_green <- bigratio_green
dat_green <- 1*dat_green
head(dat_green)
corrplot(dat_green, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_green <- simRatiogreen < 0.35
low_green <- lowratio_green
low_green <- 1*low_green
head(low_green)
corrplot(low_green, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))


## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# PINK
distTOMHumanpink <- TOMdist1(AdjMatHumanrestpink)
simTOMHumanpink = 1-distTOMHumanpink
diag(simTOMHumanpink)=0
distTOMChimppink <- TOMdist1(AdjMatChimprestpink)
simTOMChimppink = 1-distTOMChimppink
diag(simTOMChimppink)=0

simRatiopink=simTOMHumanpink/mean(simTOMHumanpink)/(simTOMHumanpink/mean(simTOMHumanpink)+ 
                                                         simTOMChimppink/mean(simTOMChimppink)+0.00001)


###### PINK ######################
corrplot(simRatiopink, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.4,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_pink <- simRatiopink > 0.65
### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat_pink <- bigratio_pink
dat_pink <- 1*dat_pink
head(dat_pink)
corrplot(dat_pink, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.4,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_pink <- simRatiopink < 0.35
low_pink <- lowratio_pink
low_pink <- 1*low_pink
head(low_pink)
corrplot(low_pink, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.4,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))


## To calculate the topological overlap for all genes within each module for human and chimp: # Turquoise

############ TURQOUISE ################

distTOMHumanturquoise <- TOMdist1(AdjMatHumanrestturquoise)
simTOMHumanturquoise = 1-distTOMHumanturquoise
diag(simTOMHumanturquoise)=0
distTOMChimpturquoise <- TOMdist1(AdjMatChimprestturquoise)
simTOMChimpturquoise = 1-distTOMChimpturquoise
diag(simTOMChimpturquoise)=0
## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# Black
simRatioturquoise=simTOMHumanturquoise/mean(simTOMHumanturquoise)/(simTOMHumanturquoise/mean(simTOMHumanturquoise)+ 
                                                                     simTOMChimpturquoise/mean(simTOMChimpturquoise)+0.00001)
# Comment: we add the negligible constant of 0.00001 to the denominator as a pre-caution to #ensure that the denominator is non-zero.

corrplot(simRatioturquoise, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_turquoise <- simRatioturquoise > 0.65
### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat_turquoise <- bigratio_turquoise
dat_turquoise <- 1*dat_turquoise
head(dat_turquoise)
corrplot(dat_turquoise, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_turquoise <- simRatioturquoise < 0.35
low_turquoise <- lowratio_turquoise
low_turquoise <- 1*low_turquoise
head(low_turquoise)
corrplot(low_turquoise, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.2,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))


## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# GREENYELLOW
distTOMHumangreenyellow <- TOMdist1(AdjMatHumanrestgreenyellow)
simTOMHumangreenyellow = 1-distTOMHumangreenyellow
diag(simTOMHumangreenyellow)=0
distTOMChimpgreenyellow <- TOMdist1(AdjMatChimprestgreenyellow)
simTOMChimpgreenyellow = 1-distTOMChimpgreenyellow
diag(simTOMChimpgreenyellow)=0

simRatiogreenyellow=simTOMHumangreenyellow/mean(simTOMHumangreenyellow)/(simTOMHumangreenyellow/mean(simTOMHumangreenyellow)+ 
                                                         simTOMChimpgreenyellow/mean(simTOMChimpgreenyellow)+0.00001)


###### GREENYELLOW ######################
corrplot(simRatiogreenyellow, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_greenyellow <- simRatiogreenyellow > 0.65
### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat_greenyellow <- bigratio_greenyellow
dat_greenyellow <- 1*dat_greenyellow
head(dat_greenyellow)
corrplot(dat_greenyellow, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_greenyellow <- simRatiogreenyellow < 0.35
low_greenyellow <- lowratio_greenyellow
low_greenyellow <- 1*low_greenyellow
head(low_greenyellow)
corrplot(low_greenyellow, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))


## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# RED
distTOMHumanred <- TOMdist1(AdjMatHumanrestred)
simTOMHumanred = 1-distTOMHumanred
diag(simTOMHumanred)=0
distTOMChimpred <- TOMdist1(AdjMatChimprestred)
simTOMChimpred = 1-distTOMChimpred
diag(simTOMChimpred)=0

simRatiored=simTOMHumanred/mean(simTOMHumanred)/(simTOMHumanred/mean(simTOMHumanred)+ 
                                                      simTOMChimpred/mean(simTOMChimpred)+0.00001)


###### red ######################
corrplot(simRatiored, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.4,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_red <- simRatiored > 0.65
### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat_red <- bigratio_red
dat_red <- 1*dat_red
head(dat_red)
corrplot(dat_red, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.4,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_red <- simRatiored < 0.35
low_red <- lowratio_red
low_red <- 1*low_red
head(low_red)
corrplot(low_red, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.4,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))

## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# PINK
distTOMHumanpink <- TOMdist1(AdjMatHumanrestpink)
simTOMHumanpink = 1-distTOMHumanpink
diag(simTOMHumanpink)=0
distTOMChimppink <- TOMdist1(AdjMatChimprestpink)
simTOMChimppink = 1-distTOMChimppink
diag(simTOMChimppink)=0

simRatiopink=simTOMHumanpink/mean(simTOMHumanpink)/(simTOMHumanpink/mean(simTOMHumanpink)+ 
                                                      simTOMChimppink/mean(simTOMChimppink)+0.00001)


###### PINK ######################
corrplot(simRatiopink, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.4,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_pink <- simRatiopink > 0.65
### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat_pink <- bigratio_pink
dat_pink <- 1*dat_pink
head(dat_pink)
corrplot(dat_pink, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.4,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_pink <- simRatiopink < 0.35
low_pink <- lowratio_pink
low_pink <- 1*low_pink
head(low_pink)
corrplot(low_pink, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.4,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))

## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# PURPLE
distTOMHumanpurple <- TOMdist1(AdjMatHumanrestpurple)
simTOMHumanpurple = 1-distTOMHumanpurple
diag(simTOMHumanpurple)=0
distTOMChimppurple <- TOMdist1(AdjMatChimprestpurple)
simTOMChimppurple = 1-distTOMChimppurple
diag(simTOMChimppurple)=0

simRatiopurple=simTOMHumanpurple/mean(simTOMHumanpurple)/(simTOMHumanpurple/mean(simTOMHumanpurple)+ 
                                                      simTOMChimppurple/mean(simTOMChimppurple)+0.00001)


###### PURPLE ######################
corrplot(simRatiopurple, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_purple <- simRatiopurple > 0.65
### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat_purple <- bigratio_purple
dat_purple <- 1*dat_purple
head(dat_purple)
corrplot(dat_purple, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_purple <- simRatiopurple < 0.35
low_purple <- lowratio_purple
low_purple <- 1*low_purple
head(low_purple)
corrplot(low_purple, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))


