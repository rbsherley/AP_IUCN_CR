#======================================
# INSTALL the latest version of JARA if you need to - requires devtools - and the other required libraries
#install.packages(devtools)
#library(devtools)
#install_github("Henning-Winker/JARA",force=T)
#install.packages(HDInterval)
#install.packages(curl)
#======================================
library(JARA)
library(HDInterval)
library(curl)

# Set your working directory
setwd("~/YourDir")

## Run this code if you want to set the random seed to reproduce the values in the published ms - if not, you can run with a new random seed to robustness check the results
curl_download("https://raw.githubusercontent.com/rbsherley/AP_IUCN_CR/master/Seed.rds", destfile = "Seed.rds")
Seed <- readRDS("Seed.rds")
.Random.seed=Seed


# Data preparation
data <- read.csv("Afr_penguin_2023_JARAinput.csv") # all counts from 1979
data.g <- data[1:nrow(data),] # subset to global
data.sa <- data[1:nrow(data),1:20] # South Africa only
data.wc <- data[1:nrow(data),1:8] # Orange River to Cape Point only 
data.sc <- data[1:nrow(data),c(1,9:14)] # Cape Point to Cape Agulhas only 
data.ec <- data[1:nrow(data),c(1,15:20)] # Algoa Bay only 
data.na <- data[1:nrow(data),c(1,21:27)] # Namibia only 

# JAGS settings
ni<-25000 
nb<-10000
nt<-5
nc<-3

####################################################################
########## FITS TO OBSERVED DATA FOR IUCN RL CRITERION A2 ##########
####################################################################

##########
# Global #
##########
# create JARA input
ap.input = build_jara(I=data.g, model.type = "census",GL=10, assessment = "Afr_penguin",scenario = "Global", proj.r="GL1")
#set folder to save outputs
output.dir<-paste0(getwd(),"/AP output/",ap.input$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
# fit JARA model
ap.fit = fit_jara(ap.input, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)

pdf(file=paste0(getwd(),"/AP output/",ap.input$settings$scenario,".pdf"), width=6,height = 8)
par(mfrow=c(2,1),family="Times",mar=c(3,3,0.1,0.1),mgp=c(2,1,0))
jrplot_iucn(ap.fit,add=T)
jrplot_poptrj(ap.fit,add=T)
dev.off()

jrplot_iucn(ap.fit, Plot=FALSE)
jrplot_state(ap.fit, ref.yr = c(1979:1981))

################
# South Africa #
################
# create JARA input
ap.input.sa = build_jara(I=data.sa, model.type = "census",GL=10, assessment = "Afr_penguin",scenario = "S.Africa", proj.r="GL1")
#set folder to save outputs
output.dir<-paste0(getwd(),"/AP output/",ap.input.sa$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
# fit JARA model
ap.fit.sa = fit_jara(ap.input.sa, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)

pdf(file=paste0(getwd(),"/AP output/",ap.input.sa$settings$scenario,".pdf"), width=6,height = 8)
par(mfrow=c(2,1),family="Times",mar=c(3,3,0.1,0.1),mgp=c(2,1,0))
jrplot_iucn(ap.fit.sa,add=T)
jrplot_poptrj(ap.fit.sa,add=T)
dev.off()

jrplot_iucn(ap.fit.sa, Plot=FALSE)
jrplot_state(ap.fit.sa, ref.yr = c(1979:1981))


##############
# West Coast #
##############
# create JARA input
ap.input.wc = build_jara(I=data.wc, model.type = "census",GL=10, assessment = "Afr_penguin",scenario = "W_Coast", proj.r="GL1")
#set folder to save outputs
output.dir<-paste0(getwd(),"/AP output/",ap.input.wc$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
# fit JARA model
ap.fit.wc = fit_jara(ap.input.wc, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)

pdf(file=paste0(getwd(),"/AP output/",ap.input.wc$settings$scenario,".pdf"), width=6,height = 8)
par(mfrow=c(2,1),family="Times",mar=c(3,3,0.1,0.1),mgp=c(2,1,0))
jrplot_iucn(ap.fit.wc,add=T)
jrplot_poptrj(ap.fit.wc,add=T)
dev.off()

jrplot_iucn(ap.fit.wc, Plot=FALSE)
jrplot_state(ap.fit.wc, ref.yr = c(1979:1981))

### INDIVIDUAL COLONY PLOTS #######
fit = ap.fit.wc
years = fit$pyr
indices = fit$indices
n.indices = length(indices)
colors <- c(rep('#e6194b',length(indices)))
Par = list(mfrow=c(round(n.indices/2),2),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.5, 0.5, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
m1 <- 0
cs = sample(seq(80,90,1))
cols=paste0("gray",cs)
col_line <- colors

pdf(file = paste0(getwd(),"/AP output/",fit$settings$scenario,"_colonies.pdf"), width=6,height = 8)
par(Par)
for(i in 1:n.indices)
{
  m2 <- max(c(fit$trj$mu[fit$trj$name==indices[i]], fit$trj$uci[fit$trj$name==indices[i]]*1.1), na.rm = TRUE)
  plot(0, 0, ylim = c(m1, m2), xlim = c(min(fit$pyr-1),max(fit$pyr+1)), ylab = "Pairs", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
  polygon(x = c(years,rev(years)), y = c(fit$trj$lci[fit$trj$name==indices[i]],rev(fit$trj$uci[fit$trj$name==indices[i]])), col = "gray", border = "gray90")
  lines(years,fit$trj$mu[fit$trj$name==indices[i]], type = "l",col=col_line[i], lwd=1)
  points(fit$inputseries$I[,1],(fit$inputseries$I[,i+1]), bg = col_line[i],pch=21)
  posl = c(max(fit$trj$mu[fit$trj$name==indices[i]][1]),max(fit$trj$mu[fit$trj$name==indices[i]][length(years)]))
  legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(runs), legend = c(paste(indices[i])),pch=c(-1), pt.bg = c(col_line[i]),col = c(col_line[i]), bty = "n", cex = 0.9,y.intersp = 0.9)
  axis(1,at=seq(min(years),max(years)+5,5),tick=seq(min(years),max(years),5))
}
dev.off()
##################################################################################################


###############
# South Coast #
###############
# create JARA input
ap.input.sc = build_jara(I=data.sc, model.type = "census",GL=10, assessment = "Afr_penguin",scenario = "S_Coast", proj.r="GL1")
#set folder to save outputs
output.dir<-paste0(getwd(),"/AP output/",ap.input.sc$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
# fit JARA model
ap.fit.sc = fit_jara(ap.input.sc, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)

pdf(file=paste0(getwd(),"/AP output/",ap.input.sc$settings$scenario,".pdf"), width=6,height = 8)
par(mfrow=c(2,1),family="Times",mar=c(3,3,0.1,0.1),mgp=c(2,1,0))
jrplot_iucn(ap.fit.sc,add=T)
jrplot_poptrj(ap.fit.sc,add=T)
dev.off()

jrplot_iucn(ap.fit.sc, Plot=FALSE)
jrplot_state(ap.fit.sc, ref.yr = c(1979:1981))

### INDIVIDUAL COLONY PLOTS #######
fit = ap.fit.sc
years = fit$pyr
indices = fit$indices
n.indices = length(indices)
colors <- c(rep("#f58231",length(indices)))
Par = list(mfrow=c(round(n.indices/2),2),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.5, 0.5, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
m1 <- 0
cs = sample(seq(80,90,1))
cols=paste0("gray",cs)
col_line <- colors

pdf(file = paste0(getwd(),"/AP output/",fit$settings$scenario,"_colonies.pdf"), width=6,height = 8)
par(Par)
for(i in 1:n.indices)
{
  m2 <- max(c(fit$trj$mu[fit$trj$name==indices[i]], fit$trj$uci[fit$trj$name==indices[i]]*1.1), na.rm = TRUE)
  plot(0, 0, ylim = c(m1, m2), xlim = c(min(fit$pyr-1),max(fit$pyr+1)), ylab = "Pairs", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
  polygon(x = c(years,rev(years)), y = c(fit$trj$lci[fit$trj$name==indices[i]],rev(fit$trj$uci[fit$trj$name==indices[i]])), col = "gray", border = "gray90")
  lines(years,fit$trj$mu[fit$trj$name==indices[i]], type = "l",col=col_line[i], lwd=1)
  points(fit$inputseries$I[,1],(fit$inputseries$I[,i+1]), bg = col_line[i],pch=21)
  posl = c(max(fit$trj$mu[fit$trj$name==indices[i]][1]),max(fit$trj$mu[fit$trj$name==indices[i]][length(years)]))
  legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(runs), legend = c(paste(indices[i])),pch=c(-1), pt.bg = c(col_line[i]),col = c(col_line[i]), bty = "n", cex = 0.9,y.intersp = 0.9)
  axis(1,at=seq(min(years),max(years)+5,5),tick=seq(min(years),max(years),5))
}
dev.off()
##################################################################################################

################
# Eastern Cape #
################
# create JARA input
ap.input.ec = build_jara(I=data.ec, model.type = "census",GL=10, assessment = "Afr_penguin",scenario = "E_Cape", proj.r="GL1")
#set folder to save outputs
output.dir<-paste0(getwd(),"/AP output/",ap.input.ec$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
# fit JARA model
ap.fit.ec = fit_jara(ap.input.ec, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)

pdf(file=paste0(getwd(),"/AP output/",ap.input.ec$settings$scenario,".pdf"), width=6,height = 8)
par(mfrow=c(2,1),family="Times",mar=c(3,3,0.1,0.1),mgp=c(2,1,0))
jrplot_iucn(ap.fit.ec,add=T)
jrplot_poptrj(ap.fit.ec,add=T)
dev.off()

jrplot_iucn(ap.fit.ec, Plot=FALSE)
jrplot_state(ap.fit.ec, ref.yr = c(1979:1981))

### INDIVIDUAL COLONY PLOTS #######
fit = ap.fit.ec
years = fit$pyr
indices = fit$indices
n.indices = length(indices)
colors <- c(rep("#f032e6",length(indices)))
Par = list(mfrow=c(round(n.indices/2),2),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.5, 0.5, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
m1 <- 0
cs = sample(seq(80,90,1))
cols=paste0("gray",cs)
col_line <- colors

pdf(file = paste0(getwd(),"/AP output/",fit$settings$scenario,"_colonies.pdf"), width=6,height = 8)
par(Par)
for(i in 1:n.indices)
{
  m2 <- max(c(fit$trj$mu[fit$trj$name==indices[i]], fit$trj$uci[fit$trj$name==indices[i]]*1.1), na.rm = TRUE)
  plot(0, 0, ylim = c(m1, m2), xlim = c(min(fit$pyr-1),max(fit$pyr+1)), ylab = "Pairs", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
  polygon(x = c(years,rev(years)), y = c(fit$trj$lci[fit$trj$name==indices[i]],rev(fit$trj$uci[fit$trj$name==indices[i]])), col = "gray", border = "gray90")
  lines(years,fit$trj$mu[fit$trj$name==indices[i]], type = "l",col=col_line[i], lwd=1)
  points(fit$inputseries$I[,1],(fit$inputseries$I[,i+1]), bg = col_line[i],pch=21)
  posl = c(max(fit$trj$mu[fit$trj$name==indices[i]][1]),max(fit$trj$mu[fit$trj$name==indices[i]][length(years)]))
  legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(runs), legend = c(paste(indices[i])),pch=c(-1), pt.bg = c(col_line[i]),col = c(col_line[i]), bty = "n", cex = 0.9,y.intersp = 0.9)
  axis(1,at=seq(min(years),max(years)+5,5),tick=seq(min(years),max(years),5))
}
dev.off()
##################################################################################################


###########
# Namibia #
###########
# create JARA input
ap.input.na = build_jara(I=data.na, model.type = "census",GL=10, assessment = "Afr_penguin",scenario = "Namibia", proj.r="GL1")
#set folder to save outputs
output.dir<-paste0(getwd(),"/AP output/",ap.input.na$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
# fit JARA model
ap.fit.na = fit_jara(ap.input.na, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)

pdf(file=paste0(getwd(),"/AP output/",ap.input.na$settings$scenario,".pdf"), width=6,height = 8)
par(mfrow=c(2,1),family="Times",mar=c(3,3,0.1,0.1),mgp=c(2,1,0))
jrplot_iucn(ap.fit.na,add=T)
jrplot_poptrj(ap.fit.na,add=T)
dev.off()

jrplot_iucn(ap.fit.na, Plot=FALSE)
jrplot_state(ap.fit.na, ref.yr = c(1985:1987))

### INDIVIDUAL COLONY PLOTS #######
fit = ap.fit.na
years = fit$pyr
indices = fit$indices
n.indices = length(indices)
colors <- c(rep("#3cb44b",length(indices)))
Par = list(mfrow=c(round(n.indices/2),2),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.5, 0.5, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7)
m1 <- 0
cs = sample(seq(80,90,1))
cols=paste0("gray",cs)
col_line <- colors

pdf(file = paste0(getwd(),"/AP output/",fit$settings$scenario,"_colonies.pdf"), width=6,height = 8)
par(Par)
for(i in 1:n.indices)
{
  m2 <- max(c(fit$trj$mu[fit$trj$name==indices[i]], fit$trj$uci[fit$trj$name==indices[i]]*1.1), na.rm = TRUE)
  plot(0, 0, ylim = c(m1, m2), xlim = c(min(fit$pyr-1),max(fit$pyr+1)), ylab = "Pairs", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
  polygon(x = c(years,rev(years)), y = c(fit$trj$lci[fit$trj$name==indices[i]],rev(fit$trj$uci[fit$trj$name==indices[i]])), col = "gray", border = "gray90")
  lines(years,fit$trj$mu[fit$trj$name==indices[i]], type = "l",col=col_line[i], lwd=1)
  points(fit$inputseries$I[,1],(fit$inputseries$I[,i+1]), bg = col_line[i],pch=21)
  posl = c(max(fit$trj$mu[fit$trj$name==indices[i]][1]),max(fit$trj$mu[fit$trj$name==indices[i]][length(years)]))
  legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(runs), legend = c(paste(indices[i])),pch=c(-1), pt.bg = c(col_line[i]),col = c(col_line[i]), bty = "n", cex = 0.9,y.intersp = 0.9)
  axis(1,at=seq(min(years),max(years)+5,5),tick=seq(min(years),max(years),5))
}
dev.off()
##################################################################################################


##########################################################
########## PROJECTIONS FOR IUCN RL CRITERION A4 ##########
##########################################################

##############################################################################
# Global plus check of RL status for 1 to 10 yr projections based on last GL #
##############################################################################
hc = list(posteiors=NULL)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit$fits[1,1],level=0,pop.change=ap.fit$posteriors[,1]))

# Check decline from 1995 to 2024 - 1 year into the future
data.p <- data[16:nrow(data),] # subset to 3GL of past data before target terminal year
ap.input.1 = build_jara(I=data.p, model.type = "census", assessment = "Afr_penguin",scenario = 1, proj.r="GL1", proj.yrs.user = 2)
output.dir = paste0(getwd(),"/AP output/",ap.input.1$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
ap.input.1$data$pyr[ap.input.1$settings$mp.assess[1]] # 1995
ap.input.1$data$pyr[ap.input.1$settings$mp.assess[2]] # 2024
ap.fit.1 = fit_jara(ap.input.1, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit.1$fits[1,1],level=ap.fit.1$fits[1,2],pop.change=ap.fit.1$posteriors[,1]))

# Check decline from 1996 to 2025 - 2 years into the future
data.p <- data[17:nrow(data),] # subset to 3GL of past data before target terminal year
ap.input.2 = build_jara(I=data.p, model.type = "census", assessment = "Afr_penguin",scenario = 2, proj.r="GL1", proj.yrs.user = 3)
output.dir = paste0(getwd(),"/AP output/",ap.input.2$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
ap.input.2$data$pyr[ap.input.2$settings$mp.assess[1]] # 1996
ap.input.2$data$pyr[ap.input.2$settings$mp.assess[2]] # 2025
ap.fit.2 = fit_jara(ap.input.2, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit.2$fits[1,1],level=ap.fit.2$fits[1,2],pop.change=ap.fit.2$posteriors[,1]))

# Check decline from 1997 to 2026 - 3 years into the future
data.p <- data[18:nrow(data),] # subset to 3GL of past data before target terminal year
ap.input.3 = build_jara(I=data.p, model.type = "census", assessment = "Afr_penguin",scenario = 3, proj.r="GL1", proj.yrs.user = 4)
output.dir = paste0(getwd(),"/AP output/",ap.input.3$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
ap.input.3$data$pyr[ap.input.3$settings$mp.assess[1]] # 1997
ap.input.3$data$pyr[ap.input.3$settings$mp.assess[2]] # 2026
ap.fit.3 = fit_jara(ap.input.3, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit.3$fits[1,1],level=ap.fit.3$fits[1,2],pop.change=ap.fit.3$posteriors[,1]))

# Check decline from 1998 to 2027 - 4 years into the future
data.p <- data[19:nrow(data),] # subset to 3GL of past data before target terminal year
ap.input.4 = build_jara(I=data.p, model.type = "census", assessment = "Afr_penguin",scenario = 4, proj.r="GL1", proj.yrs.user = 5)
output.dir = paste0(getwd(),"/AP output/",ap.input.4$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
ap.input.4$data$pyr[ap.input.4$settings$mp.assess[1]] # 1998
ap.input.4$data$pyr[ap.input.4$settings$mp.assess[2]] # 2027
ap.fit.4 = fit_jara(ap.input.4, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit.4$fits[1,1],level=ap.fit.4$fits[1,2],pop.change=ap.fit.4$posteriors[,1]))

# Check decline from 1999 to 2028 - 5 years into the future
data.p <- data[20:nrow(data),] # subset to 3GL of past data before target terminal year
ap.input.5 = build_jara(I=data.p, model.type = "census", assessment = "Afr_penguin",scenario = 5, proj.r="GL1", proj.yrs.user = 6)
output.dir = paste0(getwd(),"/AP output/",ap.input.5$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
ap.input.5$data$pyr[ap.input.5$settings$mp.assess[1]] # 1999
ap.input.5$data$pyr[ap.input.5$settings$mp.assess[2]] # 2028
ap.fit.5 = fit_jara(ap.input.5, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit.5$fits[1,1],level=ap.fit.5$fits[1,2],pop.change=ap.fit.5$posteriors[,1]))

# Check decline from 2000 to 2029 - 6 years into the future
data.p <- data[21:nrow(data),] # subset to 3GL of past data before target terminal year
ap.input.6 = build_jara(I=data.p, model.type = "census", assessment = "Afr_penguin",scenario = 6, proj.r="GL1", proj.yrs.user = 7)
output.dir = paste0(getwd(),"/AP output/",ap.input.6$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
ap.input.6$data$pyr[ap.input.6$settings$mp.assess[1]] # 2000
ap.input.6$data$pyr[ap.input.6$settings$mp.assess[2]] # 2029
ap.fit.6 = fit_jara(ap.input.6, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit.6$fits[1,1],level=ap.fit.6$fits[1,2],pop.change=ap.fit.6$posteriors[,1]))

# Check decline from 2001 to 2030 - 7 years into the future
data.p <- data[22:nrow(data),] # subset to 3GL of past data before target terminal year
ap.input.7 = build_jara(I=data.p, model.type = "census", assessment = "Afr_penguin",scenario = 7, proj.r="GL1", proj.yrs.user = 8)
output.dir = paste0(getwd(),"/AP output/",ap.input.7$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
ap.input.7$data$pyr[ap.input.7$settings$mp.assess[1]] # 2001
ap.input.7$data$pyr[ap.input.7$settings$mp.assess[2]] # 2030
ap.fit.7 = fit_jara(ap.input.7, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit.7$fits[1,1],level=ap.fit.7$fits[1,2],pop.change=ap.fit.7$posteriors[,1]))

# Check decline from 2002 to 2031 - 8 years into the future
data.p <- data[23:nrow(data),] # subset to 3GL of past data before target terminal year
ap.input.8 = build_jara(I=data.p, model.type = "census", assessment = "Afr_penguin",scenario = 8, proj.r="GL1", proj.yrs.user = 9)
output.dir = paste0(getwd(),"/AP output/",ap.input.8$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
ap.input.8$data$pyr[ap.input.8$settings$mp.assess[1]] # 2002
ap.input.8$data$pyr[ap.input.8$settings$mp.assess[2]] # 2031
ap.fit.8 = fit_jara(ap.input.8, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit.8$fits[1,1],level=ap.fit.8$fits[1,2],pop.change=ap.fit.8$posteriors[,1]))

# Check decline from 2003 to 2032 - 9 years into the future
data.p <- data[24:nrow(data),] # subset to 3GL of past data before target terminal year
ap.input.9 = build_jara(I=data.p, model.type = "census", assessment = "Afr_penguin",scenario = 9, proj.r="GL1", proj.yrs.user = 10)
output.dir = paste0(getwd(),"/AP output/",ap.input.9$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
ap.input.9$data$pyr[ap.input.9$settings$mp.assess[1]] # 2003
ap.input.9$data$pyr[ap.input.9$settings$mp.assess[2]] # 2032
ap.fit.9 = fit_jara(ap.input.9, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit.9$fits[1,1],level=ap.fit.9$fits[1,2],pop.change=ap.fit.9$posteriors[,1]))

# Check decline from 2004 to 2033 - 10 years (1 GL) into the future
data.p <- data[25:nrow(data),] # subset to 2GL of past data
ap.input.10 = build_jara(I=data.p, model.type = "census", assessment = "Afr_penguin",scenario = 10, proj.r="GL1", proj.yrs.user = 11)
output.dir = paste0(getwd(),"/AP output/",ap.input.10$settings$scenario)
dir.create(output.dir,showWarnings = F,recursive=T)
ap.input.10$data$pyr[ap.input.10$settings$mp.assess[1]] # 2004
ap.input.10$data$pyr[ap.input.10$settings$mp.assess[2]] # 2033
ap.fit.10 = fit_jara(ap.input.10, ni=ni, nb=nb, nt=nt, nc=nc,save.jara=T,save.all=T,save.csvs = T,output.dir=output.dir)
hc$posteriors= rbind(hc$posteriors,data.frame(factor=ap.fit.10$fits[1,1],level=ap.fit.10$fits[1,2],pop.change=ap.fit.10$posteriors[,1]))


###### Plot of 10 year projection scenario ######
pdf(file=paste0(getwd(),"/AP output/",ap.input.10$settings$scenario,".pdf"), width=6,height = 8)
par(mfrow=c(2,1),family="Times",mar=c(3,3,0.1,0.1),mgp=c(2,1,0))
jrplot_iucn(ap.fit.10, add=T)
jrplot_poptrj(ap.fit.10,add=T)
lines(rep(2033,2),c(0,58000),lty=2,col=4,lwd=2)
lines(rep(2013,2),c(0,58000),lty=2,col=4,lwd=2)
lines(rep(2003,2),c(0,58000),lty=2,col=3,lwd=2)
text(2033,60000,"+1G")
text(2013,60000,"1G")
text(2003,60000,"2G")
dev.off()

###### INDIVIDUAL COLONY PLOTS #######
fit = ap.fit.10
years = fit$pyr
indices = fit$indices
n.indices = length(indices)
colors <- c(rep("#f58231",length(indices)))
Par = list(mfrow=c(4,2),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.5, 0.5, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.7,family="Times")
m1 <- 0
cs = sample(seq(80,90,1))
cols=paste0("gray",cs)
col_line <- colors

pdf(file = paste0(getwd(),"/AP output/",fit$settings$scenario,"_colonies.pdf"), width=6,height = 8)
par(Par)
for(i in 1:n.indices)
{
  m2 <- max(c(fit$trj$mu[fit$trj$name==indices[i]], fit$trj$uci[fit$trj$name==indices[i]]*1.1), na.rm = TRUE)
  plot(0, 0, ylim = c(m1, m2), xlim = c(min(fit$pyr-1),max(fit$pyr+1)), ylab = "Pairs", xlab = "Year", col = "black", type = "n", lwd = 2, frame = FALSE,xaxs="i",yaxs="i",xaxt="n")
  polygon(x = c(years,rev(years)), y = c(fit$trj$lci[fit$trj$name==indices[i]],rev(fit$trj$uci[fit$trj$name==indices[i]])), col = "gray", border = "gray90")
  lines(years,fit$trj$mu[fit$trj$name==indices[i]], type = "l",col=col_line[i], lwd=1)
  points(fit$inputseries$I[,1],(fit$inputseries$I[,i+1]), bg = col_line[i],pch=21)
  posl = c(max(fit$trj$mu[fit$trj$name==indices[i]][1]),max(fit$trj$mu[fit$trj$name==indices[i]][length(years)]))
  legend(ifelse(posl[1]<posl[2],"topleft","topright"),paste(runs), legend = c(paste(indices[i])),pch=c(-1), pt.bg = c(col_line[i]),col = c(col_line[i]), bty = "n", cex = 0.9,y.intersp = 0.9)
  axis(1,at=seq(min(years),max(years)+5,5),tick=seq(min(years),max(years),5))
  
  lines(rep(2033,2),c(0,m2*.9),lty=2,col=2,lwd=2)
  lines(rep(2013,2),c(0,m2*.9),lty=2,col=4,lwd=2)
  lines(rep(2023,2),c(0,m2*.9),lty=2,col=1,lwd=2)
  lines(rep(2003,2),c(0,m2*.9),lty=2,col=3,lwd=2)
  
}
dev.off()

###### PLOT OF CHANGE IN RL STATUS OVER PROJECTION WINDOW #######
A1 = FALSE
runs = hc$peels = 0:10
hc$yr = 2023:2033
d = hc$posteriors
ymax1 = quantile(hc$posteriors$pop.change[hc$posteriors$level==0],0.995)
ymax2 = quantile(hc$posteriors$pop.change[hc$posteriors$level>0],0.985)
ymax = max(ymax1,ymax2)
ylim=c(-100,-48)
xlim=c(0.5,length(hc$peels)+0.49)
xall = hc$posteriors$pop.change
cols = c("#60C659","lightgreen","#F9E814","#FC7F3F","#D81E05")[c(1,3:5)] # green to red

pdf(file = "Projections_IUCN.pdf", width=6,height = 3)
#
#quartz(width=6,height = 3)
Par = list(mfrow=c(1,1),mar = c(3, 3, 1, 1), mgp =c(2.5,1,0),mai = c(0.5, 0.5, 0.1, 0.1),mex=0.8, tck = -0.02,cex=0.9,family="Times")
par(Par)
plot(0,0,type="n",xlab="",xlim=xlim,ylim=ylim,axes=F,xaxs = "i",yaxs="i",ylab="")  
axis(1,at=seq(0,length(hc$peels),1),labels=max(hc$yr)-rev(c(hc$peels,max(hc$peels+1))),cex.axis=0.8,mgp=c(2,0.5,0))
axis(2,at=seq(-100,-50,10),tick=seq(-100,-50,10),cex.axis=0.8,mgp=c(2,0.5,0))
mtext(side=2,"Change (%)",line=1.9)
mtext(side=1,"Assessment terminal year",line=1.9)
box()
out = NULL
for(j in 1:length(runs)){
  change = d[d$level == runs[j],]$pop.change
  change=ifelse(change>ymax,max(ymax+50,45),change)
  den = stats::density(change,adjust=1)

  #y1 = exp(den$x)-100
  y1 = den$x
  x1 = den$y/max(den$y)+j-0.5
  # get status
  categories = c("Pr.Decl","change3GL","CR","EN","VU","LC")  
  CR = round(sum(ifelse(change< ifelse(A1,-90,-80),1,0))/length(change)*100,1)
  EN = round(sum(ifelse(change> ifelse(A1,-90,-80) & change< ifelse(A1,-70,-50),1,0))/length(change)*100,1)
  VU = round(sum(ifelse(change> ifelse(A1,-70,-50) & change< ifelse(A1,-50,-30),1,0))/length(change)*100,1)
  LC = round(sum(ifelse(change> -30,1,0))/length(change)*100,1)
  Change3xGL = round(median(change),3)# round(sum(ifelse(change< 0,1,0))/length(change)*100,1)
  percentages = c(CR,EN,VU,LC)
  status= ifelse(which(percentages==max(percentages))==4 & max(percentages)<50,"NT",categories[3:6][which(percentages==max(percentages))])
  out = rbind(out,data.frame(Year=max(hc$yr)-rev(runs)[j], Change3xGL,CR,EN,VU,LC,status))
  
  lc = c(ifelse(A1,-50,-30))
  polygon(c(x1[y1>=lc],rep(min(x1),length(y1[y1>=lc]))),c(y1[y1>=lc],rev(y1[y1>=lc])),col=cols[1] ,border=cols[1])
  vu = c(ifelse(A1,-70,-50))
  polygon(c(x1[y1<=lc & y1>=vu],rep(min(x1),length(y1[y1<=lc & y1>=vu]))),c(y1[y1<lc & y1>=vu],rev(y1[y1<lc & y1>=vu])),col=cols[2],border=cols[2])
  en =ifelse(A1,-90,-80)
  polygon(c(x1[y1<vu & y1>=en],rep(min(x1),length(y1[y1<vu & y1>=en]))),c(y1[y1<vu & y1>=en],rev(y1[y1<vu & y1>=en])),col=cols[3],border=cols[3])
  polygon(c(x1[y1<en],rep(min(x1),length(y1[y1<en]))),c(y1[y1<en],rev(y1[y1<en])),col=cols[4],border=cols[4])
  polygon(c(x1,rep(min(x1),length(x1))),c(y1,rev(y1)))
  
}

text(1:length(runs),par('usr')[4],(out$status),cex=0.8,pos=1,offset = 0.2)
text(1:length(runs),par('usr')[4],round(out$CR,1),cex=0.8,pos=1,offset = 1)
text(1:length(runs),par('usr')[4],round(out$Change3xGL,1),cex=0.8,pos=1,offset = 12.5)
dev.off()


pdf(file="Global_10_project.pdf", width=4,height=4)
#quartz(width=4,height=4)
par(mfrow=c(1,1),family="Times",mar=c(2.5,2.5,0.1,0.1),mgp=c(2,.7,0),yaxt="n",las=1)
jrplot_poptrj(ap.fit.10,add=T,ylab="",xlab="")
lines(rep(2033,2),c(0,58000),lty=2,col=4,lwd=2)
lines(rep(2013,2),c(0,58000),lty=2,col=4,lwd=2)
lines(rep(2003,2),c(0,58000),lty=2,col=3,lwd=2)
text(2033,60000,"+1G")
text(2013,60000,"1G")
text(2003,60000,"2G")
axis(2,at=seq(0,60000,by=10000),labels=seq(0,60,10),yaxt="s")
mtext(side=2,"Breeding pairs x 1000",line=1.6,cex=1.2,las=3)
mtext(side=1,"Year",line=1.6,cex=1.2)
dev.off()

GIUCN<- jrplot_iucn(ap.fit,Plot = FALSE)
NIUCN<- jrplot_iucn(ap.fit.na,Plot = FALSE)
SIUCN<- jrplot_iucn(ap.fit.sa,Plot = FALSE)

pdf(file="All_IUCN_plots.pdf", width=6.9,height=8.85)
#quartz(width=6.9,height=8.85)
par(mfrow=c(3,2),family="Times",mar=c(2.8,2.7,0.9,0.1),mgp=c(2,.75,0),yaxt="n",las=1,cex=1.1)
jrplot_poptrj(ap.fit,add=T,xlab="",ylab="",plot.cex=1.5)
axis(2,at=seq(0,100000,by=10000),labels=seq(0,100,10),yaxt="s")
mtext(side=2,"Breeding pairs x 1000",line=1.7,las=3,cex=1.2)
mtext(side=1,"Year",line=1.6,las=1,cex=1.2)
mtext(side=3,"A: Global",line=0,las=1,cex=1.2)
par(mar=c(2.8,0.5,0.1,0.1))
jrplot_iucn(ap.fit,add=T,xlim=c(-100,-10),ylimadj=1.05,xlab="",print_change=FALSE)
text(paste0("Change =",GIUCN$perc.risk[2],"%"),x=-60,y=0.135)
mtext(side=1,"Change (%)",line=1.6,las=1,cex=1.2)
# Namibia
par(mar=c(2.8,2.7,0.9,0.1),mgp=c(2,.75,0),yaxt="n")
jrplot_poptrj(ap.fit.na,add=T,xlab="",ylab="",plot.cex=1.5)
axis(2,at=seq(0,14000,by=2000),labels=seq(0,14,2),yaxt="s")
mtext(side=2,"Breeding pairs x 1000",line=1.7,las=3,cex=1.2)
mtext(side=1,"Year",line=1.6,las=1,cex=1.2)
mtext(side=3,"B: Namibia",line=0,las=1,cex=1.2)
par(mar=c(2.8,0.5,0.1,0.1))
jrplot_iucn(ap.fit.na,add=T,xlim=c(-100,-10),ylimadj=1.1,xlab="",print_change=FALSE)
text(paste0("Change =",NIUCN$perc.risk[2],"%"),x=-60,y=0.234)
mtext(side=1,"Change (%)",line=1.6,las=1,cex=1.2)
# South Africa
par(mar=c(2.8,2.7,0.9,0.1),mgp=c(2,.75,0),yaxt="n")
jrplot_poptrj(ap.fit.sa,add=T,xlab="",ylab="",plot.cex=1.5)
axis(2,at=seq(0,80000,by=10000),labels=seq(0,80,10),yaxt="s")
mtext(side=2,"Breeding pairs x 1000",line=1.7,las=3,cex=1.2)
mtext(side=1,"Year",line=1.5,las=1,cex=1.2)
mtext(side=3,"C: South Africa",line=0,las=1,cex=1.2)
par(mar=c(2.8,0.5,0.1,0.1))
jrplot_iucn(ap.fit.sa,add=T,xlim=c(-100,-10),ylimadj=1.1,xlab="",print_change=FALSE)
text(paste0("Change =",SIUCN$perc.risk[2],"%"),x=-60,y=0.121)
mtext(side=1,"Change (%)",line=1.6,las=1,cex=1.2)
dev.off()

