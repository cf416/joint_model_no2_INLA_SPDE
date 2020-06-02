remove(list = ls())

#=================================================== 
### Load necessary R packages
#===================================================

library(INLA)
library(ggplot2)
library(raster)
library(spacetime)
library(scoringRules)
library(gridExtra)
library(gstat)
library(rgeos)
library(RColorBrewer)
library(gdata)
library(viridis)

## Set seed for reproducibility
set.seed(123)

## Load useful functions
source("utility.R")


#=================================================== 
### Data preparation
#===================================================

load("workspace_data.RData")

pol = "no2"

final_dataset$site.type = factor(as.character(final_dataset$site.type), ordered = T, levels = c("RUR","URB","RKS"))
final_dataset$site.type.n = as.numeric(final_dataset$site.type)
final_dataset$sitetype.idx = as.numeric(final_dataset$site.type)

##--- Visualise data

plot(shape)
points(coordinates.pcm, pch=3, cex=.2)
points(coordinates.aqum, pch=19, cex=.5, col="red")
points(monitors[monitors$site.type=="RUR", c("easting","northing")], pch=8, cex=.5, col="green3")
points(monitors[monitors$site.type=="URB", c("easting","northing")], pch=15, cex=.5, col="orange")
points(monitors[monitors$site.type=="RKS", c("easting","northing")], pch=17, cex=.5, col="blue")
lines(london.shape)
par(xpd=TRUE)
legend("topleft",inset=c(-0.15,0.1),legend=c("PCM","AQUM","Rural monitors","Urban monitors",
                                             "Roadside/ \n Kerbside monitors"), col = c("black","red","green3","orange","blue"),
       pch=c(3,19,8,15,17), cex=.7, pt.cex=1.2, bty = "n")

formula = as.formula(paste0("y ~ -1 + alpha1 + alpha2 + alpha3 + stURB.",pol," + stRKS.",pol," +  
                                  f(z2, model='rw1', hyper = rw1.aqum.prior, scale.model = TRUE, constr = TRUE) +
                            f(z23, copy='z2', fixed = FALSE, hyper=lambda23) + 
                            f(z1, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) +
                            f(z12, copy='z1', fixed = FALSE, hyper=lambda12) + 
                            f(z13, copy='z1', fixed = FALSE, hyper=lambda13) + 
                            f(date.idx.",pol,", model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)"))

#=================================================== 
### Run models and extract results
#===================================================
##--- NOTE: data_id indicates the validation set and is received from the bash script
data_id <- ifelse(nchar(Sys.getenv("DATA_ID"))>0, as.numeric(Sys.getenv("DATA_ID")), 1) # 1 to 6

if (!file.exists(file.path(getwd(),paste0("Output_",data_id)))){
  dir.create(file.path(getwd(),paste0("Output_",data_id)))
} 
setwd(paste0("Output_",data_id))
print(getwd())

valid = final_dataset[final_dataset$code %in% monitors_val[[data_id]] , ]
estim = final_dataset[!(final_dataset$code %in% monitors_val[[data_id]]) , ]

coordinates.estim<-unique(estim[,c("loc.idx","easting","northing")])
coordinates.valid<-unique(valid[,c("loc.idx","easting","northing")])

n_monitors=nrow(coordinates.y)
n_data=nrow(estim)+nrow(valid)
n_days=n_data/n_monitors

#=================================================== 
### Create mesh
#===================================================

mesh = inla.mesh.2d(rbind(coordinates.y,coordinates.aqum,coordinates.pcm),
                    loc.domain=boundary@polygons[[1]]@Polygons[[1]]@coords,
                    max.edge = c(75000,40000),
                    offset = c(10000,30000),
                    cutoff=8000)
plot(mesh, main="")
lines(shape, col="blue")
lines(london.shape, col="blue")
title("Domain triangulation")
points(coordinates.estim, col="green")
points(coordinates.valid, col="red")


#=================================================== 
### Construct the SPDE model for Matern field with some prior information obtained from the mesh or the spatial domain 
#=================================================== 

range0 <- min(c(diff(range(mesh$loc[,1])),diff(range(mesh$loc[,2]))))/5 
spde <- inla.spde2.pcmatern(mesh=mesh, alpha=2, ### mesh and smoothness parameter
                            prior.range=c(range0, 0.95), ### P(practic.range<range0)=0.95
                            prior.sigma=c(100, 0.5)) ### P(sigma>100)=0.5

#=================================================== 
### Hyperpriors 
#=================================================== 

rw1.aqum.prior = list(theta=list(prior="pc.prec", param=c(sd(aqum[,paste0(pol,"_log")]),0.01)))
ar1.time.prior = list(theta2 = list(prior='normal', 
                                    param=c(inla.models()$latent$ar1$hyper$theta2$to.theta(0.3), 0.5)))

lambda23 = list(theta = list(prior = 'normal', param = c(0.9, 0.01), initial=0.9))  # lambda_2,3
lambda12 = list(theta = list(prior = 'normal', param = c(1.1, 0.01), initial=1.1))  # lambda_1,2
lambda13 = list(theta = list(prior = 'normal', param = c(1.3, 0.01), initial=1.3)) # lambda_1,3



#=================================================== 
### Stack 
#===================================================

## ***** PCM *****
A_pcm <- inla.spde.make.A(mesh=mesh, cbind(pcm$easting, pcm$northing))
stk_pcm <- inla.stack(data=list(y=cbind(pcm[,paste0(pol,"_log")], NA, NA)),  
                      effects=list(list(alpha3=rep(1,nrow(pcm))), list(z1=1:spde$n.spde)),
                      A=list(1,A_pcm),
                      tag="est.pcm")

## ***** AQUM *****
A_aqum <- inla.spde.make.A(mesh=mesh,cbind(aqum$easting, aqum$northing))
stk_aqum <- inla.stack(data=list(y=cbind(NA, aqum[,paste0(pol,"_log")], NA)),  
                       effects=list(list(alpha2=1, z2=aqum$date.idx),
                                    list(z12=1:spde$n.spde)), 
                       A=list(1, A_aqum), 
                       tag="est.aqum")

## data stack: include all the effects
A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
stk_y_e <- inla.stack(data=list(y=cbind(NA,NA,estim[,paste0(pol,"_log")])),  
                      effects=list(list(z23=estim[,paste0("date.idx.",pol)]),
                                   list(z13=1:spde$n.spde),
                                   data.frame(alpha1=1,
                                              estim[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing")])), 
                      A=list(1, A_y_e, 1),
                      tag="est.y")

### validation scenario
A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
stk_y_v <- inla.stack(data=list(y=cbind(NA,NA,rep(NA,length(valid[,paste0(pol,"_log")])))),
                      effects=list(list(z23=valid[,paste0("date.idx.",pol)]),
                                   list(z13=1:spde$n.spde),
                                   data.frame(alpha1=1,
                                              valid[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing")])),
                      A=list(1, A_y_v, 1),
                      tag="val.y")

stack <- inla.stack(stk_aqum, 
                    stk_pcm,
                    stk_y_v, 
                    stk_y_e) 


#=================================================== 
### INLA call 
#===================================================

##--- NOTE: the INLA call can be parallelized using Pardiso - see inla.pardiso()

mod<- inla(formula,
           family=c("gaussian","gaussian","gaussian"),
           data=inla.stack.data(stack),
           control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
           control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
           control.inla = list(adaptive,int.strategy='eb'),
           verbose=TRUE)

#=================================================== 
###  Summary of Posterior distributions of parameters and hyperparameters of interest
#===================================================

##--- FIXED EFFECTS
print(round(mod$summary.fixed,4))

##--- HYPERPARAMETERS
print(round(mod$summary.hyperpar,4))



#=================================================== 
###  Check model performance
#===================================================

n.failures = sum(mod$cpo$failure, na.rm = T)
if(mean(mod$cpo$failure, na.rm = T)!=0){
  # summary(mod$cpo$failure)
  #  Two options:
  #  1. recompute using control.inla=list(int.strategy = "grid", diff.logdens = 4, strategy = "laplace", npoints = 21) see http://www.r-inla.org/faq
  #  2. run mod = inla.cpo(mod) #recomputes in an efficient way the cpo/pit for which mod$cpo$failure > 0
  mod_imp = inla.cpo(mod)
  print("Model cpo have been improved")
  logscore = -mean(log(mod_imp$cpo$cpo), na.rm=T)
  par(mfrow=c(1,2))
  hist(mod_imp$cpo$pit, breaks=100, main=paste0("PIT improved - n.failures=", n.failures))
  hist(mod_imp$cpo$cpo, breaks=100, main=paste0("CPO improved - n.failures=", n.failures))
}else{
  logscore = -mean(log(mod$cpo$cpo), na.rm=T)
  par(mfrow=c(1,2))
  hist(mod$cpo$pit, breaks=100, main=paste0("PIT - n.failures=", n.failures))
  hist(mod$cpo$cpo, breaks=100, main=paste0("CPO - n.failures=", n.failures))
}


#=================================================== 
###  Extract posterior latent fields
#===================================================

print(names(mod$summary.random))
# Index for the temporal fields
index.temp.field = which(substr(names(mod$summary.random),1,2) == "z2")
print(index.temp.field)
# Index for the spatial fields
index.spat.field = which(substr(names(mod$summary.random),1,2) == "z1")
print(index.spat.field)

##--- Time-sitetype interaction

rur.mean <- ts(mod$summary.random$date.idx$mean[1:1826], start = c(2007, 1), frequency = 365)
rur.mean.low <- ts(mod$summary.random$date.idx$`0.025quant`[1:1826], start = c(2007, 1), frequency = 365)
rur.mean.upp <- ts(mod$summary.random$date.idx$`0.975quant`[1:1826], start = c(2007, 1), frequency = 365)
urb.mean <- ts(mod$summary.random$date.idx$mean[1827:(2*1826)], start = c(2007, 1), frequency = 365)
urb.mean.low <- ts(mod$summary.random$date.idx$`0.025quant`[1827:(2*1826)], start = c(2007, 1), frequency = 365)
urb.mean.upp <- ts(mod$summary.random$date.idx$`0.975quant`[1827:(2*1826)], start = c(2007, 1), frequency = 365)
rks.mean <- ts(mod$summary.random$date.idx$mean[(2*1826+1):(3*1826)], start = c(2007, 1), frequency = 365)
rks.mean.low <- ts(mod$summary.random$date.idx$`0.025quant`[(2*1826+1):(3*1826)], start = c(2007, 1), frequency = 365)
rks.mean.upp <- ts(mod$summary.random$date.idx$`0.975quant`[(2*1826+1):(3*1826)], start = c(2007, 1), frequency = 365)

lower = min(mod$summary.random$date.idx.no2$`0.025quant`)
upper = max(mod$summary.random$date.idx.no2$`0.975quant`)
png('time_sitetype_interaction.png', width = 8, height = 8, unit="in", res=600)
par(mfrow=c(3,1))
plot(rur.mean.upp, type="l", ylab=expression(paste(mu,"g/",m^3)), xlab="Year", ylim=c(lower,upper), main="RUR", lty="twodash", col="red")
lines(rur.mean.low, lty="twodash", col="red")
lines(rur.mean)
abline(h=0,col="blue")
plot(urb.mean.upp, type="l", ylab=expression(paste(mu,"g/",m^3)), xlab="Year", ylim=c(lower,upper), main="URB", lty="twodash", col="red")
lines(urb.mean.low, lty="twodash", col="red")
lines(urb.mean)
abline(h=0,col="blue")
plot(rks.mean.upp, type="l", ylab=expression(paste(mu,"g/",m^3)), xlab="Year", ylim=c(lower,upper), main="RKS", lty="twodash", col="red")
lines(rks.mean.low, lty="twodash", col="red")
lines(rks.mean)
abline(h=0,col="blue")
dev.off()


##--- Latent temporal fields

for(i in index.temp.field){
  aqum.ts.mean <- ts(mod$summary.random[[i]]$mean, start = c(2007, 1), frequency = 365)
  aqum.ts.low <- ts(mod$summary.random[[i]]$`0.025quant`, start = c(2007, 1), frequency = 365)
  aqum.ts.upp <- ts(mod$summary.random[[i]]$`0.975quant`, start = c(2007, 1), frequency = 365)
  aqum.ts.sd <- ts(mod$summary.random[[i]]$sd, start = c(2007, 1), frequency = 365)
  png(paste0(names(mod$summary.random)[i],'.png'), width = 10, height = 8, unit="in", res=600)
  par(mfrow=c(2,1))
  plot(aqum.ts.upp,
       main=paste0(names(mod$marginals.random)[i]," - mean and CI"),lty="twodash",
       xlab="Year",ylab=expression(paste(mu,"g/",m^3)),
       ylim=c(min(aqum.ts.low),max(aqum.ts.upp))
  )
  lines(aqum.ts.low, lty="twodash")
  lines(aqum.ts.mean,  col="red")
  plot(aqum.ts.sd,
       main=paste0(names(mod$marginals.random)[i]," - SD"),
       xlab="Year",ylab=expression(paste(mu,"g/",m^3)))
  dev.off()
}


##--- Latent spatial fields

for(i in index.spat.field){
  
  length.grid.x=150
  length.grid.y=150
  
  # construct a lattice over the mesh extent
  pred.grid.lat= inla.mesh.lattice(x=seq(extent(boundary)[1],extent(boundary)[2],length.out = length.grid.x),
                                   y=seq(extent(boundary)[3],extent(boundary)[4],length.out = length.grid.y))
  
  proj = inla.mesh.projector(mesh, lattice = pred.grid.lat)
  
  spat.field = cbind(expand.grid(x=seq(extent(boundary)[1],extent(boundary)[2],length.out = length.grid.x),
                                 y=seq(extent(boundary)[3],extent(boundary)[4],length.out = length.grid.y)),
                     mean.log=as.vector(inla.mesh.project(proj,field= mod$summary.random[[i]]$mean)),
                     sd.log=as.vector(inla.mesh.project(proj,field= mod$summary.random[[i]]$sd)))
  
  coordinates(spat.field) = ~x+y
  proj4string(spat.field) = proj4string(shape)
  spat.field = subset(spat.field, over(spat.field, shape)$objectid==1)
  spat.field = as.data.frame(spat.field)
  
  
  pcm_lat_log = ggplot(spat.field) +
    #gg(mesh.pcm) +
    geom_raster(aes(x, y, fill = mean.log)) +
    scale_fill_viridis(bquote(paste("log(",.(toupper(pol)),") (",mu,"g/",m^3,")"))) +
    ggtitle(paste0(names(mod$marginals.random)[i]," - posterior mean")) +
    geom_path(data = fortify(shape), aes(group = group, x = long, y = lat)) +
    geom_path(data = fortify(london.shape), aes(group = group, x = long, y = lat), size=0.5) +
    theme(plot.title = element_text(family = "sans", color="#666666", size=14, hjust=0.5, face="bold"),
          axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
          legend.key.height = unit(1.5,"cm"), 
          legend.text = element_text(size=7, vjust=-0.5)) +
    labs(x="Easting",y="Northing")
  
  pcm_lat_sd_log = ggplot(spat.field) +
    geom_raster(aes(x, y, fill = sd.log)) +
    #gg(mesh.pcm) +
    scale_fill_viridis(bquote(paste("log(",.(toupper(pol)),") (",mu,"g/",m^3,")"))) +
    ggtitle(paste0(names(mod$marginals.random)[i]," - posterior SD")) +
    geom_path(data = fortify(shape), aes(group = group, x = long, y = lat)) +
    geom_path(data = fortify(london.shape), aes(group = group, x = long, y = lat), size=0.5) +
    theme(plot.title = element_text(family = "sans", color="#666666", size=14, hjust=0.5, face="bold"),
          axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
          legend.key.height = unit(1.5,"cm"), 
          legend.text = element_text(size=7, vjust=-0.5)) +
    labs(x="Easting",y="Northing")
  
  ggsave(paste0(names(mod$marginals.random)[i],'.png'), plot=pcm_lat_log, width = 8, height = 8, units="in", dpi=600)
  
}


#=================================================== 
###  Predictive capability
#===================================================

sample.size=50 
n.pred=100
predictions= list()

# extract sample.size values from each marginal of the fitted values
posterior.samples = inla.posterior.sample(n = sample.size, result=mod, intern = FALSE, use.improved.mean = TRUE, add.names = TRUE, seed = 0L)
# extract the posterior latent for the validation sites and reshape in order to have a matrix of nrow(valid) x sample.size
samples.fitted <- matrix(unlist(lapply(lapply(posterior.samples,"[[", "latent" ),"[", inla.stack.index(stack,"val.y")$data), use.names=FALSE), nrow = nrow(valid), byrow = F)

# extract sample.size values from marginal of the likelihood variance
# (it is always the last of the Gaussian observations, so we extract the position)
par.index = length(which(substr(rownames(mod$summary.hyperpar),1,26) == "Precision for the Gaussian"))
samples.var = 1/unlist(lapply(lapply(posterior.samples,"[[", "hyperpar" ),"[", par.index))
# this is not exactly what we want, in theory we should sample from the transformed posterior marginal
# rather than inverting the values sampled from the posterior marginal of the precision

samples.fitted = cbind(samples.fitted, obs.value=valid[,paste0(pol,"_log")])

predictive.capability = apply(samples.fitted, MARGIN = 1, FUN=extract.predicted, variance=samples.var)

# To aggregate the results we keep code, date and site type variables
predictive.capability = cbind(code=valid$code, date.idx=valid[,paste0("date.idx.",pol)], site.type=valid$site.type, as.data.frame(t(predictive.capability)))

predictive.capability$CRPS = crps_sample(y=samples.fitted[,(ncol(samples.fitted)-1)], dat=as.matrix(samples.fitted[,1:(ncol(samples.fitted)-1)]))
predictive.capability$logscore = logscore

saveRDS(predictive.capability, "predictive_capability.rds")

##--- Run this chunk once all the models results have been extracted, so we have the predictions for all the monitors (the 6 validation sets)
files = list.files(path=list.dirs(path = "..", full.names = TRUE, recursive = TRUE), pattern ="predictive_capability.rds", full.names = TRUE)
predictive.capability = do.call(rbind, lapply(files, function (x) readRDS(x)))
predictive.capability.measures = c(my.validation(z = predictive.capability[,paste0(pol,"_obs")], 
                                                 zhat = predictive.capability[,paste0(pol,"_pred")], 
                                                 penalty = predictive.capability$pmcc_penalty,
                                                 coverage = predictive.capability$coverage),
                                   LOGSCORE = predictive.capability$logscore,
                                   COV_CI = mean(predictive.capability$ci_amplitude),
                                   CRPS = mean(predictive.capability$CRPS, na.rm = T))
# NOTE: predictive.capability.measures can be computed by site-type, by day or by site, 
#       subsetting predictive.capability by site.type, date.idx or code respectively.
##---

#=================================================== 
###  Extract daily predictions 
#===================================================

# NOTE: daily predictions here are extracted for each model; 
#       this chunk can be skipped and run only after re-running the model using all data (no validation)

n.samples = 50
n.days = 1826
n.locs = nrow(pred.grid)

A_pred =inla.spde.make.A(mesh, loc=cbind(pred.grid$easting, pred.grid$northing))

contents=as.data.frame(mod$misc$configs$contents)
contents$end = contents$start + contents$length - 1

fields = lapply(1:n.samples, FUN=extract.contents)

if (!file.exists(file.path(getwd(), "predictions_by_day"))){
  dir.create(file.path(getwd(), "predictions_by_day"))
}

date = data.frame(date=seq.Date(as.Date("2007-01-01"), as.Date("2011-12-31"), "days"),
                  date.idx.no2=c(1:n_days),
                  year=as.numeric(substr(unique(aqum$date),1,4)),
                  month=as.numeric(substr(unique(aqum$date),6,7)),
                  day=as.numeric(substr(unique(aqum$date),9,10)))

invisible(lapply(1:n.days, FUN=compute.daily.predictions))

##---  predictions for days of pollution events

pred_2007_12_11 = readRDS("predictions_by_day/predictions_2007-12-11.rds")
pred_2007_12_19 = readRDS("predictions_by_day/predictions_2007-12-19.rds")
pred_2009_01_03 = readRDS("predictions_by_day/predictions_2009-01-03.rds")
pred_2010_11_16 = readRDS("predictions_by_day/predictions_2010-11-16.rds")

##---  predictions for days of low pollution

pred_2007_06_24 = readRDS("predictions_by_day/predictions_2007-06-24.rds")
pred_2008_06_22 = readRDS("predictions_by_day/predictions_2008-06-22.rds")
pred_2009_06_21 = readRDS("predictions_by_day/predictions_2009-06-21.rds") 
pred_2010_06_20 = readRDS("predictions_by_day/predictions_2010-06-20.rds")

pred_events = rbind(cbind(mean = rowMeans(pred_2007_12_11[,-c(1:3)], na.rm = T), day = "2007-12-11", pred.grid),
                    cbind(mean = rowMeans(pred_2007_12_19[,-c(1:3)], na.rm = T), day = "2007-12-19", pred.grid),
                    cbind(mean = rowMeans(pred_2009_01_03[,-c(1:3)], na.rm = T), day = "2009-01-03", pred.grid),
                    cbind(mean = rowMeans(pred_2010_11_16[,-c(1:3)], na.rm = T), day = "2010-11-16", pred.grid),
                    cbind(mean = rowMeans(pred_2007_06_24[,-c(1:3)], na.rm = T), day = "2007-06-24", pred.grid),
                    cbind(mean = rowMeans(pred_2008_06_22[,-c(1:3)], na.rm = T), day = "2008-06-22", pred.grid),
                    cbind(mean = rowMeans(pred_2009_06_21[,-c(1:3)], na.rm = T), day = "2009-06-21", pred.grid),
                    cbind(mean = rowMeans(pred_2010_06_20[,-c(1:3)], na.rm = T), day = "2010-06-20", pred.grid))

ggsave("daily_predictions.png",
       ggplot(pred_events) +
         geom_raster(aes(x=easting, y=northing, fill = mean))+
         facet_wrap(~day, ncol=4) +
         scale_fill_viridis(bquote(paste("log(",.(toupper(pol)),") (",mu,"g/",m^3,")"))) +
         ggtitle("Daily predictions on selected days with air pollution events (top row) or low concentration (bottom row)") +
         geom_path(data = fortify(london.shape), aes(group = group, x = long, y = lat), size=0.5)+
         geom_path(data = fortify(shape), aes(group = group, x = long, y = lat))+  
         geom_path(data = fortify(roads_major), aes(group = group, x = long, y = lat), size=0.3, color="grey85")+
         theme(plot.title = element_text(family = "sans", color="#666666", size=14, hjust=0.5, face="bold"),
               axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
               legend.key.height = unit(1,"cm"),
               legend.text = element_text(size=7, vjust=-0.5)),
       width=15, height=8, unit="in", dpi=600)