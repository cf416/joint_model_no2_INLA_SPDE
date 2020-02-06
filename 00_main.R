remove(list = ls())

############################################################################
########                      load R packages                      #########
############################################################################

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
library(spAir)


############################################################################
########                    useful functions                       #########
############################################################################


MSE <- function(z, zhat) {
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  x <- c(z - zhat)
  u <- x[!is.na(x)]
  # round(sqrt(sum(u^2)/length(u)), 4) # cannot be rooted now because we take the average later
  round(sum(u^2)/length(u), 4)
}
MAE <- function(z, zhat) {
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  x <- abs(c(zhat - z))
  u <- x[!is.na(x)]
  round(sum(u)/length(u), 4)
}
MAPE <- function(z, zhat) {
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  x <- abs(c(zhat - z))/z
  u <- x[!is.na(x)]
  u <- u[!is.infinite(u)]
  round(sum(u)/length(u) * 100, 4)
}
BIAS <- function(z, zhat) {
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  x <- c(zhat - z)
  u <- x[!is.na(x)]
  round(sum(u)/length(u), 4)
}
pBIAS <- function(z, zhat) {
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  x <- c(zhat - z)/z
  u <- x[!is.na(x)]
  u <- u[!is.infinite(u)]
  round(sum(u)/length(u) * 100, 4)
}
CORR <- function(z, zhat) {
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  round(cor(z,zhat,use="pairwise.complete.obs", method="pearson"), 4)
}
COV <- function(z, lower=NULL, upper=NULL) {
  if(!is.null(lower) && !is.null(upper)){
    z <- as.matrix(z)
    lower <- as.matrix(lower)
    upper <- as.matrix(upper)
    x <- z>=lower & z<=upper
    u <- x[!is.na(x)]
    round(sum(u)/length(u) * 100, 4)
  }else{NA}
}
FRAC2 <- function(z, zhat) {
  z <- as.matrix(z)
  zhat <- as.matrix(zhat)
  x <- z/zhat>=0.5 & z/zhat<=2
  u <- x[!is.na(x)]
  round(sum(u)/length(u) * 100, 4)
}

PMCC <- function(z, zhat, penalty=NULL) {
  if(!is.null(penalty)){
    z <- as.matrix(z)
    zhat <- as.matrix(zhat)
    x <- c(z-zhat)^2
    gof <- sum(x[!is.na(x)])
    penalty <- sum(penalty[!is.na(x)])
    pmcc <- gof + penalty
    round(pmcc, 4)
  }else{NA}
}

my.validation = function (z, zhat, lower=NULL, upper=NULL, penalty=NULL, names = FALSE, pred=FALSE){
  if (names == TRUE){
    cat("##\n Mean Squared Error (MSE) \n Mean Absolute Error (MAE) \n Mean Absolute Percentage Error (MAPE) \n Bias (BIAS) \n Percentage Bias (pBIAS) \n Correlation (CORR) \n Coverage (COV) \n Fraction within a factor of 2 (FRAC2) \n Predictive model choice criterion (PMCC) \n##\n")
  }
  out <- NULL
  out$MSE <- MSE(c(z), c(zhat))
  out$MAE <- MAE(c(z), c(zhat))
  out$MAPE <- MAPE(c(z), c(zhat))
  out$BIAS <- BIAS(c(z), c(zhat))
  out$pBIAS <- pBIAS(c(z), c(zhat))
  out$CORR <- CORR(c(z), c(zhat))
  out$COV <- COV(c(z), c(lower), c(upper))
  out$FRAC2 <- FRAC2(c(z), c(zhat))
  out$PMCC <- PMCC(c(z), c(zhat), c(penalty))
  return(unlist(out))
}


# this function gets in input rows of sample.fitted, each with the sample.size samples and the observed value in last position,
# and the vector of varicances of length sample.size;
# extracts n.pred predictions for each sample generating vector of n.pred x sample.size elements;
# returns predicted values.
# NOTE: correlation cannot be computed from single observed value as SD is zero; the sample mean is returned instead, to allow computation later.
# NOTE: the CI for coverage is computed at each sample.
extract.predicted = function(sample, variance){
  prediction.vector.tot = NULL
  coverage=NULL
  ci_amplitude=NULL
  for(j in 1:(length(sample)-1)){
    # generate 50 samples to get a distribution from which we calculate the quantiles to get the coverage
    prediction.vector = rnorm(n.pred, mean=sample[j], sd=sqrt(variance[j]))
    prediction.vector.tot = c(prediction.vector.tot, prediction.vector) # here we save all samples together: the mean is pol_pred
    ci = quantile(prediction.vector, probs=c(0.025,0.975)) # CI is computed for each sample
    ci_amplitude = c(ci_amplitude, ci[2]-ci[1]) # we save all the CI amplitudes to take the mean
    coverage = c(coverage, COV(sample[length(sample)], lower = ci[1], upper=ci[2]))
  }
  predicted = c(sample[length(sample)], mean(prediction.vector.tot), round(sum(coverage[!is.na(coverage)])/length(coverage[!is.na(coverage)]), 4), mean(ci_amplitude), var(prediction.vector.tot))
  names(predicted)=c(paste0(pol,"_obs"),paste0(pol,"_pred"),"coverage", "ci_amplitude", "pmcc_penalty")
  return(predicted)
}

extract.contents = function(sample){
  print(paste0("sample n. ",sample))
  ps=posterior.samples[[sample]]$latent
  mod.comp = list()
  for(i in 3:nrow(contents)){ # first 2 contents are predictors
    if(contents[i,"length"]==mesh$n){
      mod.comp[[i-2]] = drop(A_pred%*%ps[c(contents[i,"start"]:contents[i,"end"])])
    }else{ 
      mod.comp[[i-2]] = ps[c(contents[i,"start"]:contents[i,"end"])]
    } 
  }
  names(mod.comp) = contents[3:nrow(contents), "tag"]
  
  return(mod.comp)
}

compute.daily.predictions = function(day){
  print(paste0("Calculating predictions for day ",day))
  
  ## prepare first columns with geographic information
  pred_day = pred.grid[,c("easting","northing","sitetype")]
  
  pred_day_data=matrix(NA, nrow=nrow(pred.grid), ncol=nrow(contents))
  colnames(pred_day_data) = c(as.character(contents[3:nrow(contents), "tag"]),"date.idx.no2.urb","date.idx.no2.rks")
  
  for(sample in 1:n.samples){  # for each sample...
    if(sample %% 10 == 0) print(paste0("Sample: ",sample))
    ## ...create a dataframe with all the model components, for that day
    for(i in 1:(nrow(contents)-2)){
      if(contents[i+2,"length"]==mesh$n){
        pred_day_data[,i] = fields[[sample]][i][[1]]
      }else if(contents[i+2,"length"]==1){
        if(as.character(contents[i+2,"tag"])==paste0("aqum_log_",pol)){
          pred_day_data[,i] = fields[[sample]][i][[1]]*pred.grid.aqum[,day]
        }else if(as.character(contents[i+2,"tag"])==paste0("pcm_log_",pol)){
          pred_day_data[,i] = fields[[sample]][i][[1]]*pred.grid.pcm[,date[date$date.idx.no2==day,"year"]-2006]
        }else{
          pred_day_data[,i] = rep(fields[[sample]][i][[1]],nrow(pred_day))
        }
      }else if(contents[i+2,"length"]==n.days){
        pred_day_data[,i] = rep(fields[[sample]][i][[1]][day],nrow(pred_day))
      }else if(contents[i+2,"length"]==n.days*3){
        pred_day_data[,i] = rep(fields[[sample]][i][[1]][day],nrow(pred_day))
        pred_day_data[,nrow(contents)-1] = rep(fields[[sample]][i][[1]][n.days+day],nrow(pred_day))
        pred_day_data[,nrow(contents)] = rep(fields[[sample]][i][[1]][n.days+n.days+day],nrow(pred_day))
      }
    }
    ## then sum all the model components except "intercept.aqum" and "intercept.pcm" and add the column to pred_day
    if("intercept.aqum" %in% colnames(pred_day_data) & !("intercept.pcm" %in% colnames(pred_day_data))){
      pred_day = cbind(pred_day, rowSums(pred_day_data[, -which(colnames(pred_day_data) %in% c("intercept.aqum"))], na.rm = T))
    }else if(!("intercept.aqum" %in% colnames(pred_day_data)) & ("intercept.pcm" %in% colnames(pred_day_data))){
      pred_day = cbind(pred_day, rowSums(pred_day_data[,-which(colnames(pred_day_data) %in% c("intercept.pcm"))], na.rm = T))
    }else if("intercept.aqum" %in% colnames(pred_day_data) & "intercept.pcm" %in% colnames(pred_day_data)){
      pred_day = cbind(pred_day, rowSums(pred_day_data[,-which(colnames(pred_day_data) %in% c("intercept.aqum","intercept.pcm"))], na.rm = T))
    }else{
      pred_day = cbind(pred_day, rowSums(pred_day_data, na.rm = T))
    }
    ## naming the last column with the sample number
    names(pred_day)[ncol(pred_day)]=paste0("sample_",sample)
  }
  saveRDS(pred_day, paste0("predictions_by_day/predictions_",date[day,"date"],".rds"))
}


############################################################################
########                     data preparation                      #########
############################################################################


pol = "no2"

load("workspace_data.RData")

valid = final_dataset[final_dataset$code %in% monitors_val_cheat[[data_id]] , ]
estim = final_dataset[!(final_dataset$code %in% monitors_val_cheat[[data_id]]) , ]

n_monitors=nrow(coordinates.y)
n_data=nrow(estim)+nrow(valid)
n_days=n_data/n_monitors


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


#=================================================== 
### Formulas 
#=================================================== 

formulas = list()

# spatial + temporal effects for aqum
# no pcm
# beta for the copied effects varying with default prior
formulas[[1]] = as.formula(paste0("y ~ -1 + intercept + intercept.aqum + stURB.",pol," + stRKS.",pol," +  
                                  f(csi.field.time, model='rw1', hyper = rw1.aqum.prior, scale.model = TRUE, constr = TRUE) + 
                                  f(csi.field.time.copy, copy='csi.field.time', fixed = FALSE) + 
                                  f(csi.field.space, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) +  
                                  f(csi.field.space.copy, copy='csi.field.space', fixed = FALSE) + 
                                  f(date.idx.",pol,", model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)"))

# spatial effect for pcm
# no aqum
# beta for the copied effects varying with default prior
formulas[[2]] = as.formula(paste0("y ~ -1 + intercept + intercept.pcm + stURB.",pol," + stRKS.",pol," + 
                                  f(psi.field, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) + 
                                  f(psi.field.copy, copy='psi.field', fixed = FALSE) + 
                                  f(date.idx.",pol,", model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)"))

# spatial+temporal effect for AQUM
# beta for copied effects varying with strong prior 
formulas[[3]] = as.formula(paste0("y ~ -1 + intercept + intercept.aqum + intercept.pcm + stURB.",pol," + stRKS.",pol," +  
                                  f(csi.field, model='rw1', hyper = rw1.aqum.prior, scale.model = TRUE, constr = TRUE) +
                                  f(csi.field.copy, copy='csi.field', fixed = FALSE, hyper=lambda.aqum.time.prior) + 
                                  f(psi.field, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) +
                                  f(psi.field.aqum.copy, copy='psi.field', fixed = FALSE, hyper=lambda.aqum.space.prior) + 
                                  f(psi.field.copy, copy='psi.field', fixed = FALSE, hyper=lambda.pcm.prior) + 
                                  f(date.idx.",pol,", model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)"))

# temporal effect for AQUM
# beta for copied effects varying with strong prior 
formulas[[4]] = as.formula(paste0("y ~ -1 + intercept + intercept.aqum + intercept.pcm + stURB.",pol," + stRKS.",pol," +  
                                  f(csi.field, model='rw1', hyper = rw1.aqum.prior, scale.model = TRUE, constr = TRUE) +
                                  f(csi.field.copy, copy='csi.field', fixed = FALSE, hyper=lambda.aqum.time.prior) + 
                                  f(psi.field, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) +
                                  f(psi.field.copy, copy='psi.field', fixed = FALSE, hyper=lambda.pcm.prior) + 
                                  f(date.idx.",pol,", model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)"))


# linear effects for AQUM and PCM  
# residual structured spatial effect and temporal effect
formulas[[5]] = as.formula(paste0("y ~ -1 + intercept + stURB.",pol," + stRKS.",pol," +  aqum_log_",pol," + pcm_log_",pol," + 
                                  f(csi.field, model='rw1', hyper = rw1.aqum.prior, scale.model = TRUE, constr = TRUE) +
                                  f(psi.field, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) +
                                  f(date.idx.",pol,", model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)"))

# linear effects for AQUM and PCM plugged in from separate models (kriging)
# residual structured spatial effect and temporal effect
formulas[[6]] = as.formula(paste0("y ~ -1 + intercept + stURB.",pol," + stRKS.",pol," +  aqum_pred + pcm_pred + 
                                  f(csi.field, model='rw1', hyper = rw1.aqum.prior, scale.model = TRUE, constr = TRUE) +
                                  f(psi.field, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) +
                                  f(date.idx.",pol,", model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)"))

# Model from Mukhopadhyay,2017
# aqum only 
# non stationary 
formulas[[7]] = as.formula(paste0(pol," ~ aqum_log_",pol," + aqum_log_URB + stURB.",pol," + aqum_log_RKS + stRKS.",pol))

# Model from Mukhopadhyay,2017
# aqum + pcm
# non stationary
formulas[[8]] = as.formula(paste0(pol," ~ aqum_log_",pol," + pcm_log_",pol," + aqum_log_URB + pcm_log_URB + stURB.",pol," + aqum_log_RKS + pcm_log_RKS + stRKS.",pol))

# Model from Mukhopadhyay,2017
# aqum + pcm
# stationary
formulas[[9]] = as.formula(paste0(pol," ~ aqum_log_",pol," + pcm_log_",pol," + aqum_log_URB + pcm_log_URB + stURB.",pol," + aqum_log_RKS + pcm_log_RKS + stRKS.",pol))


############################################################################
########        run scripts for each dataset and formula           #########
############################################################################

##--- run models and extract results
for(data_id in 1:length(monitors_val)){
  for(formula_id in 1:length(formulas)){
    source("01_models.R")
    source("02_results.R")
  }
}

