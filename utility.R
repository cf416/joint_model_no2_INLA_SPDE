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
COV <- function(z, lower=NULL, upper=NULL, coverage=NULL) {
  if(!is.null(lower) && !is.null(upper)){
    z <- as.matrix(z)
    lower <- as.matrix(lower)
    upper <- as.matrix(upper)
    x <- z>=lower & z<=upper
    u <- x[!is.na(x)]
    round(sum(u)/length(u) * 100, 4)
  }else if(!is.null(coverage)){
    round(mean(coverage, na.rm = T),4)
  }
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

my.validation = function (z, zhat, lower=NULL, upper=NULL, penalty=NULL, coverage=NULL, names = FALSE, pred=FALSE){
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
  out$COV <- COV(c(z), c(lower), c(upper), c(coverage))
  out$FRAC2 <- FRAC2(c(z), c(zhat))
  out$PMCC <- PMCC(c(z), c(zhat), c(penalty))
  return(unlist(out))
}


# this function gets in input rows of samples.fitted, each with the sample.size samples and the observed value in last position,
# and the vector of varicances of length sample.size;
# extracts n.pred predictions for each sample generating vector of n.pred x sample.size elements;
# returns predicted values.
# NOTE: correlation cannot be computed from single observed value as SD is zero; the sample mean is returned instead, to allow for computation later.
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

# This function takes in input the inla.posterior.sample id and returns the list of sampled model components 
# posterior.samples and contents must be in the global environment.
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

# This function takes in input the index of the day and saves a file "predictions_by_day/predictions_yyyy_mm_dd.rds" containing 
# a data.frame with n.samples columns, one for each sample, and a number of rows equal to the number of locations in the prediction grid.
# n.sample and pred.grid must be in the global environment.
compute.daily.predictions = function(day){
  print(paste0("Calculating predictions for day ",day))
  
  ## prepare first columns with geographic information
  pred_day = pred.grid[,c("easting","northing","sitetype")]
  
  pred_day_data=matrix(NA, nrow=nrow(pred.grid), ncol=nrow(contents)-2)
  colnames(pred_day_data) = as.character(contents[3:nrow(contents), "tag"])
  
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
        pred_day_data[,i] = rep(fields[[sample]][i][[1]][day],nrow(pred_day))*pred.grid$pred.stRUR +
          rep(fields[[sample]][i][[1]][n.days+day],nrow(pred_day))*pred.grid$pred.stURB +
          rep(fields[[sample]][i][[1]][n.days+n.days+day],nrow(pred_day))*pred.grid$pred.stRKS
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
