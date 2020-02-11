remove(list = ls())

load("workspace_results.RData")

set.seed(123)

wd = getwd()


print("#############################################################################################")
print(paste0("Extracting results for model ",formula_id, ", dataset ",data_id))
print("#############################################################################################")

print(formulas[[formula_id]])

##--- Create empty list to save all results
results_mod=list()

if(formula_id<=6){
  
  coordinates.estim<-unique(estim[,c("loc.idx","easting","northing")])
  coordinates.valid<-unique(valid[,c("loc.idx","easting","northing")])

  results_mod$mod_size=object.size(mod)
  results_mod$exec_time = mod$cpu.used
  
  
  ######################################################
  print("#---  Summary of Posterior distributions of       ---#")
  print("#--- parameters and hyperparameters of interest   ---#")
  ######################################################
  
  png('betas.png', width = 15, height = 15, unit="in", res=600)
  par(mfrow=c(4,5))
  
  
  ##--- FIXED EFFECTS
  
  results_mod$summary.fixed = round(mod$summary.fixed,4)
  for(i in 1:nrow(mod$summary.fixed)){
    fixed = rownames(mod$summary.fixed)[i]
    plot(mod$marginals.fixed[[i]], type="l", main=paste0(fixed," (mean=",round(mod$summary.fixed[i,1],3),")"))
    abline(v=mod$summary.fixed[i,1], col="red")
  }
  
  ##--- HYPERPARAMETERS
  
  results_mod$summary.hyperpar = round(mod$summary.hyperpar,4)
  for(i in 1:nrow(mod$summary.hyperpar)){
    hyperpar = rownames(mod$summary.hyperpar)[i]
    plot(mod$marginals.hyperpar[[i]], type="l", main=paste0(hyperpar," (mean=",round(mod$summary.hyperpar[i,1],3),")"))
    abline(v=mod$summary.hyperpar[i,1], col="red")
  }
  
  print(results_mod$summary.fixed)
  print(results_mod$summary.hyperpar)
  
  
  ######################################################
  print("#--- goodness of fit  ---#")
  ######################################################
  
  # extract summary fitted values
  fitted=mod$summary.fitted.values[inla.stack.index(stack,"val.y")$data,]
  
  validation= list()
  
  validation$data = data.frame(valid[,c("code","site.type",paste0("date.idx.",pol),pol,paste0(pol,"_log"),"easting","northing")],fitted=fitted)
  validation$data$code = as.character(validation$data$code)
  
  validation$dic = mod$dic$dic
  validation$waic = mod$waic$waic
  validation$mean.deviance = mod$dic$mean.deviance
  validation$pD = mod$dic$p.eff
  
  validation$data$res.log=valid[,paste0(pol,"_log")]-fitted$mean
  validation$data$stdres.log = (valid[,paste0(pol,"_log")] - fitted$mean) / fitted$sd
  validation$measures = my.validation(z=valid[,paste0(pol,"_log")], zhat=fitted$mean, names=TRUE)
  
  
  ##--- CHECK MODEL PERFORMANCE
  
  n.failures = sum(mod$cpo$failure, na.rm = T)
  if(mean(mod$cpo$failure, na.rm = T)!=0){
    # summary(mod$cpo$failure)
    #  Two options:
    #  1. recompute using control.inla=list(int.strategy = "grid", diff.logdens = 4, strategy = "laplace", npoints = 21) see http://www.r-inla.org/faq
    #  2. run mod = inla.cpo(mod) #recomputes in an efficient way the cpo/pit for which mod$cpo$failure > 0
    mod_imp = inla.cpo(mod)
    print("Model cpo have been improved")
    validation$logscore = -mean(log(mod_imp$cpo$cpo), na.rm=T)
    hist(mod_imp$cpo$pit, breaks=100, main=paste0("PIT improved - n.failures=", n.failures))
    hist(mod_imp$cpo$cpo, breaks=100, main=paste0("CPO improved - n.failures=", n.failures))
  }else{
    validation$logscore = -mean(log(mod$cpo$cpo), na.rm=T)
    hist(mod$cpo$pit, breaks=100, main=paste0("PIT - n.failures=", n.failures))
    hist(mod$cpo$cpo, breaks=100, main=paste0("CPO - n.failures=", n.failures))
  }
  
  
  dev.off()
  
  results_mod$validation = validation
  
  
  
  ######################################################
  print("#--- Correlation observed vs expected (on log scale) ---#")
  ######################################################
  
  png('correlation_val.png', width = 10, height = 6, unit="in", res=600)
  par(mfrow=c(1,2))
  plot(y=validation$data[,paste0(pol,"_log")],x=validation$data$fitted.mean,asp=1,xlab='Expected values', ylab='Observed values',
       main=paste("Observed vs Expected \ncor = ",validation$measures[6],sep=""))
  abline(0:1, col="red")
  plot(y=(validation$data[,paste0(pol,"_log")]/validation$data$fitted.mean),x=validation$data$fitted.mean,
       xlab='Expected values', ylab='Observed/Expected',
       main="Observed/Expected vs Expected")
  abline(h=1, col="red")
  dev.off()
  
  png('correlation_val_bysite.png', width = 15, height = 15, unit="in", res=600)
  par(mfrow=c(5,5))
  unique.site = unique(validation$data$code)[1:25] # ifelse(length(unique(validation$data$code))>25, unique(validation$data$code)[1:25], unique(validation$data$code))
  for(site in unique.site){
    index=which(validation$data$code==site)
    plot(y=validation$data[,paste0(pol,"_log")][index],x=validation$data$fitted.mean[index],asp=1,xlab='Expected values', ylab='Observed values',
         main=paste(site,": cor = ",CORR(validation$data[,paste0(pol,"_log")][index],validation$data$fitted.mean[index]),sep=""))
    abline(0:1, col="red")
  }
  dev.off()
  
  png('correlation_val_bysitetype.png', width = 15, height = 5, unit="in", res=600)
  par(mfrow=c(1,3))
  for(sitetype in unique(valid$site.type)){
    index=which(valid$site.type==sitetype)
    plot(y=validation$data[,paste0(pol,"_log")][index],x=validation$data$fitted.mean[index],asp=1,xlab='Expected values', ylab='Observed values',
         main=paste(sitetype,": cor = ",CORR(validation$data[,paste0(pol,"_log")][index],validation$data$fitted.mean[index]),sep=""))
    abline(0:1, col="red")
  }
  dev.off()
  
  month=as.numeric(substr(valid$date,6,7))
  png('correlation_val_bymonth.png', width = 15, height = 15, unit="in", res=600)
  par(mfrow=c(3,4))
  for(m in 1:12){
    index=which(month==m)
    plot(y=validation$data[,paste0(pol,"_log")][index],x=validation$data$fitted.mean[index],asp=1,xlab='Expected values', ylab='Observed values',
         main=paste("Month ",m,": cor = ",CORR(validation$data[,paste0(pol,"_log")][index],validation$data$fitted.mean[index]),sep=""))
    abline(0:1, col="red")
  }
  dev.off()
  
  
  ######################################################
  print("#--- check spatial and temporal trend of residuals  ---#")
  ######################################################
  
  print("# residual temporal trend")
  
  z <- ts(matrix(validation$data$res.log, ncol=nrow(coordinates.valid),byrow=T), start = c(2007, 1), frequency = 365)
  
  png('temp_res_fix.png', width = 15, height = 10, unit="in", res=600)
  par(mfrow=c(5,5), mar=c(3,2,6,2), oma=c(0,0,1,0))
  for(i in 1:25){ 
    plot(z[,i], main=paste0(unique(valid$code)[i]," - ",valid[valid$code==unique(valid$code)[i],"site.type"][1]), ylab="", xlab="", ylim=c(min(z, na.rm = T),max(z,na.rm = T)))
    abline(h=0, col="red")
  }
  title("Temporal trend of residuals on log scale, by monitor", outer=TRUE,cex.main=2)
  dev.off()
  
  png('temp_res.png', width = 15, height = 10, unit="in", res=600)
  par(mfrow=c(5,5), mar=c(3,2,6,2), oma=c(0,0,1,0))
  for(i in 1:25){ 
    plot(z[,i], main=paste0(unique(valid$code[i])," - ",valid[valid$code==unique(valid$code)[i],"site.type"][1]), ylab="", xlab="")
    abline(h=0, col="red")
  }
  title("Temporal trend of residuals on log scale, by monitor", outer=TRUE,cex.main=2)
  dev.off()
  

  print("variogram")

  # #Create a SpatialPointsDataFrame in UTM coordinates
  res=validation$data
  coordinates(res) <- ~easting+northing
  proj4string(res) <- crs(london.shape)
  # Create STFDF objest (spatial, temporal and data components)
  res.SP <- SpatialPoints(unique(res@coords),crs(london.shape))
  res.TM <- seq(from=as.POSIXct(strptime("2007-01-01", "%Y-%m-%d")), to=as.POSIXct(strptime("2011-12-31", "%Y-%m-%d")), by="day")
  res.DF <- data.frame(res=res@data$res)
  res.STFDF <- STFDF(sp=res.SP,time=res.TM,data=res.DF) #data with spatial index moving faster
  # Variogram
  var <- gstat::variogramST(res~1, data=res.STFDF, assumeRegular=T)
  var.plot1=plot(var)
  var.plot2=plot(var, map=F)
  var.plot3=plot(var, wireframe=T)
  png("variogram.png", width =10, height = 7, unit="in", res=600)
  gridExtra::grid.arrange(var.plot1, var.plot2, var.plot3, ncol=2)
  dev.off()
  
  print("# residual spatial patterns")
  
  spat_res = ggplot(validation$data[!is.na(validation$data$res.log),]) +
    geom_boxplot(aes(x=code,y=res.log, fill = site.type), outlier.size = .5) +
    labs(x="Monitors",y="Residuals (log scale)") +
    theme(axis.title = element_text(size=8),axis.text = element_text(size=6), legend.title = element_text(size=8))
  ggsave('spat_res.png', plot=spat_res, width = 10, height = 3, units="in", dpi=600, scale=3)
  
  
  
  ######################################################
  print("#--- Extracting posterior means of latent fields ---#")
  ######################################################
  
  print(names(mod$summary.random))
  index.temp.field = which(substr(names(mod$summary.random),1,9) == "csi.field" & substr(names(mod$summary.random),1,15) != "csi.field.space" )
  index.spat.field = which(substr(names(mod$summary.random),1,9) == "psi.field" | substr(names(mod$summary.random),1,15) == "csi.field.space")
  print(index.temp.field)
  print(index.spat.field)
  
  print("# time-sitetype interaction")
  
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
  
  print("# latent temporal trend ")
  
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
  

  print("# latent spatial field")
  
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
    
    ggsave(paste0(names(mod$marginals.random)[i],'.png'), 
           plot = ggplot(spat.field) +
             geom_raster(aes(x, y, fill = mean.log)) +
             #gg(mesh.pcm) +  # this requires package inlabru and plots mesh over the mean.log layer
             scale_fill_viridis(expression(paste(log(NO[2])," (",mu,"g/",m^3,")"))) +
             ggtitle(paste0(names(mod$marginals.random)[i]," - posterior mean")) +
             geom_path(data = fortify(shape), aes(group = group, x = long, y = lat)) +
             geom_path(data = fortify(london.shape), aes(group = group, x = long, y = lat), size=0.5) +
             theme(plot.title = element_text(family = "sans", color="#666666", size=14, hjust=0.5, face="bold"),
                   axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
                   legend.key.height = unit(1.5,"cm"), 
                   legend.text = element_text(size=7, vjust=-0.5)) +
             labs(x="Easting",y="Northing"), 
      width = 8, height = 8, units="in", dpi=600)
    
  }

  ######################################################
  print("#--- predictive capability  ---#")
  ######################################################
  
  sample.size=50 
  n.pred=50
  pred.perf = list()
  predictions= list()

  # extract sample.size values from each marginal of the fitted values
  # posterior.samples.all = inla.posterior.sample(n = sample.size.all, result=mod, intern = FALSE, use.improved.mean = TRUE, add.names = TRUE, seed = 0L)
  # posterior.samples = posterior.samples.all[1:sample.size]
  posterior.samples = inla.posterior.sample(n = sample.size, result=mod, intern = FALSE, use.improved.mean = TRUE, add.names = TRUE, seed = 0L)
  print("Posterior samples extracted.")
  saveRDS(posterior.samples, "posterior_samples.rds")
  print("Posterior samples saved")

  samples.fitted = unlist(lapply(lapply(posterior.samples,"[[", "latent" ),"[", inla.stack.index(stack,"val.y")$data))

  #reshape in order to have a matrix of nrow(valid) x sample.size
  samples.fitted = matrix(samples.fitted, nrow=nrow(valid), byrow = F)

  # extract sample.size values from marginal of the likelihood variance
  # (it is always the last of the Gaussian observations, so we extract the position)
  par.index = length(which(substr(rownames(mod$summary.hyperpar),1,26) == "Precision for the Gaussian"))
  samples.var = 1/unlist(lapply(lapply(posterior.samples,"[[", "hyperpar" ),"[", par.index))
  # this is not exactly what we want, in theory we should sample from the transformed posterior marginal
  # rather than inverting the values sampled from the posterior marginal of the precision

  samples.fitted = cbind(samples.fitted, obs.value=valid[,paste0(pol,"_log")])

  predictive.capability=apply(samples.fitted, MARGIN = 1, FUN=extract.predicted, variance=samples.var)

  predictive.capability = as.data.frame(t(predictive.capability))
  predictive.capability = cbind(code=valid$code, date.idx=valid$date.idx, site.type=valid$site.type, predictive.capability)

  predictive.capability$CRPS = crps_sample(y=samples.fitted[,(ncol(samples.fitted)-1)], dat=as.matrix(samples.fitted[,1:(ncol(samples.fitted)-1)]))

  results_mod$predictive_capability = predictive.capability

  results_mod$samples.fitted = samples.fitted

  saveRDS(results_mod, "results.rds")


  ######################################################
  print("#--- extract daily predictions ---#")
  ######################################################

  n.samples = 50
  n.days = 1826
  n.locs = nrow(pred.grid)

  A_pred =inla.spde.make.A(mesh, loc=cbind(pred.grid$easting, pred.grid$northing))
  
  print(mod$misc$configs$contents)
  contents=as.data.frame(mod$misc$configs$contents)
  contents$end = contents$start + contents$length - 1

  fields = lapply(1:n.samples, FUN=extract.contents)
  
  saveRDS(fields, "posterior_samples_fields.rds")

  # Create folder to save daily predictions
  if (!file.exists(file.path(getwd(), "predictions_by_day"))){
    dir.create(file.path(getwd(), "predictions_by_day"))
  }

  date = data.frame(date=seq.Date(as.Date("2007-01-01"), as.Date("2011-12-31"), "days"),
                    date.idx.no2=c(1:n_days),
                    year=as.numeric(substr(unique(aqum$date),1,4)),
                    month=as.numeric(substr(unique(aqum$date),6,7)),
                    day=as.numeric(substr(unique(aqum$date),9,10)))

  # This can be edited to run in parallel with mclapply() if machine supports multi-core 
  lapply(1:n.days, FUN=compute.daily.predictions)


  ###############################################################################
  print("#---  compute monthly and annual means from daily prediction files ---#")
  ###############################################################################
  
  setwd("predictions_by_day")
  
  monthly.mean=pred.grid[,1:3]
  annual.mean=pred.grid[,1:3]
  for(year in c("2007":"2011")){
    y=as.character(year)
    for(month in c("-01":"-12")){
      m=ifelse(nchar(as.character(month))==3, as.character(month), paste0("-0",substr(as.character(month),2,2)))
      file.list = list.files(pattern = paste0("predictions_",y,m,"-"), full.names = TRUE)
      print(paste0("Calculating prediction mean for ",y,m))
      pred.mean=rep(0,nrow(pred.grid))
      for(file in file.list){
        pred=readRDS(file)
        pred.mean=pred.mean+rowMeans(pred[,which(substr(colnames(pred),1,7)=="sample_")])
      }
      pred.mean=pred.mean/length(file.list)
      monthly.mean=cbind(monthly.mean,pred.mean)
      names(monthly.mean)[ncol(monthly.mean)]=paste0("no2_",y,m)
    }
    pred.mean=rowMeans(monthly.mean[,(ncol(monthly.mean)-11):ncol(monthly.mean)])
    annual.mean=cbind(annual.mean,pred.mean)
    names(annual.mean)[ncol(annual.mean)]=paste0("no2_",y)
  }

  print(head(monthly.mean))
  print(head(annual.mean))

  setwd("..")
  saveRDS(monthly.mean, "pred_mean_monthly.rds")
  saveRDS(annual.mean, "pred_mean_annual.rds")
  print("Saved annual and monthly average of predictions")
  
  ###############################################################################
  print("#---  plot monthly and annual means  ---#")
  ###############################################################################
  
  ##--- reshape to long format
  
  annual.mean.long = melt(annual.mean, id.vars = c("easting","northing","sitetype"), variable.name = "year", value.name = pol)
  annual.mean.long$year = as.factor(substr(as.character(annual.mean.long$year),5,8))
  
  monthly.mean.long = melt(monthly.mean, id.vars = c("easting","northing","sitetype"), variable.name = "month", value.name = pol)
  monthly.mean.long$year = as.factor(substr(as.character(monthly.mean.long$month),5,8))
  monthly.mean.long$yearmonth = as.factor(substr(as.character(monthly.mean.long$month),5,11))
  monthly.mean.long$month = as.factor(substr(as.character(monthly.mean.long$month),10,11))
  
  ggsave('annual.mean.png', 
         plot=ggplot(annual.mean.long) +
           geom_raster(aes(x=easting,y=northing, fill = no2)) +
           facet_wrap("year",ncol=3) +
           scale_fill_viridis(expression(paste(log(NO[2])," (",mu,"g/",m^3,")"))) +
           ggtitle(expression(paste("Average ",log(NO[2])," predictions by year"))) +
           geom_path(data = fortify(shape), aes(group = group, x = long, y = lat))+
           geom_path(data = fortify(london.shape), aes(group = group, x = long, y = lat), size=0.5)+
           theme(plot.title = element_text(family = "sans", color="#666666", size=14, hjust=0.5, face="bold"),
                 axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
                 legend.key.height = unit(1.5,"cm"), #legend.title = element_text(margin = 10), 
                 legend.text = element_text(size=7, vjust=-0.5)) +
           labs(x="Easting",y="Northing"), 
         width = 9, height = 6, unit="in", dpi=600)
  
  
  ggsave('monthly.mean.png', 
         plot=ggplot(monthly.mean.long) + 
           geom_raster(aes(x=easting,y=northing, fill = no2))+
           facet_wrap("yearmonth", ncol=12) +
           scale_fill_viridis(expression(paste(log(NO[2])," (",mu,"g/",m^3,")"))) +
           ggtitle(expression(paste("Average ",log(NO[2])," predictions by year and month"))) +
           geom_path(data = fortify(shape), aes(group = group, x = long, y = lat))+
           geom_path(data = fortify(london.shape), aes(group = group, x = long, y = lat), size=0.5)+
           theme(plot.title = element_text(family = "sans", color="#666666", size=14, hjust=0.5, face="bold"),
                 axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
                 legend.key.height = unit(1.5,"cm"), #legend.title = element_text(margin = 10), 
                 legend.text = element_text(size=7, vjust=-0.5)) +
           labs(x="Easting",y="Northing"), 
         width = 15, height = 10, unit="in", dpi=600)
  
  
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
           scale_fill_viridis(expression(paste(log(NO[2])," (",mu,"g/",m^3,")"))) +
           ggtitle("Daily predictions on selected days with air pollution events (top row) or low concentration (bottom row)") +
           geom_path(data = fortify(london.shape), aes(group = group, x = long, y = lat), size=0.5)+
           geom_path(data = fortify(shape), aes(group = group, x = long, y = lat))+  
           geom_path(data = fortify(roads_major), aes(group = group, x = long, y = lat), size=0.3, color="grey85")+
           theme(plot.title = element_text(family = "sans", color="#666666", size=14, hjust=0.5, face="bold"),
                 axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
                 legend.key.height = unit(1,"cm"),
                 legend.text = element_text(size=7, vjust=-0.5)),
         width=15, height=8, unit="in", dpi=600)
  

  } else if(formula_id>=7){

  ######################################################
  print("#---  Summary of Posterior distributions of       ---#")
  print("#--- parameters and hyperparameters of interest   ---#")
  ######################################################

  results_mod$summary.parameters = mod$parameters

  marginals.parameters = read.table("OutGPP_Values_Parameter.txt", stringsAsFactors=F, header=F)
  colnames(marginals.parameters) = rownames(mod$parameters)

  png('betas.png', width = 15, height = 15, unit="in", res=600)
  par(mfrow=c(4,5))

  for(i in 1:ncol(marginals.parameters)){
    plot(density(marginals.parameters[,i]),
         main=paste0(colnames(marginals.parameters)[i]," (mean=",round(results_mod$summary.parameters$Mean[i],3),")"))
    abline(v=results_mod$summary.parameters$Mean[i], col="red")
  }

  dev.off()

  ######################################################
  print("#--- Correlation observed vs expected  ---#")
  ######################################################

  validation= list()

  validation$data = data.frame(valid[,c("code","site.type",paste0("date.idx.",pol),pol,paste0(pol,"_log"),"easting","northing")],fitted=mod$prediction)
  validation$data$code = as.character(validation$data$code)

  validation$data$res.log=valid[,pol]-mod$prediction$Mean
  validation$measures = my.validation(z=valid[,pol], zhat=mod$prediction$Mean, names=TRUE)

  results_mod$validation = validation

  png('correlation_val.png', width = 10, height = 6, unit="in", res=600)
  par(mfrow=c(1,2))
  plot(y=validation$data[,pol],x=validation$data$fitted.Mean,asp=1,xlab='Expected values', ylab='Observed values',
       main=paste("Observed vs Expected \ncor = ",validation$measures[6],sep=""))
  abline(0:1, col="red")
  plot(y=(validation$data[,pol]/validation$data$fitted.Mean),x=validation$data$fitted.Mean,
       xlab='Expected values', ylab='Observed/Expected',
       main="Observed/Expected vs Expected")
  abline(h=1, col="red")
  dev.off()

  png('correlation_val_bysite.png', width = 15, height = 15, unit="in", res=600)
  par(mfrow=c(5,5))
  unique.site = unique(validation$data$code)#[1:25] # ifelse(length(unique(validation$data$code))>25, unique(validation$data$code)[1:25], unique(validation$data$code))
  for(site in unique.site){
    index=which(validation$data$code==site)
    plot(y=validation$data[,pol][index],x=validation$data$fitted.Mean[index],asp=1,xlab='Expected values', ylab='Observed values',
         main=paste(site,": cor = ",CORR(validation$data[,pol][index],validation$data$fitted.Mean[index]),sep=""))
    abline(0:1, col="red")
  }
  dev.off()

  png('correlation_val_bysitetype.png', width = 15, height = 5, unit="in", res=600)
  par(mfrow=c(1,3))
  for(sitetype in unique(valid$site.type)){
    index=which(valid$site.type==sitetype)
    plot(y=validation$data[,pol][index],x=validation$data$fitted.Mean[index],asp=1,xlab='Expected values', ylab='Observed values',
         main=paste(sitetype,": cor = ",CORR(validation$data[,pol][index],validation$data$fitted.Mean[index]),sep=""))
    abline(0:1, col="red")
  }
  dev.off()

  month=as.numeric(substr(valid$date,6,7))
  png('correlation_val_bymonth.png', width = 15, height = 15, unit="in", res=600)
  par(mfrow=c(3,4))
  for(m in 1:12){
    index=which(month==m)
    plot(y=validation$data[,pol][index],x=validation$data$fitted.Mean[index],asp=1,xlab='Expected values', ylab='Observed values',
         main=paste("Month ",m,": cor = ",CORR(validation$data[,pol][index],validation$data$fitted.Mean[index]),sep=""))
    abline(0:1, col="red")
  }
  dev.off()

  
  ######################################################
  print("#--- predictive capability  ---#")
  ######################################################

  n.pred = 50
  n.samples=60

  pred.mcmc = read.table("Prediction_sites_mcmc1.txt", stringsAsFactors=F, header=F)

  print(mod$parameters)

  # we need samples of sigma2eps: instead of exctracting samples from its posterior marginal,
  # we use rnorm because we do not have the posterior marginal
  samples.var = rnorm(n.samples, mean=mod$parameters$Mean[which(rownames(mod$parameters)=="sig2eps")], 
                      sd=mod$parameters$SD[which(rownames(mod$parameters)=="sig2eps")])

  pred.mcmc = pred.mcmc[sample(nrow(pred.mcmc), size=n.samples),]

  pred.mcmc = rbind(pred.mcmc,valid$no2)


  predictive.capability=apply(pred.mcmc, MARGIN = 2, FUN=extract.validation, variance=samples.var)

  predictive.capability = as.data.frame(t(predictive.capability))
  predictive.capability = cbind(code=valid$code,
                                date.idx=valid$date.idx,
                                site.type=valid$site.type,
                                predictive.capability)

  pred.mcmc = as.data.frame(t(pred.mcmc))
  predictive.capability$CRPS = crps_sample(y=pred.mcmc[,(ncol(pred.mcmc)-1)], dat=as.matrix(pred.mcmc[,1:(ncol(pred.mcmc)-1)]))

  results_mod$predictive.capability=predictive.capability
  saveRDS(results, "results.rds")

  
  ######################################################
  print("#--- extract daily predictions ---#")
  ######################################################

  coordinates.pred = pred.grid[,c("longitude","latitude")]
  
  pred.day = predict.spT(mod,newdata = pred.grid, newcoords = coordinates.pred[,c("longitude","latitude")])
  
  ###############################################################################
  print("#---  plot monthly and annual means  ---#")
  ###############################################################################

}

setwd("..")
q(save = "no")
