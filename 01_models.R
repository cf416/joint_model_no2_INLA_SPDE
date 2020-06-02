print("####################################################################################")
print(paste0("Output_",formula_id,"_",pol,"_",data_id))
print("####################################################################################")

if (!file.exists(file.path(getwd(),paste0("Output_",formula_id,"_",pol,"_",data_id)))){
  dir.create(file.path(getwd(),paste0("Output_",formula_id,"_",pol,"_",data_id)))
} 
setwd(paste0("Output_",formula_id,"_",pol,"_",data_id))
print(getwd())

if (!file.exists(file.path(getwd(),"/workspace_results.RData"))){
  
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
  
  lambda.aqum.time.prior = list(theta = list(prior = 'normal', param = c(0.9, 0.01), initial=0.9))  # lambda_2,3
  lambda.aqum.space.prior = list(theta = list(prior = 'normal', param = c(1.1, 0.01), initial=1.1))  # lambda_1,2
  lambda.pcm.prior = list(theta = list(prior = 'normal', param = c(1.3, 0.01), initial=1.3)) # lambda_1,3
  
  
  ##--- NOTE: the INLA call can be parallelized using Pardiso - see inla.pardiso()
  
  ############################################################################
  ########                JOINT MODEL WITH AQUM ONLY                #########
  ############################################################################
  
  if(formula_id==1){
    
    A_aqum <- inla.spde.make.A(mesh=mesh,cbind(aqum$easting, aqum$northing))
    stk_aqum <- inla.stack(data=list(y=cbind(aqum[,paste0(pol,"_log")], NA)),  
                           A=list(1,A_aqum),
                           effects=list(list(intercept.aqum=rep(1,nrow(aqum)),
                                             csi.field.time = aqum$date.idx),
                                        list(csi.field.space = c(1:spde$n.spde))),
                           tag="est.aqum")
    
    ## data stack: include all the effects
    A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
    stk_y_e <- inla.stack(data=list(y=cbind(NA,estim[,paste0(pol,"_log")])),  
                          effects=list(list(csi.field.time.copy=estim[,paste0("date.idx.",pol)]),
                                       list(csi.field.space.copy=c(1:spde$n.spde)),
                                       data.frame(intercept=1,
                                                  estim[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing")])), 
                          A=list(1, A_y_e, 1),
                          tag="est.y")
    
    ### validation scenario
    A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
    stk_y_v <- inla.stack(data=list(y=cbind(NA,rep(NA,length(valid[,paste0(pol,"_log")])))),
                          effects=list(list(csi.field.time.copy=valid[,paste0("date.idx.",pol)]),
                                       list(csi.field.space.copy=c(1:spde$n.spde)),
                                       data.frame(intercept=1,
                                                  valid[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing")])),
                          A=list(1, A_y_v, 1),
                          tag="val.y")
    
    stack <- inla.stack(stk_aqum, 
                        stk_y_v, 
                        stk_y_e) 
    
    mod<- inla(formulas[[formula_id]],
               family=c("gaussian","gaussian"),
               data=inla.stack.data(stack),
               control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
               #control.fixed = list(mean.intercept=3, prec.intercept=0.001),
               control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
               control.inla = list(adaptive,int.strategy='eb'),
               verbose=TRUE)
    
    # saveRDS(mod, "mod_",formula_id,".rds")  
  }
  
  
  
  ############################################################################
  ########                JOINT MODEL WITH PCM ONLY                 #########
  ############################################################################
  
  if(formula_id==2){
    
    
    ## ***** pcm *****
    ### spatial only index
    ### projection matrix: use mesh and pcm data coordinates 
    A_pcm <- inla.spde.make.A(mesh=mesh, cbind(pcm$easting, pcm$northing))
    ### the pcm data stack
    stk_pcm <- inla.stack(data=list(y=cbind(pcm[,paste0(pol,"_log")], NA)),  
                          effects=list(list(intercept.pcm=rep(1,nrow(pcm))), list(psi.field=1:spde$n.spde)),
                          A=list(1,A_pcm),
                          tag="est.pcm")
    
    ## data stack: include all the effects
    A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
    stk_y_e <- inla.stack(data=list(y=cbind(NA,estim[,paste0(pol,"_log")])),  
                          effects=list(list(psi.field.copy=1:spde$n.spde),
                                       data.frame(intercept=1,
                                                  estim[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing")])), 
                          A=list(A_y_e, 1),
                          tag="est.y")
    
    ### validation scenario
    A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
    stk_y_v <- inla.stack(data=list(y=cbind(NA,rep(NA,length(valid[,paste0(pol,"_log")])))),
                          effects=list(list(psi.field.copy=1:spde$n.spde),
                                       data.frame(intercept=1,
                                                  valid[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing")])),
                          A=list(A_y_v, 1),
                          tag="val.y")
    
    stack <- inla.stack(stk_pcm,
                        stk_y_v, 
                        stk_y_e) 
    
    mod<- inla(formulas[[formula_id]],
               family=c("gaussian","gaussian"),
               data=inla.stack.data(stack),
               control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
               # control.fixed = list(mean.intercept=3, prec.intercept=0.001),
               control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
               control.inla = list(adaptive,int.strategy='eb'),
               verbose=TRUE)
    
    # saveRDS(mod, "mod_",formula_id,".rds")  
  }
  
  
  ############################################################################
  ########                JOINT MODEL WITH ALL DATA                  #########
  ############################################################################
  
  
  if(formula_id==3){
    
    
    ## ***** pcm *****
    A_pcm <- inla.spde.make.A(mesh=mesh, cbind(pcm$easting, pcm$northing))
    stk_pcm <- inla.stack(data=list(y=cbind(pcm[,paste0(pol,"_log")], NA, NA)),  
                          effects=list(list(intercept.pcm=rep(1,nrow(pcm))), list(psi.field=1:spde$n.spde)),
                          A=list(1,A_pcm),
                          tag="est.pcm")
    
    ## ***** Aqum *****
    A_aqum <- inla.spde.make.A(mesh=mesh,cbind(aqum$easting, aqum$northing))
    stk_aqum <- inla.stack(data=list(y=cbind(NA, aqum[,paste0(pol,"_log")], NA)),  
                           effects=list(list(intercept.aqum=1, csi.field=aqum$date.idx),
                                        list(psi.field.aqum.copy=1:spde$n.spde)), 
                           A=list(1, A_aqum), 
                           tag="est.aqum")
    
    ## data stack: include all the effects
    A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
    stk_y_e <- inla.stack(data=list(y=cbind(NA,NA,estim[,paste0(pol,"_log")])),  
                          effects=list(list(csi.field.copy=estim[,paste0("date.idx.",pol)]),
                                       list(psi.field.copy=1:spde$n.spde),
                                       data.frame(intercept=1,
                                                  estim[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing")])), 
                          A=list(1, A_y_e, 1),
                          tag="est.y")
    
    ### validation scenario
    A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
    stk_y_v <- inla.stack(data=list(y=cbind(NA,NA,rep(NA,length(valid[,paste0(pol,"_log")])))),
                          effects=list(list(csi.field.copy=valid[,paste0("date.idx.",pol)]),
                                       list(psi.field.copy=1:spde$n.spde),
                                       data.frame(intercept=1,
                                                  valid[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing")])),
                          A=list(1, A_y_v, 1),
                          tag="val.y")
    
    stack <- inla.stack(stk_aqum, 
                        stk_pcm,
                        stk_y_v, 
                        stk_y_e) 
    
    mod<- inla(formulas[[formula_id]],
               family=c("gaussian","gaussian","gaussian"),
               data=inla.stack.data(stack),
               control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
               # control.fixed = list(mean.intercept=3, prec.intercept=0.001),
               # control.mode = list(theta = c(
               #   log(1 / sd(pcm$no2_log, na.rm=T)^2), # theta[0] = [Log precision for the Gaussian observations]
               #   log(1 / sd(aqum$no2_log, na.rm=T)^2), # theta[1] = [Log precision for the Gaussian observations[2]]
               #   log(1 / sd(estim$no2_log, na.rm=T)^2), # theta[2] = [Log precision for the Gaussian observations[3]]
               #   log(1 / 0.1), # theta[3] = [Log precision for z2]
               #   log(range0), # theta[4] = [log(Range) for z1]
               #   log(1 / 1000), # theta[5] = [log(Stdev) for z1]
               #   log(1 / 0.1), # theta[6] = [Log precision for date.idx.no2]
               #   inla.models()$latent$ar1$hyper$theta2$to.theta(0.3), # theta[7] = [Rho_intern for date.idx.no2]
               #   1, # theta[8] = [Beta_intern for z2.copy]
               #   1), # theta[9] = [Beta_intern for z1.copy]
               #   restart = TRUE),
               control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
               control.inla = list(adaptive,int.strategy='eb'),
               verbose=TRUE)
    
    # saveRDS(mod, "mod_",formula_id,".rds")  
  }
  
  if(formula_id==4){
    
    ## ***** pcm *****
    A_pcm <- inla.spde.make.A(mesh=mesh, cbind(pcm$easting, pcm$northing))
    stk_pcm <- inla.stack(data=list(y=cbind(pcm[,paste0(pol,"_log")], NA, NA)),  
                          effects=list(list(intercept.pcm=rep(1,nrow(pcm))), list(psi.field=1:spde$n.spde)),
                          A=list(1,A_pcm),
                          tag="est.pcm")
    
    ## ***** Aqum *****
    stk_aqum <- inla.stack(data=list(y=cbind(NA, aqum[,paste0(pol,"_log")], NA)),  
                           effects=list(list(intercept.aqum=1, csi.field=aqum$date.idx)), 
                           A=list(1), 
                           tag="est.aqum")
    
    ## data stack: include all the effects
    A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
    stk_y_e <- inla.stack(data=list(y=cbind(NA,NA,estim[,paste0(pol,"_log")])),  
                          effects=list(list(csi.field.copy=estim[,paste0("date.idx.",pol)]),
                                       list(psi.field.copy=1:spde$n.spde),
                                       data.frame(intercept=1,
                                                  estim[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing")])), 
                          A=list(1, A_y_e, 1),
                          tag="est.y")
    
    ### validation scenario
    A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
    stk_y_v <- inla.stack(data=list(y=cbind(NA,NA,rep(NA,length(valid[,paste0(pol,"_log")])))),
                          effects=list(list(csi.field.copy=valid[,paste0("date.idx.",pol)]),
                                       list(psi.field.copy=1:spde$n.spde),
                                       data.frame(intercept=1,
                                                  valid[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing")])),
                          A=list(1, A_y_v, 1),
                          tag="val.y")
    
    stack <- inla.stack(stk_aqum, 
                        stk_pcm,
                        stk_y_v, 
                        stk_y_e) 
    
    mod<- inla(formulas[[formula_id]],
               family=c("gaussian","gaussian","gaussian"),
               data=inla.stack.data(stack),
               control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
               # control.fixed = list(mean.intercept=3, prec.intercept=0.001),
               control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
               control.inla = list(adaptive,int.strategy='eb'),
               verbose=TRUE)
    
    # saveRDS(mod, "mod_",formula_id,".rds")  
  }
  
  
  ############################################################################
  ########                  BILINEAR INTERPOLATION                   #########
  ############################################################################
  
  
  if(formula_id==5){
    
    ## data stack: include all the effects
    A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
    stk_y_e <- inla.stack(data=list(y=estim[,paste0(pol,"_log")]),  
                          effects=list(list(csi.field=estim[,paste0("date.idx.",pol)]),
                                       list(psi.field=1:spde$n.spde),
                                       data.frame(intercept=1,
                                                  estim[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing",
                                                           paste0("aqum_log_",pol),paste0("pcm_log_",pol))])), 
                          A=list(1,A_y_e,1),
                          tag="est.y")
    
    ### validation scenario
    A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
    stk_y_v <- inla.stack(data=list(y=rep(NA,length(valid[,paste0(pol,"_log")]))),
                          effects=list(list(csi.field=valid[,paste0("date.idx.",pol)]),
                                       list(psi.field=1:spde$n.spde),
                                       data.frame(intercept=1,
                                                  valid[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing",
                                                           paste0("aqum_log_",pol),paste0("pcm_log_",pol))])),
                          A=list(1,A_y_v,1),
                          tag="val.y")
    
    stack <- inla.stack(stk_y_v, 
                        stk_y_e) 
    
    mod<- inla(formulas[[formula_id]],
               family="gaussian",
               data=inla.stack.data(stack),
               control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
               # control.fixed = list(mean.intercept=3, prec.intercept=0.001),
               control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
               control.inla = list(adaptive,int.strategy='eb'),
               verbose=TRUE)
    
    # saveRDS(mod, "mod_",formula_id,".rds")
    
  }
  
  
  ############################################################################
  ########                   PLUG-IN (kriging)                       #########
  ############################################################################
  
  if(formula_id==6){
    
    A_aqum_est <- inla.spde.make.A(mesh=mesh,
                                   loc=as.matrix(aqum[,c("easting","northing")]))
    stk_aqum_no2_est <- inla.stack(data=list(y=aqum$no2_log),
                                   A=list(1,A_aqum_est),
                                   effects=list(list(intercept=rep(1,nrow(aqum)),
                                                     csi.field.no2.time = aqum$date.idx),
                                                list(csi.field.no2.space = c(1:spde$n.spde))),
                                   tag="est.aqum.no2")
    A_aqum_pred <- inla.spde.make.A(mesh=mesh,
                                    loc=as.matrix(final_dataset[,c("easting","northing")]))
    stk_aqum_no2_pred <- inla.stack(data=list(y=rep(NA, nrow(final_dataset))),
                                    A=list(1,A_aqum_pred),
                                    effects=list(list(intercept=rep(1,nrow(final_dataset)),
                                                      csi.field.no2.time = final_dataset$date.idx),
                                                 list(csi.field.no2.space = c(1:spde$n.spde))),
                                    tag="pred.aqum.no2")
    
    stk_aqum_no2 = inla.stack(stk_aqum_no2_est,stk_aqum_no2_pred)
    
    
    if(!file.exists("../mod_aqum.rds")){
      # aqum model 
      
      u=sd(aqum$no2_log)
      
      
      formula.aqum = y ~ -1 + intercept + f(csi.field.no2.time , model='rw1', scale.model = TRUE, constr = TRUE, 
                                            hyper = list(theta = list(prior="pc.prec", param=c(u,0.01)))) + f(csi.field.no2.space , model=spde)
      
      mod.aqum=inla(formula=formula.aqum,
                    family="gaussian",
                    data=inla.stack.data(stk_aqum_no2),
                    control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk_aqum_no2)),
                    # control.fixed = list(mean.intercept=2.5, prec.intercept=0.01),
                    control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
                    control.inla = list(adaptive,int.strategy='eb'),
                    #control.inla=list(strategy="laplace", int.strategy ="eb", h=1e-03),
                    verbose=TRUE)
      
      final_dataset$aqum_pred = mod.aqum$summary.fitted.values[inla.stack.index(stk_aqum_no2,"pred.aqum.no2")$data,"mean"] 
      
      saveRDS(mod.aqum,"mod_aqum.rds")
      rm(mod.aqum)
      
    }else{
      mod.aqum = readRDS("../mod_aqum.rds")
      final_dataset$aqum_pred = mod.aqum$summary.fitted.values[inla.stack.index(stk_aqum_no2,"pred.aqum.no2")$data,"mean"] 
      rm(mod.aqum)
    }
    
    A_pcm_est <- inla.spde.make.A(mesh=mesh, cbind(pcm$easting, pcm$northing))
    stk_pcm_no2_est <- inla.stack(data=list(y=pcm$no2_log),
                                  effects=list(list(intercept=rep(1,nrow(pcm)),
                                                    psi.field.no2.time = pcm$date.idx), 
                                               list(psi.field.no2=1:spde$n.spde)),
                                  A=list(1,A_pcm_est),
                                  tag="est.pcm.no2")
    
    
    A_pcm_pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(final_dataset[,c("easting","northing")]))
    stk_pcm_no2_pred <- inla.stack(data=list(y=rep(NA, nrow(final_dataset))),
                                   effects=list(list(intercept=rep(1,nrow(final_dataset)),
                                                     psi.field.no2.time = as.numeric(final_dataset$year)), 
                                                list(psi.field.no2=1:spde$n.spde)),
                                   A=list(1,A_pcm_pred),
                                   tag="pred.pcm.no2")
    
    stk_pcm_no2 = inla.stack(stk_pcm_no2_est,stk_pcm_no2_pred)
    
    if(!file.exists("../mod_pcm.rds")){
      
      formula.pcm = y ~ -1 + intercept + f(psi.field.no2 , model=spde) + f(psi.field.no2.time, model="iid")
      
      final_dataset$pcm_pred = mod.pcm$summary.fitted.values[inla.stack.index(stk_pcm_no2,"pred.pcm.no2")$data,"mean"] 
      
      saveRDS(mod.pcm,"mod_pcm.rds")
      rm(mod.pcm)
      
    }else{
      mod.pcm = readRDS("../mod_pcm.rds")
      final_dataset$pcm_pred = mod.pcm$summary.fitted.values[inla.stack.index(stk_pcm_no2,"pred.pcm.no2")$data,"mean"] 
      rm(mod.pcm)
    }
    
    valid = final_dataset[final_dataset$code %in% monitors_val_cheat[[data_id]] , ]
    estim = final_dataset[!(final_dataset$code %in% monitors_val_cheat[[data_id]]) , ]
    
    ## data stack: include all the effects
    A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
    stk_y_e <- inla.stack(data=list(y=estim[,paste0(pol,"_log")]),  
                          effects=list(list(csi.field=estim[,paste0("date.idx.",pol)]),
                                       list(psi.field=1:spde$n.spde),
                                       data.frame(intercept=1,
                                                  estim[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing",
                                                           "aqum_pred","pcm_pred")])), 
                          A=list(1,A_y_e,1),
                          tag="est.y")
    
    ### validation scenario
    A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
    stk_y_v <- inla.stack(data=list(y=rep(NA,length(valid[,paste0(pol,"_log")]))),
                          effects=list(list(csi.field=valid[,paste0("date.idx.",pol)]),
                                       list(psi.field=1:spde$n.spde),
                                       data.frame(intercept=1,
                                                  valid[,c(paste0("date.idx.",pol),"sitetype.idx","loc.idx","code",paste0("stURB.",pol),paste0("stRKS.",pol),"easting","northing",
                                                           "aqum_pred","pcm_pred")])),
                          A=list(1,A_y_v,1),
                          tag="val.y")
    
    stack <- inla.stack(stk_y_v, 
                        stk_y_e) 
    
    mod<- inla(formulas[[formula_id]],
               family="gaussian",
               data=inla.stack.data(stack),
               control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
               # control.fixed = list(mean.intercept=3, prec.intercept=0.001),
               control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
               control.inla = list(adaptive,int.strategy='eb'),
               verbose=TRUE)
    
    # saveRDS(mod, "mod_",formula_id,".rds")
    
  }
  
  
  ############################################################################
  ########          MODEL FROM Mukhopadhyay and Sahu (2017)          #########
  ############################################################################
  
  ##--- See http://www.soton.ac.uk/~sks/pollution_estimates/  for original code
  ##--- Our configuration is based on the authors' best model results
  
  # aqum only
  # non-stationarity
  if(formula_id>=7){
    
    # create necessary columns 
    
    final_dataset$month = as.numeric(substr(final_dataset$date,6,7))
    final_dataset$day = as.numeric(substr(final_dataset$date,9,10))
    
    final_dataset$aqum_log_URB = ifelse(final_dataset[,paste0("stURB.",pol)]==1, final_dataset[,paste0("aqum_log_",pol)], 0)
    final_dataset$aqum_log_RKS = ifelse(final_dataset[,paste0("stRKS.",pol)]==1, final_dataset[,paste0("aqum_log_",pol)], 0)
    final_dataset$pcm_log_URB = ifelse(final_dataset[,paste0("stURB.",pol)]==1, final_dataset[,paste0("pcm_log_",pol)], 0)
    final_dataset$pcm_log_RKS = ifelse(final_dataset[,paste0("stRKS.",pol)]==1, final_dataset[,paste0("pcm_log_",pol)], 0)
    
    # generate probability surface on a grid (same for both datasets)
    cand_coords=as.matrix(pop.grid[,1:2]); 
    cand_len = length(cand_coords[,1]);
    
    cand_coords1=rep(0,(2*cand_len));
    for(i in 1:cand_len)
    {
      cand_coords1[(2*(i-1))+1] = cand_coords[i,1];  # assign longitudes to odd indeces
      cand_coords1[(2*(i-1))+2] = cand_coords[i,2];  # assign latitudes to even indeces
    }
    grid_prob=grid[,3];
    
    #  MCMC via Gibbs not using default choices  
    
    # hyper-parameters for the prior distributions
    priors<-spT.priors(model="GPP",inv.var.prior=Gamm(2,1),beta.prior=Norm(0,10^10));
    
    # initial values for the model parameters
    initials<-spT.initials(model="GPP", sig2eps=0.01,sig2eta=0.5, beta=NULL, phi=0.005)#, sig2l=0.01);
    
    # input for spatial decay
    spatial.decay<-spT.decay(distribution="FIXED",value=0.001);
    
    # Define scale (response variable must be on natural scale)
    scale = "LOG"
    
    # select estim and valid subsets
    valid = final_dataset[final_dataset$code %in% monitors_val_cheat[[data_id]] , ]
    estim = final_dataset[!(final_dataset$code %in% monitors_val_cheat[[data_id]]) , ]
    
    # longlat coordinates for fitting and validation sets
    
    coordinates.estim<-unique(round(estim[,c("loc.idx","longitude","latitude")],5))
    coordinates.valid<-unique(round(valid[,c("loc.idx","longitude","latitude")],5))
    coords.all = rbind(coordinates.estim,coordinates.valid)
    
    set.seed(100)
    
    # Define knots (same for set1 and set2)
    knots.coords<-spT.grid.coords(Longitude=c(max(coords.all[,"longitude"]),min(coords.all[,"longitude"])),
                                  Latitude=c(max(coords.all[,"latitude"]),min(coords.all[,"latitude"])), 
                                  by=c(5,5)); # sort of mesh
    
    # run model  
    
    posT <- c(0, 0)
    
    XY <- Formula.matrix(formula=formulas[[formula_id]],
                         data=estim);
    Y <- XY[[1]]
    X <- as.matrix(XY[[2]])
    time.data<-list(1,length(Y)/length(coordinates.estim[,1]));
    
    data_time = unique(estim[,c("year","month","day")]);
    data_time$loc = 1:time.data[[2]];
    pred_time = unique(valid[,c("year","month","day")]);
    pos = which(data_time$year==pred_time$year[1] & data_time$month==pred_time$month[1] & data_time$day==pred_time$day[1]);
    posT[1] = pos - 1
    
    end <- dim(pred_time)[1]
    pos <- which(data_time$year==pred_time$year[end] & data_time$month==pred_time$month[end] & data_time$day==pred_time$day[end])
    posT[2] = pos -1
    
    if(formula_id==7 || formula_id==8){ # use new function in spAir for non-stationarity
      mod <- spAir::spT.Gibbs(formula=formulas[[formula_id]],
                              data=estim,
                              model="GPP",
                              coords=coordinates.estim[,c("longitude","latitude")],
                              knots.coords=knots.coords,
                              newcoords=coordinates.valid[,c("longitude","latitude")],
                              newdata=valid,
                              nItr=50000,
                              nBurn=45000,
                              tol.dist=0.0001,
                              distance.method="geodetic:km",
                              initials=initials,
                              priors=priors,
                              scale.transform=scale,
                              spatial.decay=spatial.decay,
                              cand_coords=cand_coords1,
                              cand_len=cand_len,
                              grid_prob=grid_prob,
                              report=1000,
                              predloc=posT);
    }
    if(formula_id==9){ # use old function in spTimer for stationary model
      mod <- spTimer::spT.Gibbs(formula=formulas[[formula_id]],
                                data=estim,
                                model="GPP",
                                coords=coordinates.estim[,c("longitude","latitude")],
                                knots.coords=knots.coords,
                                newcoords=coordinates.valid[,c("longitude","latitude")],
                                newdata=valid,
                                nItr=50000,
                                nBurn=45000,
                                tol.dist=0.0001,
                                distance.method="geodetic:km",
                                initials=initials,
                                priors=priors,
                                scale.transform=scale,
                                spatial.decay=spatial.decay,
                                report=1000,
                                predloc=posT);
    }
    
  }
  
  
  
}





save.image("workspace_results.RData")


q(save="no")
