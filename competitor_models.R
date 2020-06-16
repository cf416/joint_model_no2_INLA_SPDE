remove(list = ls())


##--- NOTE:  
## This script is set up to be run as an array job. 
## data_id indicates the validation set and is received from the bash script;
## model_id indicates the model to run and is received from the bash script:
### 1. JOINT MODEL WITH AQUM ONLY 
### 2. JOINT MODEL WITH PCM ONLY 
### 3. BILINEAR INTERPOLATION
### 4. PLUG-IN (KRIGING)
### 5. NON-STATIONARY MODEL FROM Mukhopadhyay and Sahu (2017) WITH AQUM ONLY
### 6. NON-STATIONARY MODEL FROM Mukhopadhyay and Sahu (2017) WITH AQUM AND PCM
### 7. STATIONARY MODEL FROM Mukhopadhyay and Sahu (2017) WITH AQUM AND PCM
## If not received as parameters, values are set to 1
data_id <- ifelse(nchar(Sys.getenv("DATA_ID"))>0, as.numeric(Sys.getenv("DATA_ID")), 1) # 1 to 6
model_id <- ifelse(nchar(Sys.getenv("MODEL_ID"))>0, as.numeric(Sys.getenv("MODEL_ID")), 1) # 1 to 7

if (!file.exists(file.path(getwd(),paste0("Output_",model_id,"_no2_",data_id)))){
  dir.create(file.path(getwd(),paste0("Output_",model_id,"_no2_",data_id)))
} 
setwd(paste0("Output_",model_id,"_no2_",data_id))
print(getwd())

if(file.exists(file.path(getwd(),"/workspace_results.RData"))){
  print(paste0("Model ",model_id," has already run."))
}else{
  
  #=================================================== 
  ### Load necessary R packages
  #===================================================
  
  library(INLA)
  library(spAir) # download from http://www.soton.ac.uk/~sks/pollution_estimates/
  library(spTimer)
  
  ## Set seed for reproducibility
  set.seed(123)
  
  
  #=================================================== 
  ### Data preparation
  #===================================================
  
  load("../workspace_data.RData")
  
  pol = "no2"
  
  final_dataset$site.type = factor(as.character(final_dataset$site.type), ordered = T, levels = c("RUR","URB","RKS"))
  final_dataset$site.type.n = as.numeric(final_dataset$site.type)
  final_dataset$sitetype.idx = as.numeric(final_dataset$site.type)
  
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
  
  #=================================================== 
  ### Construct the SPDE model 
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
  
  # lambda23 = list(theta = list(prior = 'normal', param = c(0.9, 0.01), initial=0.9))  # lambda_2,3
  # lambda12 = list(theta = list(prior = 'normal', param = c(1.1, 0.01), initial=1.1))  # lambda_1,2
  # lambda13 = list(theta = list(prior = 'normal', param = c(1.3, 0.01), initial=1.3)) # lambda_1,3
  
  
  
  ##--- NOTE: the INLA call can be parallelized using Pardiso - see inla.pardiso()
  
  ############################################################################
  ########                JOINT MODEL WITH AQUM ONLY                 #########
  ############################################################################
  
  if(model_id==1){
    
    formula <- y ~ -1 + alpha1 + alpha2 + betaURB + betaRKS +  
      f(z1, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) +  
      f(z12, copy='z1', fixed = FALSE) + 
      f(z2, model='rw1', hyper = rw1.aqum.prior, scale.model = TRUE, constr = TRUE) + 
      f(z22, copy='z2', fixed = FALSE) + 
      f(z3, model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)
    
    A_aqum <- inla.spde.make.A(mesh=mesh,cbind(aqum$easting, aqum$northing))
    stk_aqum <- inla.stack(data=list(y=cbind(aqum[,paste0(pol,"_log")], NA)),  
                           A=list(1,A_aqum),
                           effects=list(list(alpha2=rep(1,nrow(aqum)),
                                             z2 = aqum$date.idx),
                                        list(z1 = c(1:spde$n.spde))),
                           tag="est.aqum")
    
    ## data stack: include all the effects
    A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
    stk_y_e <- inla.stack(data=list(y=cbind(NA,estim[,paste0(pol,"_log")])),  
                          effects=list(list(z22=estim[,paste0("date.idx.",pol)]),
                                       list(z12=c(1:spde$n.spde)),
                                       data.frame(alpha1=1,
                                                  betaURB=estim[,c(paste0("stURB.",pol))],
                                                  betaRKS=estim[,c(paste0("stRKS.",pol))],
                                                  z3=estim[,c(paste0("date.idx.",pol))],
                                                  estim[,c("sitetype.idx","loc.idx","code","easting","northing")])), 
                          A=list(1, A_y_e, 1),
                          tag="est.y")
    
    ### validation scenario
    A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
    stk_y_v <- inla.stack(data=list(y=cbind(NA,rep(NA,length(valid[,paste0(pol,"_log")])))),
                          effects=list(list(z22=valid[,paste0("date.idx.",pol)]),
                                       list(z12=c(1:spde$n.spde)),
                                       data.frame(alpha1=1,
                                                  betaURB=valid[,c(paste0("stURB.",pol))],
                                                  betaRKS=valid[,c(paste0("stRKS.",pol))],
                                                  z3=valid[,c(paste0("date.idx.",pol))],
                                                  valid[,c("sitetype.idx","loc.idx","code","easting","northing")])), 
                          A=list(1, A_y_v, 1),
                          tag="val.y")
    
    stack <- inla.stack(stk_aqum, 
                        stk_y_v, 
                        stk_y_e) 
    
    mod<- inla(formula,
               family=c("gaussian","gaussian"),
               data=inla.stack.data(stack),
               control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
               control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
               control.inla = list(strategy="adaptive",int.strategy='eb'),
               verbose=TRUE)
    
  }
  
  
  
  ############################################################################
  ########                JOINT MODEL WITH PCM ONLY                 #########
  ############################################################################
  
  if(model_id==2){
    
    formula <- y ~ -1 + alpha1 + alpha2 + betaURB + betaRKS + 
      f(z1, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) + 
      f(z12, copy='z1', fixed = FALSE) + 
      f(z3, model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)
    
    ## ***** pcm *****
    ### spatial only index
    ### projection matrix: use mesh and pcm data coordinates 
    A_pcm <- inla.spde.make.A(mesh=mesh, cbind(pcm$easting, pcm$northing))
    ### the pcm data stack
    stk_pcm <- inla.stack(data=list(y=cbind(pcm[,paste0(pol,"_log")], NA)),  
                          effects=list(list(alpha2=rep(1,nrow(pcm))), 
                                       list(z1=1:spde$n.spde)),
                          A=list(1,A_pcm),
                          tag="est.pcm")
    
    ## data stack: include all the effects
    A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
    stk_y_e <- inla.stack(data=list(y=cbind(NA,estim[,paste0(pol,"_log")])),  
                          effects=list(list(z12=1:spde$n.spde),
                                       data.frame(alpha1=1,
                                                  betaURB=estim[,c(paste0("stURB.",pol))],
                                                  betaRKS=estim[,c(paste0("stRKS.",pol))],
                                                  z3=estim[,c(paste0("date.idx.",pol))],
                                                  estim[,c("sitetype.idx","loc.idx","code","easting","northing")])), 
                          A=list(A_y_e, 1),
                          tag="est.y")
    
    ### validation scenario
    A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
    stk_y_v <- inla.stack(data=list(y=cbind(NA,rep(NA,length(valid[,paste0(pol,"_log")])))),
                          effects=list(list(psi.field.copy=1:spde$n.spde),
                                       data.frame(alpha1=1,
                                                  betaURB=valid[,c(paste0("stURB.",pol))],
                                                  betaRKS=valid[,c(paste0("stRKS.",pol))],
                                                  z3=valid[,c(paste0("date.idx.",pol))],
                                                  valid[,c("sitetype.idx","loc.idx","code","easting","northing")])), 
                          A=list(A_y_v, 1),
                          tag="val.y")
    
    stack <- inla.stack(stk_pcm,
                        stk_y_v, 
                        stk_y_e) 
    
    mod<- inla(formula,
               family=c("gaussian","gaussian"),
               data=inla.stack.data(stack),
               control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
               control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
               control.inla = list(strategy="adaptive",int.strategy='eb'),
               verbose=TRUE)
  }
  
  
  
  ############################################################################
  ########                  BILINEAR INTERPOLATION                   #########
  ############################################################################
  
  
  if(model_id==3){
    
    formula <- y ~ -1 + beta0 + betaURB + betaRKS +  aqum_log_no2 + pcm_log_no2 + 
      f(z1, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) +
      f(z2, model='rw1', hyper = rw1.aqum.prior, scale.model = TRUE, constr = TRUE) +
      f(z3, model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)
    
    ## data stack: include all the effects
    A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
    stk_y_e <- inla.stack(data=list(y=estim[,paste0(pol,"_log")]),  
                          effects=list(list(z1=1:spde$n.spde),
                                       list(z2=estim[,paste0("date.idx.",pol)]),
                                       data.frame(beta0=1,
                                                  betaURB=estim[,c(paste0("stURB.",pol))],
                                                  betaRKS=estim[,c(paste0("stRKS.",pol))],
                                                  z3=estim[,c(paste0("date.idx.",pol))],
                                                  estim[,c("sitetype.idx","loc.idx","code","easting","northing")])), 
                          A=list(A_y_e,1,1),
                          tag="est.y")
    
    ### validation scenario
    A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
    stk_y_v <- inla.stack(data=list(y=rep(NA,length(valid[,paste0(pol,"_log")]))),
                          effects=list(list(z1=1:spde$n.spde),
                                       list(z2=valid[,paste0("date.idx.",pol)]),
                                       data.frame(beta0=1,
                                                  betaURB=valid[,c(paste0("stURB.",pol))],
                                                  betaRKS=valid[,c(paste0("stRKS.",pol))],
                                                  z3=valid[,c(paste0("date.idx.",pol))],
                                                  valid[,c("sitetype.idx","loc.idx","code","easting","northing")])), 
                          A=list(A_y_v,1,1),
                          tag="val.y")
    
    stack <- inla.stack(stk_y_v, 
                        stk_y_e) 
    
    mod<- inla(formula,
               family="gaussian",
               data=inla.stack.data(stack),
               control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
               control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
               control.inla = list(strategy="adaptive",int.strategy='eb'),
               verbose=TRUE)
    
  }
  
  
  ############################################################################
  ########                   PLUG-IN (kriging)                       #########
  ############################################################################
  
  ## This model needs AQUM and PCM predictions at the monitor locations. 
  ## These are produced from the additive spatio-temporal models for AQUM and PCM respectively
  
  if(model_id==4){
    
    y ~ -1 + beta0 + betaURB + betaRKS +  aqum_pred + pcm_pred + 
      f(z1, model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) +
      f(z2, model='rw1', hyper = rw1.aqum.prior, scale.model = TRUE, constr = TRUE) +
      f(z3, model='ar1', hyper=ar1.time.prior, replicate = sitetype.idx, constr=TRUE, rankdef = 1)
    
    if(!file.exists("../mod_aqum_pred.rds")){
      # aqum model 
      
      u=sd(aqum$no2_log)
      
      formula.aqum = y ~ -1 + beta0 + 
        f(z1 , model=spde, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) + 
        f(z2 , model='rw1', scale.model = TRUE, constr = TRUE, 
          hyper = list(theta = list(prior="pc.prec", param=c(u,0.01))))
      
      A_aqum_est <- inla.spde.make.A(mesh=mesh,
                                     loc=as.matrix(aqum[,c("easting","northing")]))
      stk_aqum_no2_est <- inla.stack(data=list(y=aqum$no2_log),
                                     A=list(1,A_aqum_est),
                                     effects=list(list(beta0=rep(1,nrow(aqum)),
                                                       z2 = aqum$date.idx),
                                                  list(z1 = c(1:spde$n.spde))),
                                     tag="est.aqum.no2")
      A_aqum_pred <- inla.spde.make.A(mesh=mesh,
                                      loc=as.matrix(final_dataset[,c("easting","northing")]))
      stk_aqum_no2_pred <- inla.stack(data=list(y=rep(NA, nrow(final_dataset))),
                                      A=list(1,A_aqum_pred),
                                      effects=list(list(beta0=rep(1,nrow(final_dataset)),
                                                        z2 = final_dataset$date.idx),
                                                   list(z1 = c(1:spde$n.spde))),
                                      tag="pred.aqum.no2")
      
      stk_aqum_no2 = inla.stack(stk_aqum_no2_est,stk_aqum_no2_pred)
      
      
      
      mod.aqum=inla(formula=formula.aqum,
                    family="gaussian",
                    data=inla.stack.data(stk_aqum_no2),
                    control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk_aqum_no2)),
                    control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
                    control.inla = list(strategy="adaptive",int.strategy='eb'),
                    verbose=TRUE)
      
      final_dataset$aqum_pred = mod.aqum$summary.fitted.values[inla.stack.index(stk_aqum_no2,"pred.aqum.no2")$data,"mean"] 
      
      saveRDS(mod.aqum,"../mod_aqum_pred.rds")
      rm(mod.aqum)
      
    }else{
      mod.aqum = readRDS("../mod_aqum_pred.rds")
      final_dataset$aqum_pred = mod.aqum$summary.fitted.values[inla.stack.index(stk_aqum_no2,"pred.aqum.no2")$data,"mean"] 
      rm(mod.aqum)
    }
    
    if(!file.exists("../mod_pcm_pred.rds")){
      
      formula.pcm = y ~ -1 + beta0 + f(z1, model=spde.pcm, extraconstr = list(A=matrix(1,ncol=mesh$n,nrow=1), e=matrix(0,ncol=1))) + 
        f(z2 , model='rw1', scale.model=TRUE, constr=TRUE, hyper=list(theta=list(prior="pc.prec", param=c(q,0.1)))) 
     
      q=sd(pcm$no2_log)
      
      A_pcm_est <- inla.spde.make.A(mesh=mesh, cbind(pcm$easting, pcm$northing))
      stk_pcm_no2_est <- inla.stack(data=list(y=pcm$no2_log),
                                    effects=list(list(beta0=rep(1,nrow(pcm)),
                                                      z2 = pcm$date.idx), 
                                                 list(z1=1:spde$n.spde)),
                                    A=list(1,A_pcm_est),
                                    tag="est.pcm.no2")
      
      
      A_pcm_pred <- inla.spde.make.A(mesh=mesh, loc=as.matrix(final_dataset[,c("easting","northing")]))
      stk_pcm_no2_pred <- inla.stack(data=list(y=rep(NA, nrow(final_dataset))),
                                     effects=list(list(beta0=rep(1,nrow(final_dataset)),
                                                       z2 = as.numeric(final_dataset$year)), 
                                                  list(z1=1:spde$n.spde)),
                                     A=list(1,A_pcm_pred),
                                     tag="pred.pcm.no2")
      
      stk_pcm_no2 = inla.stack(stk_pcm_no2_est,stk_pcm_no2_pred)
      
      
      final_dataset$pcm_pred = mod.pcm$summary.fitted.values[inla.stack.index(stk_pcm_no2,"pred.pcm.no2")$data,"mean"] 
      
      saveRDS(mod.pcm,"../mod_pcm_pred.rds")
      rm(mod.pcm)
      
    }else{
      mod.pcm = readRDS("../mod_pcm_pred.rds")
      final_dataset$pcm_pred = mod.pcm$summary.fitted.values[inla.stack.index(stk_pcm_no2,"pred.pcm.no2")$data,"mean"] 
      rm(mod.pcm)
    }
    
    
    ## data stack: include all the effects
    A_y_e <- inla.spde.make.A(mesh=mesh, cbind(estim$easting,estim$northing))
    stk_y_e <- inla.stack(data=list(y=estim[,paste0(pol,"_log")]),  
                          effects=list(list(z2=estim[,paste0("date.idx.",pol)]),
                                       list(z1=1:spde$n.spde),
                                       data.frame(beta0=1,
                                                  betaURB=estim[,c(paste0("stURB.",pol))],
                                                  betaRKS=estim[,c(paste0("stRKS.",pol))],
                                                  z3=estim[,c(paste0("date.idx.",pol))],
                                                  estim[,c("sitetype.idx","loc.idx","code","easting","northing")])),  
                          A=list(1,A_y_e,1),
                          tag="est.y")
    
    ### validation scenario
    A_y_v <- inla.spde.make.A(mesh=mesh, cbind(valid$easting, valid$northing))
    stk_y_v <- inla.stack(data=list(y=rep(NA,length(valid[,paste0(pol,"_log")]))),
                          effects=list(list(z2=valid[,paste0("date.idx.",pol)]),
                                       list(z1=1:spde$n.spde),
                                       data.frame(beta0=1,
                                                  betaURB=valid[,c(paste0("stURB.",pol))],
                                                  betaRKS=valid[,c(paste0("stRKS.",pol))],
                                                  z3=valid[,c(paste0("date.idx.",pol))],
                                                  valid[,c("sitetype.idx","loc.idx","code","easting","northing")])), 
                          A=list(1,A_y_v,1),
                          tag="val.y")
    
    stack <- inla.stack(stk_y_v, 
                        stk_y_e) 
    
    mod<- inla(formula,
               family="gaussian",
               data=inla.stack.data(stack),
               control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stack)),
               control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
               control.inla = list(strategy="adaptive",int.strategy='eb'),
               verbose=TRUE)
    
  }
  
  
  ############################################################################
  ########          MODEL FROM Mukhopadhyay and Sahu (2017)          #########
  ############################################################################
  
  ##--- See http://www.soton.ac.uk/~sks/pollution_estimates/  for original code
  ##--- Our configuration is based on the authors' best model results
  
  if(model_id>=5){
    
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
    valid = final_dataset[final_dataset$code %in% monitors_val[[data_id]] , ]
    estim = final_dataset[!(final_dataset$code %in% monitors_val[[data_id]]) , ]
    
    # longlat coordinates for fitting and validation sets
    
    coordinates.estim<-unique(round(estim[,c("loc.idx","longitude","latitude")],5))
    coordinates.valid<-unique(round(valid[,c("loc.idx","longitude","latitude")],5))
    coords.all = rbind(coordinates.estim,coordinates.valid)
    
    set.seed(100)
    
    # Define knots (same for set1 and set2)
    knots.coords<-spT.grid.coords(Longitude=c(max(coords.all[,"longitude"]),min(coords.all[,"longitude"])),
                                  Latitude=c(max(coords.all[,"latitude"]),min(coords.all[,"latitude"])), 
                                  by=c(5,5)); # sort of mesh
    
    
    posT <- c(0, 0)
    
    XY <- Formula.matrix(formula=formulas[[model_id]],
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
    
    if(model_id==5){ # use new function in spAir for non-stationarity
      
      formula <- no2 ~ aqum_log_no2 + aqum_log_URB + stURB.no2 + aqum_log_RKS + stRKS.no2
      
      mod <- spAir::spT.Gibbs(formula=formula,
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
    
    if(model_id==6){ # use new function in spAir for non-stationarity
      
      formula <- pol ~ aqum_log_no2 + pcm_log_no2 + aqum_log_URB + pcm_log_URB + stURB.no2 + aqum_log_RKS + pcm_log_RKS + stRKS.no2
      
      mod <- spAir::spT.Gibbs(formula=formula,
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
    
    if(model_id==7){ # use old function in spTimer for stationary model
      
      formula <- pol ~ aqum_log_no2 + pcm_log_no2 + aqum_log_URB + pcm_log_URB + stURB.no2 + aqum_log_RKS + pcm_log_RKS + stRKS.no2
      
      mod <- spTimer::spT.Gibbs(formula=formula,
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
