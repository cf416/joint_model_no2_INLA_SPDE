remove(list = ls())

library(INLA)


load("workspace_data.RData")
pol = "no2"

n_days=1826
u=sd(aqum[,paste0(pol,"_log")])
q=sd(pcm[,paste0(pol,"_log")])


#=================================================== 
### Create mesh
#===================================================

knots <- seq(from=1, to=n_days, by=30) # about one knot per month
length(knots)
mesh1d <- inla.mesh.1d(loc=knots)

mesh.aqum = inla.mesh.2d(rbind(coordinates.aqum),
                         loc.domain=boundary@polygons[[1]]@Polygons[[1]]@coords,
                         max.edge = c(75000,40000),
                         offset = c(10000,30000),
                         cutoff=8000)
mesh.aqum$n
plot(mesh.aqum)

mesh.pcm = inla.mesh.2d(coordinates.pcm,
                        loc.domain=boundary@polygons[[1]]@Polygons[[1]]@coords,
                        max.edge = c(75000,25000),
                        offset = c(2000,20000),
                        cutoff=4000)

mesh.pcm$n 
plot(mesh.pcm)



#=================================================== 
### Construct the SPDE model 
#=================================================== 

range0.aqum <- min(c(diff(range(mesh.aqum$loc[, 1])), diff(range(mesh.aqum$loc[, 2])))) /5
spde.aqum <- inla.spde2.pcmatern(mesh=mesh.aqum, alpha=2, ### mesh and smoothness parameter
                                 prior.range=c(range0.aqum, 0.9), ### P(practic.range<range0)=0.9
                                 prior.sigma=c(10, 0.1)) ### P(sigma>10)=0.1

range0.pcm <- min(c(diff(range(mesh.pcm$loc[, 1])), diff(range(mesh.pcm$loc[, 2])))) /5
spde.pcm <- inla.spde2.pcmatern(mesh=mesh.pcm, alpha=2, ### mesh and smoothness parameter
                                prior.range=c(range0.pcm, 0.9), ### P(practic.range<range0)=0.9
                                prior.sigma=c(100, 0.5)) ### P(sigma>100)=0.5





#===================================================
# model PCM space only
#===================================================

A_pcm <- inla.spde.make.A(mesh=mesh.pcm, cbind(pcm$easting, pcm$northing))
stk_pcm <- inla.stack(data=list(y=pcm[,paste0(pol,"_log")]),
                      effects=list(list(intercept=rep(1,nrow(pcm))), list(z1=1:spde.pcm$n.spde)),
                      A=list(1,A_pcm),
                      tag="est.pcm")

formula = y ~ -1 + intercept + f(z1, model=spde.pcm)

mod=inla(formula=formula,
         family="gaussian",
         data=inla.stack.data(stk_pcm),
         control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk_pcm)),
         control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
         control.inla = list(strategy='adaptive',int.strategy='eb'),
         verbose=TRUE)

saveRDS(mod,"mod_PCM.rds")

print(summary(mod))




#===================================================
# model PCM space + time
#===================================================

A_pcm <- inla.spde.make.A(mesh=mesh.pcm, cbind(pcm$easting, pcm$northing))
stk_pcm <- inla.stack(data=list(y=pcm[,paste0(pol,"_log")]),
                      A=list(1,A_pcm),
                      effects=list(list(intercept = rep(1, nrow(pcm)),
                                        z2 = pcm$date.idx),
                                   list(z1 = c(1:spde.pcm$n.spde))),
                      tag="est.pcm")

n <- mesh.pcm$n
A <- matrix(1, ncol = n, nrow = 1)
e <- matrix(0, ncol = 1)

formula = y ~ -1 + intercept + 
  f(z1, model=spde.pcm, extraconstr = list(A=A, e=e)) + 
  f(z2 , model='rw1', scale.model=TRUE, constr=TRUE, hyper=list(theta=list(prior="pc.prec", param=c(q,0.1)))) 

mod=inla(formula=formula,
         family="gaussian",
         data=inla.stack.data(stk_pcm),
         control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk_pcm)),
         control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
         control.inla = list(strategy='adaptive',int.strategy='eb'),
         verbose=TRUE)

saveRDS(mod,"mod_PCM_st.rds")

print(summary(mod))



#===================================================
# model PCM space*time interaction
#===================================================

st_index <- inla.spde.make.index(name="z3",
                                 n.spde=spde.pcm$n.spde,
                                 n.group=max(pcm$date.idx))
A_pcm <- inla.spde.make.A(mesh=mesh.pcm,
                          loc=as.matrix(pcm[,c("easting","northing")]),
                          group=pcm$date.idx,
                          n.group=max(pcm$date.idx))

stk_pcm <- inla.stack(data=list(y=pcm[,paste0(pol,"_log")]),
                      A=list(1,A_pcm),
                      effects=list(data.frame(intercept = rep(1, nrow(pcm))),
                                   st_index),
                      tag="est.pcm")

n <- mesh.pcm$n
A <- matrix(1, ncol = n, nrow = 1)
e <- matrix(0, ncol = 1)

formula = y ~ -1 + intercept + f(z3, model=spde.pcm, extraconstr = list(A=A, e=e),
                                 group=z3.group, control.group=list(model="rw1",scale.model = TRUE, 
                                 hyper=list(theta = list(prior="pc.prec", param=c(q,0.1)))))

mod=inla(formula=formula,
         family="gaussian",
         data=inla.stack.data(stk_pcm),
         control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk_pcm)),
         control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
         control.inla = list(strategy='adaptive',int.strategy='eb'),
         verbose=TRUE)

saveRDS(mod,"mod_PCM_st_interaction.rds")

print(summary(mod))




#===================================================
# model AQUM time only
#===================================================
stk_aqum <- inla.stack(data=list(y=aqum[,paste0(pol,"_log")]),
                       effects=list(list(intercept=rep(1,length(aqum$date.idx)),
                                         z2=aqum$date.idx)),
                       A=list(1),
                       tag="est.aqum")

formula = y ~ -1 + intercept + f(z2 , model='rw1', scale.model = TRUE, constr = TRUE, 
                                 hyper = list(theta = list(prior="pc.prec", param=c(u,0.01)))) 

mod=inla(formula=formula,
         family="gaussian",
         data=inla.stack.data(stk_aqum),
         control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk_aqum)),
         control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
         control.inla = list(strategy='adaptive',int.strategy='eb'),
         verbose=TRUE)

saveRDS(mod,"mod_AQUM.rds")

print(summary(mod))



#===================================================
# model AQUM space+time
#===================================================

A_aqum <- inla.spde.make.A(mesh=mesh.aqum,
                           loc=as.matrix(aqum[,c("easting","northing")]))
stk_aqum <- inla.stack(data=list(y=aqum[,paste0(pol,"_log")]),
                       A=list(1,A_aqum),
                       effects=list(list(intercept=rep(1,nrow(aqum)),
                                         z2 = aqum$date.idx),
                                    list(z1 = c(1:spde.aqum$n.spde))),
                       tag="est.aqum")
n <- mesh.aqum$n
A <- matrix(1, ncol = n, nrow = 1)
e <- matrix(0, ncol = 1)

formula = y ~ -1 + intercept + f(z1, model=spde.aqum, extraconstr=list(A=A, e=e)) + 
              f(z2, model='rw1', scale.model=TRUE, constr=TRUE, hyper=list(theta=list(prior="pc.prec", param=c(u,0.01)))) 

mod=inla(formula=formula,
         family="gaussian",
         data=inla.stack.data(stk_aqum),
         control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk_aqum)),
         control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
         control.inla = list(strategy='adaptive',int.strategy='eb'),
         verbose=TRUE)

saveRDS(mod,"mod_AQUM_st.rds")

print(summary(mod))




#===================================================
# model AQUM space*time interaction
#===================================================

st_index <- inla.spde.make.index(name="z3",
                                 n.spde=spde.aqum$n.spde,
                                 n.group=mesh1d$n)
A_aqum <- inla.spde.make.A(mesh=mesh.aqum,
                           loc=as.matrix(aqum[,c("easting","northing")]),
                           group=aqum$date.idx,
                           group.mesh = mesh1d)

stk_aqum <- inla.stack(data=list(y=aqum[,paste0(pol,"_log")]),
                       A=list(1,A_aqum),
                       effects=list(intercept=rep(1,nrow(aqum)),st_index),
                       tag="est.aqum")

n <- mesh.aqum$n
A <- matrix(1, ncol = n, nrow = 1)
e <- matrix(0, ncol = 1)

formula = y ~ -1 + intercept + f(z3, model=spde.aqum, 
                                 group=z3.group, extraconstr=list(A=A, e=e), 
                                 control.group=list(model='rw1', scale.model=TRUE, 
                                                    hyper=list(theta=list(prior="pc.prec", param=c(u,0.01)))) )

mod = inla(formula=formula,
           family="gaussian",
           data=inla.stack.data(stk_aqum),
           control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk_aqum)),
           control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
           control.inla = list(strategy='adaptive',int.strategy='eb'),
           verbose=TRUE)

saveRDS(mod,"mod_AQUM_st_interaction.rds")

print(summary(mod))



