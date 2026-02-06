###main
#generate data
###################reasonable
#install.packages("glmtrans")
rm(list = ls())
setwd("C:/model_average/model_average/Final")
library(glmtrans)
library(glmnet)
library(MASS)
library(lars)
library(ncvreg)
library(matrans)
source("Q.agg.r")
source("method_function.r")
source("MMA.r")
source("Algorithm.r")
source("Opt_Tools.r")
source("PMA.r")
source("help2.r")
nt=150
ns=200
n.test=250
p=500
TIMES=200

#######################################################################
#### R_2=0.8
M=7
mse.est.mat=matrix(0,13,TIMES)
mse.pre.mat=matrix(0,13,TIMES)
weight.mat=matrix(0,M,TIMES)

for (t in 1:TIMES) {
  
  set.seed(5*t)
  ##conduct beta0
  betaT=c(rep(1,5),rep(0,p-5))
  
  ##design 1: homo-data
  ####design1.1:construction difference: exist;
  ##############numeric difference: small;
  K=8
  num_varible=1:5
  betaS=matrix(0,p,K)
  for (i in 1:5) {
    betaS[1:num_varible[i],i]=.8
    if(num_varible[i]<10){
      betaS[(num_varible[i]+1):10,i]=0
    }
    
  }
  betaS[,6]=c(rep(0,10),rep(.25,20),rep(0,p-30))
  betaS[,7]=c(rep(0,10),rep(.1,50),rep(0,p-60))
  betaS[,8]=c(rep(0,10),rep(0.5,10),rep(0,p-20))
  #betaS[,9]=betaT
  #betaS[,10]=betaT
  #for (i in 5:K) {
  #  betaS[1:10,i]=rev(betaS[1:10,(i-4)])
  #}
  
  ##conduct sample:X,Y
  SIGMA=matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      SIGMA[i,j]=.5^abs(i-j)
    }
  }
  Xt=as.matrix(mvrnorm(nt,rep(0,p),SIGMA))
  var_signal=var(as.vector(Xt%*%betaT))
  R2_targets=0.8
  sigma2_values=var_signal * (1 - R2_targets) / R2_targets
  Yt=crossprod(t(Xt),betaT)+rnorm(nt,0,sqrt(sigma2_values))
  Xs.list=list()
  Ys.list=list()
  for (k in 1:K) {
    Xs.list[[k]]=as.matrix(mvrnorm(ns,rep(0,p),SIGMA))
    Ys.list[[k]]=crossprod(t(Xs.list[[k]]),betaS[,k])+rnorm(ns,0,sqrt(sigma2_values))
  }
  ##text sample
  X.test=as.matrix(mvrnorm(n.test,rep(0,p),SIGMA))
  Y.test=crossprod(t(X.test),betaT)+rnorm(n.test,0,sqrt(sigma2_values))
  
  ##standardize
  Xt=scale(Xt)
  Yt=Yt-mean(Yt)
  for (k in 1:K) {
    Xs.list[[k]]=scale(Xs.list[[k]])
    Ys.list[[k]]=Ys.list[[k]]-mean(Ys.list[[k]])
  }
  #X.test=scale(X.test)
  #Y.test=Y.test-mean(Y.test)
  ####number of candidate model
  
  beta.est=matrix(0, nrow = p, ncol = 13)
  mse.est=rep(0,13)
  mse.pre=rep(0,13)
  ################################################################
  ################################################################
  ####method:DEMA
  res_DEMA=DEMA(Xt,Yt,Xs.list,Ys.list,M,c=NULL)
  beta.est[,1]=res_DEMA$beta.DEMA
  mse.est[1]=mse.fun(betaT,beta.est[,1],X.test)$est.err
  mse.pre[1]=mean((Y.test-X.test%*%beta.est[,1])^2)
  weight.mat[,t]=apply(res_DEMA$weight,1,mean)
  ################################################################
  ################################################################
  ####method:Trans-Fusion
  DATASET=list(X0=Xt,y0=Yt,X=Xs.list,y=Ys.list)
  res_TFU=TransFusion(dataset = DATASET,debias = T)
  beta.est[,2]=res_TFU$xplus
  mse.est[2]=mse.fun(betaT,beta.est[,2],X.test)$est.err
  mse.pre[2]=mean((Y.test-X.test%*%beta.est[,2])^2)
  ################################################################
  ################################################################
  ####method:Trans-Lasso
  res_TLA=Trans_lasso(DATASET)
  
  beta.est[,3]=res_TLA$x
  mse.est[3]=mse.fun(betaT,beta.est[,3],X.test)$est.err
  mse.pre[3]=mean((Y.test-X.test%*%beta.est[,3])^2)
  #############################################################
  ################################################################
  ####method:PMA(2020,zhang)
  res_PMA=PMA(Xt,Yt)
  beta.est[,4]=res_PMA$beta_est
  mse.est[4]=mse.fun(betaT,beta.est[,4],X.test)$est.err
  mse.pre[4]=mean((Y.test-X.test%*%beta.est[,4])^2)
  
  
  ################################################################
  ################################################################
  ####method:MCP-TAR
  res_MCP=cv.ncvreg(X=Xt,y=Yt,penalty="MCP")
  beta.est[,5]=as.vector(coef(res_MCP)[-1])
  mse.est[5]=mse.fun(betaT,beta.est[,5],X.test)$est.err
  mse.pre[5]=mean((Y.test-X.test%*%beta.est[,5])^2)
  
  ################################################################
  ################################################################
  ####method:SCAD-TAR
  res_SCAD=cv.ncvreg(Xt,Yt,penalty="SCAD")
  beta.est[,6]=as.vector(coef(res_SCAD)[-1])
  mse.est[6]=mse.fun(betaT,beta.est[,6],X.test)$est.err
  mse.pre[6]=mean((Y.test-X.test%*%beta.est[,6])^2)
  
  ################################################################
  ################################################################
  ####method:LASSO-TAR
  cv.las=cv.glmnet(Xt,Yt)
  lambda.min=cv.las$lambda.min
  res_LASSO=glmnet(Xt,Yt,lambda = lambda.min)
  beta.est[,7]=as.vector(coef(res_LASSO)[-1])
  mse.est[7]=mse.fun(betaT,beta.est[,7],X.test)$est.err
  mse.pre[7]=mean((Y.test-X.test%*%beta.est[,7])^2)
  
  ################################################################
  ################################################################
  ####method:MCP-COM
  X.com=Xs.list
  X.com[[K+1]]=Xt
  X.com.mat=do.call(rbind,X.com)
  Y.com=Ys.list
  Y.com[[K+1]]=Yt
  Y.com.mat=do.call(rbind,Y.com)
  
  
  res_MCP_com=cv.ncvreg(X=X.com.mat,y=Y.com.mat,penalty="MCP")
  beta.est[,8]=as.vector(coef(res_MCP_com)[-1])
  mse.est[8]=mse.fun(betaT,beta.est[,8],X.test)$est.err
  mse.pre[8]=mean((Y.test-X.test%*%beta.est[,8])^2)
  
  ################################################################
  ################################################################
  ####method:LASSO-COM
  cv.las=cv.glmnet(X.com.mat,Y.com.mat)
  lambda.min=cv.las$lambda.min
  res_LASSO_com=glmnet(X.com.mat,Y.com.mat,lambda = lambda.min)
  beta.est[,9]=as.vector(coef(res_LASSO_com)[-1])
  mse.est[9]=mse.fun(betaT,beta.est[,9],X.test)$est.err
  mse.pre[9]=mean((Y.test-X.test%*%beta.est[,9])^2)
  
  #############################################################
  ################################################################
  #######method:SCAD.com
  res_MCP_com=cv.ncvreg(X=X.com.mat,y=Y.com.mat,penalty="SCAD")
  beta.est[,10]=as.vector(coef(res_MCP_com)[-1])
  mse.est[10]=mse.fun(betaT,beta.est[,10],X.test)$est.err
  mse.pre[10]=mean((Y.test-X.test%*%beta.est[,10])^2)
  
  ################################################################
  ################################################################
  ####method:MMA
  res_MMA=MMA(Xt,Yt,M)
  beta.est[,11]=as.vector(res_MMA$beta.MMA)
  mse.est[11]=mse.fun(betaT,beta.est[,11],X.test)$est.err
  mse.pre[11]=mean((Y.test-X.test%*%beta.est[,11])^2)
  ################################################################
  ################################################################
  ####method:Trans-SMAP(Hu. 2023)
  data.y=c(list(Yt),Ys.list)
  data.x=c(list(Xt),Xs.list)
  data.smap=list(data.y=data.y,data.x=data.x)
  fit_SMAP=trans.smap(train.data = data.smap,nfold = 3,
                      bs.para =NULL,
                      if.penalty = T,
                      lm.set = 1:p,
                      pen.para = list(pen.nfold = 10, pen.lambda = NULL))
  res_SMAP=pred.transsmap(fit_SMAP,newdata = list(data.x=X.test),
                          bs.para = NULL,
                          if.lm = T,if.penalty = T)
  beta.est[,12]=res_SMAP$beta.ma
  mse.est[12]=mse.fun(betaT,beta.est[,12],X.test)$est.err
  mse.pre[12]=mean((Y.test-X.test%*%beta.est[,12])^2)
  
  ################################################################
  ################################################################
  ####method:Utrans
  Utrans_target=list(x=Xt,y=as.vector(Yt))
  Utrans_source=Map(function(x_mat, y_vec) {
    list(x = x_mat, y = as.vector(y_vec))
  }, Xs.list, Ys.list)
  D.utrans=list(target=Utrans_target,source=Utrans_source)
  res_UTRANS=UTrans(D.utrans,family = "gaussian",
                    valid.nfolds = 5,mode = "data",type = "lasso")
  beta.est[,13]=res_UTRANS$b1[-1]
  mse.est[13]=mse.fun(betaT,beta.est[,13],X.test)$est.err
  mse.pre[13]=mean((Y.test-X.test%*%beta.est[,13])^2)
  
  
  
  mse.est.mat[,t]=mse.est
  mse.pre.mat[,t]=mse.pre
  print(paste(t,"times finished"))
  print(mse.est)
}

MSE.EST.8=as.matrix(apply(mse.est.mat,1,mean),nrow=13)
MSE.EST.SD.8=as.matrix(apply(mse.est.mat,1,sd),nrow=13)
MSE.PRE.8=as.matrix(apply(mse.pre.mat,1,mean),nrow=13)
MSE.PRE.SD.8=as.matrix(apply(mse.pre.mat,1,sd),nrow=13)
weight.8=as.matrix(apply(weight.mat,1,mean),nrow=M)

#######################################################################
#######################################################################
#### (2)R_2=0.6
M=7
mse.est.mat=matrix(0,13,TIMES)
mse.pre.mat=matrix(0,13,TIMES)
weight.mat=matrix(0,M,TIMES)

for (t in 1:TIMES) {
  
  set.seed(5*t)
  ##conduct beta0
  betaT=c(rep(1,5),rep(0,p-5))
  
  ##design 1: homo-data
  ####design1.1:construction difference: exist;
  ##############numeric difference: small;
  K=8
  num_varible=1:5
  betaS=matrix(0,p,K)
  for (i in 1:5) {
    betaS[1:num_varible[i],i]=.8
    if(num_varible[i]<10){
      betaS[(num_varible[i]+1):10,i]=0
    }
    
  }
  betaS[,6]=c(rep(0,10),rep(.25,20),rep(0,p-30))
  betaS[,7]=c(rep(0,10),rep(.1,50),rep(0,p-60))
  betaS[,8]=c(rep(0,10),rep(0.5,10),rep(0,p-20))
  #betaS[,9]=betaT
  #betaS[,10]=betaT
  #for (i in 5:K) {
  #  betaS[1:10,i]=rev(betaS[1:10,(i-4)])
  #}
  
  ##conduct sample:X,Y
  SIGMA=matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      SIGMA[i,j]=.5^abs(i-j)
    }
  }
  Xt=as.matrix(mvrnorm(nt,rep(0,p),SIGMA))
  var_signal=var(as.vector(Xt%*%betaT))
  R2_targets=0.6
  sigma2_values=var_signal * (1 - R2_targets) / R2_targets
  Yt=crossprod(t(Xt),betaT)+rnorm(nt,0,sqrt(sigma2_values))
  Xs.list=list()
  Ys.list=list()
  for (k in 1:K) {
    Xs.list[[k]]=as.matrix(mvrnorm(ns,rep(0,p),SIGMA))
    Ys.list[[k]]=crossprod(t(Xs.list[[k]]),betaS[,k])+rnorm(ns,0,sqrt(sigma2_values))
  }
  ##text sample
  X.test=as.matrix(mvrnorm(n.test,rep(0,p),SIGMA))
  Y.test=crossprod(t(X.test),betaT)+rnorm(n.test,0,sqrt(sigma2_values))
  
  ##standardize
  Xt=scale(Xt)
  Yt=Yt-mean(Yt)
  for (k in 1:K) {
    Xs.list[[k]]=scale(Xs.list[[k]])
    Ys.list[[k]]=Ys.list[[k]]-mean(Ys.list[[k]])
  }
  #X.test=scale(X.test)
  #Y.test=Y.test-mean(Y.test)
  ####number of candidate model
  
  beta.est=matrix(0, nrow = p, ncol = 13)
  mse.est=rep(0,13)
  mse.pre=rep(0,13)
  ################################################################
  ################################################################
  ####method:DEMA
  res_DEMA=DEMA(Xt,Yt,Xs.list,Ys.list,M,c=NULL)
  beta.est[,1]=res_DEMA$beta.DEMA
  mse.est[1]=mse.fun(betaT,beta.est[,1],X.test)$est.err
  mse.pre[1]=mean((Y.test-X.test%*%beta.est[,1])^2)
  weight.mat[,t]=apply(res_DEMA$weight,1,mean)
  ################################################################
  ################################################################
  ####method:Trans-Fusion
  DATASET=list(X0=Xt,y0=Yt,X=Xs.list,y=Ys.list)
  res_TFU=TransFusion(dataset = DATASET,debias = T)
  beta.est[,2]=res_TFU$xplus
  mse.est[2]=mse.fun(betaT,beta.est[,2],X.test)$est.err
  mse.pre[2]=mean((Y.test-X.test%*%beta.est[,2])^2)
  ################################################################
  ################################################################
  ####method:Trans-Lasso
  res_TLA=Trans_lasso(DATASET)
  
  beta.est[,3]=res_TLA$x
  mse.est[3]=mse.fun(betaT,beta.est[,3],X.test)$est.err
  mse.pre[3]=mean((Y.test-X.test%*%beta.est[,3])^2)
  #############################################################
  ################################################################
  ####method:PMA(2020,zhang)
  res_PMA=PMA(Xt,Yt)
  beta.est[,4]=res_PMA$beta_est
  mse.est[4]=mse.fun(betaT,beta.est[,4],X.test)$est.err
  mse.pre[4]=mean((Y.test-X.test%*%beta.est[,4])^2)
  
  
  ################################################################
  ################################################################
  ####method:MCP-TAR
  res_MCP=cv.ncvreg(X=Xt,y=Yt,penalty="MCP")
  beta.est[,5]=as.vector(coef(res_MCP)[-1])
  mse.est[5]=mse.fun(betaT,beta.est[,5],X.test)$est.err
  mse.pre[5]=mean((Y.test-X.test%*%beta.est[,5])^2)
  
  ################################################################
  ################################################################
  ####method:SCAD-TAR
  res_SCAD=cv.ncvreg(Xt,Yt,penalty="SCAD")
  beta.est[,6]=as.vector(coef(res_SCAD)[-1])
  mse.est[6]=mse.fun(betaT,beta.est[,6],X.test)$est.err
  mse.pre[6]=mean((Y.test-X.test%*%beta.est[,6])^2)
  
  ################################################################
  ################################################################
  ####method:LASSO-TAR
  cv.las=cv.glmnet(Xt,Yt)
  lambda.min=cv.las$lambda.min
  res_LASSO=glmnet(Xt,Yt,lambda = lambda.min)
  beta.est[,7]=as.vector(coef(res_LASSO)[-1])
  mse.est[7]=mse.fun(betaT,beta.est[,7],X.test)$est.err
  mse.pre[7]=mean((Y.test-X.test%*%beta.est[,7])^2)
  
  ################################################################
  ################################################################
  ####method:MCP-COM
  X.com=Xs.list
  X.com[[K+1]]=Xt
  X.com.mat=do.call(rbind,X.com)
  Y.com=Ys.list
  Y.com[[K+1]]=Yt
  Y.com.mat=do.call(rbind,Y.com)
  
  
  res_MCP_com=cv.ncvreg(X=X.com.mat,y=Y.com.mat,penalty="MCP")
  beta.est[,8]=as.vector(coef(res_MCP_com)[-1])
  mse.est[8]=mse.fun(betaT,beta.est[,8],X.test)$est.err
  mse.pre[8]=mean((Y.test-X.test%*%beta.est[,8])^2)
  
  ################################################################
  ################################################################
  ####method:LASSO-COM
  cv.las=cv.glmnet(X.com.mat,Y.com.mat)
  lambda.min=cv.las$lambda.min
  res_LASSO_com=glmnet(X.com.mat,Y.com.mat,lambda = lambda.min)
  beta.est[,9]=as.vector(coef(res_LASSO_com)[-1])
  mse.est[9]=mse.fun(betaT,beta.est[,9],X.test)$est.err
  mse.pre[9]=mean((Y.test-X.test%*%beta.est[,9])^2)
  
  #############################################################
  ################################################################
  #######method:SCAD.com
  res_MCP_com=cv.ncvreg(X=X.com.mat,y=Y.com.mat,penalty="SCAD")
  beta.est[,10]=as.vector(coef(res_MCP_com)[-1])
  mse.est[10]=mse.fun(betaT,beta.est[,10],X.test)$est.err
  mse.pre[10]=mean((Y.test-X.test%*%beta.est[,10])^2)
  
  ################################################################
  ################################################################
  ####method:MMA
  res_MMA=MMA(Xt,Yt,M)
  beta.est[,11]=as.vector(res_MMA$beta.MMA)
  mse.est[11]=mse.fun(betaT,beta.est[,11],X.test)$est.err
  mse.pre[11]=mean((Y.test-X.test%*%beta.est[,11])^2)
  ################################################################
  ################################################################
  ####method:Trans-SMAP(Hu. 2023)
  data.y=c(list(Yt),Ys.list)
  data.x=c(list(Xt),Xs.list)
  data.smap=list(data.y=data.y,data.x=data.x)
  fit_SMAP=trans.smap(train.data = data.smap,nfold = 3,
                      bs.para =NULL,
                      if.penalty = T,
                      lm.set = 1:p,
                      pen.para = list(pen.nfold = 10, pen.lambda = NULL))
  res_SMAP=pred.transsmap(fit_SMAP,newdata = list(data.x=X.test),
                          bs.para = NULL,
                          if.lm = T,if.penalty = T)
  beta.est[,12]=res_SMAP$beta.ma
  mse.est[12]=mse.fun(betaT,beta.est[,12],X.test)$est.err
  mse.pre[12]=mean((Y.test-X.test%*%beta.est[,12])^2)
  
  ################################################################
  ################################################################
  ####method:Utrans
  Utrans_target=list(x=Xt,y=as.vector(Yt))
  Utrans_source=Map(function(x_mat, y_vec) {
    list(x = x_mat, y = as.vector(y_vec))
  }, Xs.list, Ys.list)
  D.utrans=list(target=Utrans_target,source=Utrans_source)
  res_UTRANS=UTrans(D.utrans,family = "gaussian",
                    valid.nfolds = 5,mode = "data",type = "lasso")
  beta.est[,13]=res_UTRANS$b1[-1]
  mse.est[13]=mse.fun(betaT,beta.est[,13],X.test)$est.err
  mse.pre[13]=mean((Y.test-X.test%*%beta.est[,13])^2)
  
  
  
  mse.est.mat[,t]=mse.est
  mse.pre.mat[,t]=mse.pre
  print(paste(t,"times finished"))
  print(mse.est)
}

MSE.EST.6=as.matrix(apply(mse.est.mat,1,mean),nrow=13)
MSE.EST.SD.6=as.matrix(apply(mse.est.mat,1,sd),nrow=13)
MSE.PRE.6=as.matrix(apply(mse.pre.mat,1,mean),nrow=13)
MSE.PRE.SD.6=as.matrix(apply(mse.pre.mat,1,sd),nrow=13)
weight.6=as.matrix(apply(weight.mat,1,mean),nrow=M)


##############################################################
#######################################################################
#### (3) R_2=0.4
M=7
mse.est.mat=matrix(0,13,TIMES)
mse.pre.mat=matrix(0,13,TIMES)
weight.mat=matrix(0,M,TIMES)

for (t in 1:TIMES) {
  
  set.seed(5*t)
  ##conduct beta0
  betaT=c(rep(1,5),rep(0,p-5))
  
  ##design 1: homo-data
  ####design1.1:construction difference: exist;
  ##############numeric difference: small;
  K=8
  num_varible=1:5
  betaS=matrix(0,p,K)
  for (i in 1:5) {
    betaS[1:num_varible[i],i]=.8
    if(num_varible[i]<10){
      betaS[(num_varible[i]+1):10,i]=0
    }
    
  }
  betaS[,6]=c(rep(0,10),rep(.25,20),rep(0,p-30))
  betaS[,7]=c(rep(0,10),rep(.1,50),rep(0,p-60))
  betaS[,8]=c(rep(0,10),rep(0.5,10),rep(0,p-20))
  #betaS[,9]=betaT
  #betaS[,10]=betaT
  #for (i in 5:K) {
  #  betaS[1:10,i]=rev(betaS[1:10,(i-4)])
  #}
  
  ##conduct sample:X,Y
  SIGMA=matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      SIGMA[i,j]=.5^abs(i-j)
    }
  }
  Xt=as.matrix(mvrnorm(nt,rep(0,p),SIGMA))
  var_signal=var(as.vector(Xt%*%betaT))
  R2_targets=0.4
  sigma2_values=var_signal * (1 - R2_targets) / R2_targets
  Yt=crossprod(t(Xt),betaT)+rnorm(nt,0,sqrt(sigma2_values))
  Xs.list=list()
  Ys.list=list()
  for (k in 1:K) {
    Xs.list[[k]]=as.matrix(mvrnorm(ns,rep(0,p),SIGMA))
    Ys.list[[k]]=crossprod(t(Xs.list[[k]]),betaS[,k])+rnorm(ns,0,sqrt(sigma2_values))
  }
  ##text sample
  X.test=as.matrix(mvrnorm(n.test,rep(0,p),SIGMA))
  Y.test=crossprod(t(X.test),betaT)+rnorm(n.test,0,sqrt(sigma2_values))
  
  ##standardize
  Xt=scale(Xt)
  Yt=Yt-mean(Yt)
  for (k in 1:K) {
    Xs.list[[k]]=scale(Xs.list[[k]])
    Ys.list[[k]]=Ys.list[[k]]-mean(Ys.list[[k]])
  }
  #X.test=scale(X.test)
  #Y.test=Y.test-mean(Y.test)
  ####number of candidate model
  
  beta.est=matrix(0, nrow = p, ncol = 13)
  mse.est=rep(0,13)
  mse.pre=rep(0,13)
  ################################################################
  ################################################################
  ####method:DEMA
  res_DEMA=DEMA(Xt,Yt,Xs.list,Ys.list,M,c=NULL)
  beta.est[,1]=res_DEMA$beta.DEMA
  mse.est[1]=mse.fun(betaT,beta.est[,1],X.test)$est.err
  mse.pre[1]=mean((Y.test-X.test%*%beta.est[,1])^2)
  weight.mat[,t]=apply(res_DEMA$weight,1,mean)
  ################################################################
  ################################################################
  ####method:Trans-Fusion
  DATASET=list(X0=Xt,y0=Yt,X=Xs.list,y=Ys.list)
  res_TFU=TransFusion(dataset = DATASET,debias = T)
  beta.est[,2]=res_TFU$xplus
  mse.est[2]=mse.fun(betaT,beta.est[,2],X.test)$est.err
  mse.pre[2]=mean((Y.test-X.test%*%beta.est[,2])^2)
  ################################################################
  ################################################################
  ####method:Trans-Lasso
  res_TLA=Trans_lasso(DATASET)
  
  beta.est[,3]=res_TLA$x
  mse.est[3]=mse.fun(betaT,beta.est[,3],X.test)$est.err
  mse.pre[3]=mean((Y.test-X.test%*%beta.est[,3])^2)
  #############################################################
  ################################################################
  ####method:PMA(2020,zhang)
  res_PMA=PMA(Xt,Yt)
  beta.est[,4]=res_PMA$beta_est
  mse.est[4]=mse.fun(betaT,beta.est[,4],X.test)$est.err
  mse.pre[4]=mean((Y.test-X.test%*%beta.est[,4])^2)
  
  
  ################################################################
  ################################################################
  ####method:MCP-TAR
  res_MCP=cv.ncvreg(X=Xt,y=Yt,penalty="MCP")
  beta.est[,5]=as.vector(coef(res_MCP)[-1])
  mse.est[5]=mse.fun(betaT,beta.est[,5],X.test)$est.err
  mse.pre[5]=mean((Y.test-X.test%*%beta.est[,5])^2)
  
  ################################################################
  ################################################################
  ####method:SCAD-TAR
  res_SCAD=cv.ncvreg(Xt,Yt,penalty="SCAD")
  beta.est[,6]=as.vector(coef(res_SCAD)[-1])
  mse.est[6]=mse.fun(betaT,beta.est[,6],X.test)$est.err
  mse.pre[6]=mean((Y.test-X.test%*%beta.est[,6])^2)
  
  ################################################################
  ################################################################
  ####method:LASSO-TAR
  cv.las=cv.glmnet(Xt,Yt)
  lambda.min=cv.las$lambda.min
  res_LASSO=glmnet(Xt,Yt,lambda = lambda.min)
  beta.est[,7]=as.vector(coef(res_LASSO)[-1])
  mse.est[7]=mse.fun(betaT,beta.est[,7],X.test)$est.err
  mse.pre[7]=mean((Y.test-X.test%*%beta.est[,7])^2)
  
  ################################################################
  ################################################################
  ####method:MCP-COM
  X.com=Xs.list
  X.com[[K+1]]=Xt
  X.com.mat=do.call(rbind,X.com)
  Y.com=Ys.list
  Y.com[[K+1]]=Yt
  Y.com.mat=do.call(rbind,Y.com)
  
  
  res_MCP_com=cv.ncvreg(X=X.com.mat,y=Y.com.mat,penalty="MCP")
  beta.est[,8]=as.vector(coef(res_MCP_com)[-1])
  mse.est[8]=mse.fun(betaT,beta.est[,8],X.test)$est.err
  mse.pre[8]=mean((Y.test-X.test%*%beta.est[,8])^2)
  
  ################################################################
  ################################################################
  ####method:LASSO-COM
  cv.las=cv.glmnet(X.com.mat,Y.com.mat)
  lambda.min=cv.las$lambda.min
  res_LASSO_com=glmnet(X.com.mat,Y.com.mat,lambda = lambda.min)
  beta.est[,9]=as.vector(coef(res_LASSO_com)[-1])
  mse.est[9]=mse.fun(betaT,beta.est[,9],X.test)$est.err
  mse.pre[9]=mean((Y.test-X.test%*%beta.est[,9])^2)
  
  #############################################################
  ################################################################
  #######method:SCAD.com
  res_MCP_com=cv.ncvreg(X=X.com.mat,y=Y.com.mat,penalty="SCAD")
  beta.est[,10]=as.vector(coef(res_MCP_com)[-1])
  mse.est[10]=mse.fun(betaT,beta.est[,10],X.test)$est.err
  mse.pre[10]=mean((Y.test-X.test%*%beta.est[,10])^2)
  
  ################################################################
  ################################################################
  ####method:MMA
  res_MMA=MMA(Xt,Yt,M)
  beta.est[,11]=as.vector(res_MMA$beta.MMA)
  mse.est[11]=mse.fun(betaT,beta.est[,11],X.test)$est.err
  mse.pre[11]=mean((Y.test-X.test%*%beta.est[,11])^2)
  ################################################################
  ################################################################
  ####method:Trans-SMAP(Hu. 2023)
  data.y=c(list(Yt),Ys.list)
  data.x=c(list(Xt),Xs.list)
  data.smap=list(data.y=data.y,data.x=data.x)
  fit_SMAP=trans.smap(train.data = data.smap,nfold = 3,
                      bs.para =NULL,
                      if.penalty = T,
                      lm.set = 1:p,
                      pen.para = list(pen.nfold = 10, pen.lambda = NULL))
  res_SMAP=pred.transsmap(fit_SMAP,newdata = list(data.x=X.test),
                          bs.para = NULL,
                          if.lm = T,if.penalty = T)
  beta.est[,12]=res_SMAP$beta.ma
  mse.est[12]=mse.fun(betaT,beta.est[,12],X.test)$est.err
  mse.pre[12]=mean((Y.test-X.test%*%beta.est[,12])^2)
  
  ################################################################
  ################################################################
  ####method:Utrans
  Utrans_target=list(x=Xt,y=as.vector(Yt))
  Utrans_source=Map(function(x_mat, y_vec) {
    list(x = x_mat, y = as.vector(y_vec))
  }, Xs.list, Ys.list)
  D.utrans=list(target=Utrans_target,source=Utrans_source)
  res_UTRANS=UTrans(D.utrans,family = "gaussian",
                    valid.nfolds = 5,mode = "data",type = "lasso")
  beta.est[,13]=res_UTRANS$b1[-1]
  mse.est[13]=mse.fun(betaT,beta.est[,13],X.test)$est.err
  mse.pre[13]=mean((Y.test-X.test%*%beta.est[,13])^2)
  
  
  
  mse.est.mat[,t]=mse.est
  mse.pre.mat[,t]=mse.pre
  print(paste(t,"times finished"))
  print(mse.est)
}

MSE.EST.4=as.matrix(apply(mse.est.mat,1,mean),nrow=13)
MSE.EST.SD.4=as.matrix(apply(mse.est.mat,1,sd),nrow=13)
MSE.PRE.4=as.matrix(apply(mse.pre.mat,1,mean),nrow=13)
MSE.PRE.SD.4=as.matrix(apply(mse.pre.mat,1,sd),nrow=13)
weight.4=as.matrix(apply(weight.mat,1,mean),nrow=M)
##############################################################
################################################################
#######################################################################
#### R_2=0.2
M=7
mse.est.mat=matrix(0,13,TIMES)
mse.pre.mat=matrix(0,13,TIMES)
weight.mat=matrix(0,M,TIMES)

for (t in 1:TIMES) {
  
  set.seed(5*t)
  ##conduct beta0
  betaT=c(rep(1,5),rep(0,p-5))
  
  ##design 1: homo-data
  ####design1.1:construction difference: exist;
  ##############numeric difference: small;
  K=8
  num_varible=1:5
  betaS=matrix(0,p,K)
  for (i in 1:5) {
    betaS[1:num_varible[i],i]=.8
    if(num_varible[i]<10){
      betaS[(num_varible[i]+1):10,i]=0
    }
    
  }
  betaS[,6]=c(rep(0,10),rep(.25,20),rep(0,p-30))
  betaS[,7]=c(rep(0,10),rep(.1,50),rep(0,p-60))
  betaS[,8]=c(rep(0,10),rep(0.5,10),rep(0,p-20))
  #betaS[,9]=betaT
  #betaS[,10]=betaT
  #for (i in 5:K) {
  #  betaS[1:10,i]=rev(betaS[1:10,(i-4)])
  #}
  
  ##conduct sample:X,Y
  SIGMA=matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      SIGMA[i,j]=.5^abs(i-j)
    }
  }
  Xt=as.matrix(mvrnorm(nt,rep(0,p),SIGMA))
  var_signal=var(as.vector(Xt%*%betaT))
  R2_targets=0.2
  sigma2_values=var_signal * (1 - R2_targets) / R2_targets
  Yt=crossprod(t(Xt),betaT)+rnorm(nt,0,sqrt(sigma2_values))
  Xs.list=list()
  Ys.list=list()
  for (k in 1:K) {
    Xs.list[[k]]=as.matrix(mvrnorm(ns,rep(0,p),SIGMA))
    Ys.list[[k]]=crossprod(t(Xs.list[[k]]),betaS[,k])+rnorm(ns,0,sqrt(sigma2_values))
  }
  ##text sample
  X.test=as.matrix(mvrnorm(n.test,rep(0,p),SIGMA))
  Y.test=crossprod(t(X.test),betaT)+rnorm(n.test,0,sqrt(sigma2_values))
  
  ##standardize
  Xt=scale(Xt)
  Yt=Yt-mean(Yt)
  for (k in 1:K) {
    Xs.list[[k]]=scale(Xs.list[[k]])
    Ys.list[[k]]=Ys.list[[k]]-mean(Ys.list[[k]])
  }
  #X.test=scale(X.test)
  #Y.test=Y.test-mean(Y.test)
  ####number of candidate model
  
  beta.est=matrix(0, nrow = p, ncol = 13)
  mse.est=rep(0,13)
  mse.pre=rep(0,13)
  ################################################################
  ################################################################
  ####method:DEMA
  res_DEMA=DEMA(Xt,Yt,Xs.list,Ys.list,M,c=NULL)
  beta.est[,1]=res_DEMA$beta.DEMA
  mse.est[1]=mse.fun(betaT,beta.est[,1],X.test)$est.err
  mse.pre[1]=mean((Y.test-X.test%*%beta.est[,1])^2)
  weight.mat[,t]=apply(res_DEMA$weight,1,mean)
  ################################################################
  ################################################################
  ####method:Trans-Fusion
  DATASET=list(X0=Xt,y0=Yt,X=Xs.list,y=Ys.list)
  res_TFU=TransFusion(dataset = DATASET,debias = T)
  beta.est[,2]=res_TFU$xplus
  mse.est[2]=mse.fun(betaT,beta.est[,2],X.test)$est.err
  mse.pre[2]=mean((Y.test-X.test%*%beta.est[,2])^2)
  ################################################################
  ################################################################
  ####method:Trans-Lasso
  res_TLA=Trans_lasso(DATASET)
  
  beta.est[,3]=res_TLA$x
  mse.est[3]=mse.fun(betaT,beta.est[,3],X.test)$est.err
  mse.pre[3]=mean((Y.test-X.test%*%beta.est[,3])^2)
  #############################################################
  ################################################################
  ####method:PMA(2020,zhang)
  res_PMA=PMA(Xt,Yt)
  beta.est[,4]=res_PMA$beta_est
  mse.est[4]=mse.fun(betaT,beta.est[,4],X.test)$est.err
  mse.pre[4]=mean((Y.test-X.test%*%beta.est[,4])^2)
  
  
  ################################################################
  ################################################################
  ####method:MCP-TAR
  res_MCP=cv.ncvreg(X=Xt,y=Yt,penalty="MCP")
  beta.est[,5]=as.vector(coef(res_MCP)[-1])
  mse.est[5]=mse.fun(betaT,beta.est[,5],X.test)$est.err
  mse.pre[5]=mean((Y.test-X.test%*%beta.est[,5])^2)
  
  ################################################################
  ################################################################
  ####method:SCAD-TAR
  res_SCAD=cv.ncvreg(Xt,Yt,penalty="SCAD")
  beta.est[,6]=as.vector(coef(res_SCAD)[-1])
  mse.est[6]=mse.fun(betaT,beta.est[,6],X.test)$est.err
  mse.pre[6]=mean((Y.test-X.test%*%beta.est[,6])^2)
  
  ################################################################
  ################################################################
  ####method:LASSO-TAR
  cv.las=cv.glmnet(Xt,Yt)
  lambda.min=cv.las$lambda.min
  res_LASSO=glmnet(Xt,Yt,lambda = lambda.min)
  beta.est[,7]=as.vector(coef(res_LASSO)[-1])
  mse.est[7]=mse.fun(betaT,beta.est[,7],X.test)$est.err
  mse.pre[7]=mean((Y.test-X.test%*%beta.est[,7])^2)
  
  ################################################################
  ################################################################
  ####method:MCP-COM
  X.com=Xs.list
  X.com[[K+1]]=Xt
  X.com.mat=do.call(rbind,X.com)
  Y.com=Ys.list
  Y.com[[K+1]]=Yt
  Y.com.mat=do.call(rbind,Y.com)
  
  
  res_MCP_com=cv.ncvreg(X=X.com.mat,y=Y.com.mat,penalty="MCP")
  beta.est[,8]=as.vector(coef(res_MCP_com)[-1])
  mse.est[8]=mse.fun(betaT,beta.est[,8],X.test)$est.err
  mse.pre[8]=mean((Y.test-X.test%*%beta.est[,8])^2)
  
  ################################################################
  ################################################################
  ####method:LASSO-COM
  cv.las=cv.glmnet(X.com.mat,Y.com.mat)
  lambda.min=cv.las$lambda.min
  res_LASSO_com=glmnet(X.com.mat,Y.com.mat,lambda = lambda.min)
  beta.est[,9]=as.vector(coef(res_LASSO_com)[-1])
  mse.est[9]=mse.fun(betaT,beta.est[,9],X.test)$est.err
  mse.pre[9]=mean((Y.test-X.test%*%beta.est[,9])^2)
  
  #############################################################
  ################################################################
  #######method:SCAD.com
  res_MCP_com=cv.ncvreg(X=X.com.mat,y=Y.com.mat,penalty="SCAD")
  beta.est[,10]=as.vector(coef(res_MCP_com)[-1])
  mse.est[10]=mse.fun(betaT,beta.est[,10],X.test)$est.err
  mse.pre[10]=mean((Y.test-X.test%*%beta.est[,10])^2)
  
  ################################################################
  ################################################################
  ####method:MMA
  res_MMA=MMA(Xt,Yt,M)
  beta.est[,11]=as.vector(res_MMA$beta.MMA)
  mse.est[11]=mse.fun(betaT,beta.est[,11],X.test)$est.err
  mse.pre[11]=mean((Y.test-X.test%*%beta.est[,11])^2)
  ################################################################
  ################################################################
  ####method:Trans-SMAP(Hu. 2023)
  data.y=c(list(Yt),Ys.list)
  data.x=c(list(Xt),Xs.list)
  data.smap=list(data.y=data.y,data.x=data.x)
  fit_SMAP=trans.smap(train.data = data.smap,nfold = 3,
                      bs.para =NULL,
                      if.penalty = T,
                      lm.set = 1:p,
                      pen.para = list(pen.nfold = 10, pen.lambda = NULL))
  res_SMAP=pred.transsmap(fit_SMAP,newdata = list(data.x=X.test),
                          bs.para = NULL,
                          if.lm = T,if.penalty = T)
  beta.est[,12]=res_SMAP$beta.ma
  mse.est[12]=mse.fun(betaT,beta.est[,12],X.test)$est.err
  mse.pre[12]=mean((Y.test-X.test%*%beta.est[,12])^2)
  
  ################################################################
  ################################################################
  ####method:Utrans
  Utrans_target=list(x=Xt,y=as.vector(Yt))
  Utrans_source=Map(function(x_mat, y_vec) {
    list(x = x_mat, y = as.vector(y_vec))
  }, Xs.list, Ys.list)
  D.utrans=list(target=Utrans_target,source=Utrans_source)
  res_UTRANS=UTrans(D.utrans,family = "gaussian",
                    valid.nfolds = 5,mode = "data",type = "lasso")
  beta.est[,13]=res_UTRANS$b1[-1]
  mse.est[13]=mse.fun(betaT,beta.est[,13],X.test)$est.err
  mse.pre[13]=mean((Y.test-X.test%*%beta.est[,13])^2)
  
  
  
  mse.est.mat[,t]=mse.est
  mse.pre.mat[,t]=mse.pre
  print(paste(t,"times finished"))
  print(mse.est)
}

MSE.EST.2=as.matrix(apply(mse.est.mat,1,mean),nrow=13)
MSE.EST.SD.2=as.matrix(apply(mse.est.mat,1,sd),nrow=13)
MSE.PRE.2=as.matrix(apply(mse.pre.mat,1,mean),nrow=13)
MSE.PRE.SD.2=as.matrix(apply(mse.pre.mat,1,sd),nrow=13)
weight.2=as.matrix(apply(weight.mat,1,mean),nrow=M)

#############################################################
###############################################################
#############################################################
############################################################
#############################################################
MSE.EST=cbind(MSE.EST.8, MSE.EST.6, MSE.EST.4, MSE.EST.2)
MSE.EST.SD=cbind(MSE.EST.SD.8,MSE.EST.SD.6,MSE.EST.SD.4,MSE.EST.SD.2)
MSE.PRE=cbind(MSE.PRE.8, MSE.PRE.6, MSE.PRE.4, MSE.PRE.2)
MSE.PRE.SD=cbind(MSE.PRE.SD.8,MSE.PRE.SD.6,MSE.PRE.SD.4,MSE.PRE.SD.2)
WEIGHT=cbind(weight.8, weight.6, weight.4, weight.2)


row.names(MSE.EST)=c("DEMA","TranFusion","Trans-Lasso","PMA",
                     "MCP","SCAD","LASSO",
                     "MCV.COM","LASSO.COM","SCAD.COM",
                     "MMA","Trans-SMAP","Utrans")
row.names(MSE.EST.SD)=c("DEMA","TranFusion","Trans-Lasso","PMA",
                        "MCP","SCAD","LASSO",
                        "MCV.COM","LASSO.COM","SCAD.COM",
                        "MMA","Trans-SMAP","Utrans")
row.names(MSE.PRE)=c("DEMA","TranFusion","Trans-Lasso","PMA",
                     "MCP","SCAD","LASSO",
                     "MCV.COM","LASSO.COM","SCAD.COM",
                     "MMA","Trans-SMAP","Utrans")

row.names(MSE.PRE.SD)=c("DEMA","TranFusion","Trans-Lasso","PMA",
                        "MCP","SCAD","LASSO",
                        "MCV.COM","LASSO.COM","SCAD.COM",
                        "MMA","Trans-SMAP","Utrans")




colnames(MSE.EST)=c("R2=0.8", "R2=0.6", "R2=0.4","R2=0.2")
colnames(MSE.EST.SD)=c("R2=0.8", "R2=0.6", "R2=0.4","R2=0.2")
colnames(MSE.PRE)=c("R2=0.8", "R2=0.6", "R2=0.4","R2=0.2")
colnames(MSE.PRE.SD)=c("R2=0.8", "R2=0.6", "R2=0.4","R2=0.2")
colnames(WEIGHT)=c("R2=0.8", "R2=0.6", "R2=0.4","R2=0.2")


write.csv(MSE.EST,"est1_2.csv")
write.csv(MSE.EST.SD,"estsd1_2.csv")
write.csv(MSE.PRE,"pre1_2.CSV")
write.csv(MSE.PRE.SD,"presd1_2.csv")
write.csv(WEIGHT,"weight1_2.csv")

