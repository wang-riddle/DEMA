###data analysis
rm(list = ls())
##PRIMARY:LUSC
#install.packages("glmtrans")
library(glmtrans)
library(glmnet)
library(MASS)
library(lars)
library(ncvreg)
library(matrans)
source("Q.agg.r")
source("method_function.r")
source("PMA.r")
source("Opt_Tools.r")
source("Algorithm.r")
source("MMA.r")
load("LUAD_clin_GE_KEGG.Rdata")
load("LUSC_clin_GE_KEGG.Rdata")
load("SKCM_GE_Breslow.Rdata")
load("SKCM_GE_KEGG.Rdata")
#LUAD_clin_clean(142)------GE.KEGG.LUAD(142*4243)
#LUSC_clin_clean(89)------GE.KEGG.LUSC(89*4243)
row.names(LUAD_clin_clean)==row.names(GE.KEGG.LUAD)
row.names(LUSC_clin_clean)==row.names(GE.KEGG.LUSC)
Y1=LUAD_clin_clean[,1]
X1=GE.KEGG.LUAD
Y2=LUSC_clin_clean[,1]
X2=GE.KEGG.LUSC
###source and target
Xt=X2
Yt=Y2
Xs.list=list(X1)
Ys.list=list(Y1)
p=ncol(Xt)
K=1
safe_scale <- function(x) {
  col_means <- colMeans(x, na.rm = TRUE)
  col_sds <- apply(x, 2, sd, na.rm = TRUE)
  col_sds[col_sds == 0] <- 1  # 处理标准差为0的情况
  scaled <- sweep(x, 2, col_means, "-")
  scaled <- sweep(scaled, 2, col_sds, "/")
  scaled
}

  ##standardize
Xt=safe_scale(Xt)
Yt=scale(Yt)
for (k in 1:K) {
    Xs.list[[k]]=safe_scale(Xs.list[[k]])
    Ys.list[[k]]=scale(Ys.list[[k]])
}

Times=100
M=7
#mse.est.mat=matrix(0,13,Times)
mse.pre.mat=matrix(0,12,Times)
weight.mat=matrix(0,M,Times)


for (t in 47:Times) {
  set.seed(t)
  sam=sample(1:nrow(Xt),floor(nrow(Xt)/3))
  
  Xt.train=Xt[-sam,]
  Xt.test=Xt[sam,]
  Yt.train=Yt[-sam,]
  Yt.test=Yt[sam,]
  ####number of candidate model
  beta.est=matrix(0, nrow = p, ncol = 12)
  mse.est=rep(0,12)
  mse.pre=rep(0,12)
  ################################################################
  ################################################################
  ####method:DEMA
  res_DEMA=DEMA(Xt.train,Yt.train,Xs.list,Ys.list,M,alpha = 0.95,c=30)
  beta.est[,1]=res_DEMA$beta.DEMA
  #mse.est[1]=mse.fun(betaT,beta.est[,1],Xt.test)$est.err
  mse.pre[1]=mean((Yt.test-Xt.test%*%beta.est[,1])^2)
  weight.mat[,1]=apply(res_DEMA$weight,1,mean)
  
  ################################################################
  ################################################################
  ####method:Trans-Fusion
  DATASET=list(X0=Xt.train,y0=Yt.train,X=Xs.list,y=Ys.list)
  res_TFU=TransFusion(dataset = DATASET,debias = T)
  beta.est[,2]=res_TFU$xplus
  #mse.est[2]=mse.fun(betaT,beta.est[,2],X.test)$est.err
  mse.pre[2]=mean((Yt.test-Xt.test%*%beta.est[,2])^2)
  ################################################################
  ################################################################
  ####method:Trans-Lasso
  res_TLA=Trans_lasso(DATASET)
  beta.est[,3]=res_TLA$x
  #mse.est[3]=mse.fun(betaT,beta.est[,3],X.test)$est.err
  mse.pre[3]=mean((Yt.test-Xt.test%*%beta.est[,3])^2)
  #############################################################
  ################################################################
  ####method:PMA(2020,zhang)
  res_PMA=PMA(Xt.train,Yt.train)
  beta.est[,4]=res_PMA$beta_est
  #mse.est[4]=mse.fun(betaT,beta.est[,4],X.test)$est.err
  mse.pre[4]=mean((Yt.test-Xt.test%*%beta.est[,4])^2)
  
  
  ################################################################
  ################################################################
  ####method:MCP-TAR
  res_MCP=cv.ncvreg(X=Xt.train,y=Yt.train,penalty="MCP")
  beta.est[,5]=as.vector(coef(res_MCP)[-1])
  #mse.est[5]=mse.fun(betaT,beta.est[,5],X.test)$est.err
  mse.pre[5]=mean((Yt.test-Xt.test%*%beta.est[,5])^2)
  
  ################################################################
  ################################################################
  ####method:SCAD-TAR
  res_SCAD=cv.ncvreg(Xt.train,Yt.train,penalty="SCAD")
  beta.est[,6]=as.vector(coef(res_SCAD)[-1])
  #mse.est[6]=mse.fun(betaT,beta.est[,6],X.test)$est.err
  mse.pre[6]=mean((Yt.test-Xt.test%*%beta.est[,6])^2)
  
  ################################################################
  ################################################################
  ####method:LASSO-TAR
  cv.las=cv.glmnet(Xt.train,Yt.train)
  lambda.min=cv.las$lambda.min
  res_LASSO=glmnet(Xt.train,Yt.train,lambda = lambda.min)
  beta.est[,7]=as.vector(coef(res_LASSO)[-1])
  #mse.est[7]=mse.fun(betaT,beta.est[,7],X.test)$est.err
  mse.pre[7]=mean((Yt.test-Xt.test%*%beta.est[,7])^2)
  
  ################################################################
  ################################################################
  ####method:MCP-COM
  X.com=Xs.list
  X.com[[K+1]]=Xt.train
  X.com.mat=do.call(rbind,X.com)
  Y.com=Ys.list
  Y.com[[K+1]]=as.matrix(Yt.train)
  Y.com.mat=do.call(rbind,Y.com)
  
  
  res_MCP_com=cv.ncvreg(X=X.com.mat,y=Y.com.mat,penalty="MCP")
  beta.est[,8]=as.vector(coef(res_MCP_com)[-1])
  #mse.est[8]=mse.fun(betaT,beta.est[,8],X.test)$est.err
  mse.pre[8]=mean((Yt.test-Xt.test%*%beta.est[,8])^2)
  
  ################################################################
  ################################################################
  ####method:LASSO-COM
  cv.las=cv.glmnet(X.com.mat,Y.com.mat)
  lambda.min=cv.las$lambda.min
  res_LASSO_com=glmnet(X.com.mat,Y.com.mat,lambda = lambda.min)
  beta.est[,9]=as.vector(coef(res_LASSO_com)[-1])
  #mse.est[9]=mse.fun(betaT,beta.est[,9],X.test)$est.err
  mse.pre[9]=mean((Yt.test-Xt.test%*%beta.est[,9])^2)
  
  #############################################################
  ################################################################
  #######method:SCAD.com
  res_MCP_com=cv.ncvreg(X=X.com.mat,y=Y.com.mat,penalty="SCAD")
  beta.est[,10]=as.vector(coef(res_MCP_com)[-1])
  #mse.est[10]=mse.fun(betaT,beta.est[,10],X.test)$est.err
  mse.pre[10]=mean((Yt.test-Xt.test%*%beta.est[,10])^2)
  
  ################################################################
  ################################################################
  ####method:MMA
  res_MMA=MMA(Xt.train,Yt.train,M)
  beta.est[,11]=as.vector(res_MMA$beta.MMA)
  #mse.est[11]=mse.fun(betaT,beta.est[,11],X.test)$est.err
  mse.pre[11]=mean((Yt.test-Xt.test%*%beta.est[,11])^2)
  ################################################################
  ################################################################
  ####method:Trans-SMAP(Hu. 2023)
  data.y=c(list(as.matrix(Yt.train)),Ys.list)
  data.x=c(list(Xt.train),Xs.list)
  data.smap=list(data.y=data.y,data.x=data.x)
  fit_SMAP=trans.smap(train.data = data.smap,nfold = 2,
                      bs.para =NULL,
                      if.penalty = T,
                      lm.set = 1:p,
                      pen.para = list(pen.nfold = 10, pen.lambda = NULL))
  res_SMAP=pred.transsmap(fit_SMAP,newdata = list(data.x=Xt.test),
                          bs.para = NULL,
                          if.lm = T,if.penalty = T)
  beta.est[,12]=res_SMAP$beta.ma
  #mse.est[12]=mse.fun(betaT,beta.est[,12],X.test)$est.err
  mse.pre[12]=mean((Yt.test-Xt.test%*%beta.est[,12])^2)
  
  
  
  #mse.est.mat[,t]=mse.est
  mse.pre.mat[,t]=mse.pre
  print(paste(t,"times finished"))
  print(mse.pre)
}
mse.pre.mat_clean=mse.pre.mat
mse.pre.mat_clean[is.na(mse.pre.mat)| mse.pre.mat == 0]=NA
#MSE.EST=as.matrix(apply(mse.est.mat,1,mean),nrow=9)
MSE.PRE.MEAN=as.matrix(apply(mse.pre.mat_clean,1,mean, na.rm = TRUE),nrow=12)
MSE.PRE.SD=as.matrix(apply(mse.pre.mat_clean,1,sd, na.rm = TRUE),nrow=12)
#mse.pre.mat1=mse.pre.mat[-8,]
row.names(MSE.PRE.MEAN)=c("DEMA","TranFusion","Trans-Lasso","PMA",
                          "MCP","SCAD","LASSO",
                          "MCV.COM","LASSO.COM","SCAD.COM",
                          "MMA","Trans-SMAP")

row.names(MSE.PRE.SD)=c("DEMA","TranFusion","Trans-Lasso","PMA",
                        "MCP","SCAD","LASSO",
                        "MCV.COM","LASSO.COM","SCAD.COM",
                        "MMA","Trans-SMAP")



weight=as.matrix(apply(weight.mat,1,mean),nrow=M)
write.csv(MSE.PRE.MEAN,"realdata_LUSC_mean.csv")
write.csv(MSE.PRE.SD,"realdata_LUSC_sd.csv")
  
  
  








