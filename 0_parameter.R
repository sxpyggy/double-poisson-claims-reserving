library(data.table)
library("ChainLadder")
library("CASdatasets")
library(mvtnorm)
### y can be negative
rm(list=ls())
deviance_ddp<-function(y,mu){
   ifelse(y<0,(y-mu)^2/mu,2*(mu-y+y*ifelse(y==0,0,log(abs(y)/mu))))
}

pearson_poi<-function(y,mu){
  (y-mu)^2/mu
}

textwidth<-7.3 #inch

res_lim<-c(-2.5,2.5)

DP<-function(data_mat,disp,fix_mean,K,verbose){
  ## data_mat is a data frame with accident year (acc), development year (dev),
  ## incremental claims (inc)
  ## disp="~ 1" would lead an ODP model
  data_mat$acc<-as.factor(data_mat$acc)
  data_mat$dev<-as.factor(data_mat$dev)
  data_mat$j<-as.numeric(data_mat$j)
  data_mat$inc<-as.numeric(data_mat$inc)
  odp<-glm(inc ~ acc + dev, family = quasipoisson,
           data=data_mat[which(data_mat$inc>0),])
  alpha_odp<-odp$coef
  # if(length(alpha_odp)<2*max(data_mat$j)-1) {
  #   alpha_odp<-c(alpha_odp,alpha_odp[length(alpha_odp)])
  #   }
  YY<-data_mat$inc
  TT<-model.matrix( ~ acc + dev , data=data_mat)
  SS<-model.matrix( as.formula(disp) , data=data_mat)
  beta<-c(log(summary(odp)$dispersion),rep(0,ncol(SS)-1))
  alpha<-alpha_odp
  eta<-TT%*%alpha
  lambda<-SS%*%beta
  theta<-exp(lambda)
  mu<-exp(eta)
  for (k in 1:K){
    DD<-abs(deviance_ddp(YY,mu))
    W1<-diag(c(mu/theta))
    HH<-diag(W1^0.5%*%TT%*%solve(t(TT)%*%W1%*%(TT))%*%t(TT)%*%W1^0.5)
    HH<-ifelse(round(HH,1)==1,0.99,HH)
    W2<-diag(c((1-HH)/2))
    Z2<-DD/(1-HH)/theta-1+lambda
    beta<-solve(t(SS)%*%W2%*%(SS))%*%t(SS)%*%W2%*%Z2
    lambda<-SS%*%beta
    theta<-exp(lambda)
    if (fix_mean==F) {
      W1<-diag(c(mu/theta))
      Z1<-(YY-mu)/mu+eta
      alpha<-solve(t(TT)%*%W1%*%(TT))%*%t(TT)%*%W1%*%Z1
      eta<-TT%*%alpha
      mu<-exp(eta)
    } 
    if (verbose==T){
    print(mean(theta))
    print(sum((alpha_odp-alpha)^2)) 
    }
  }
  fisher_alpha_dp<-t(TT)%*%W1%*%(TT)
  fisher_beta_dp<-t(SS)%*%W2%*%(SS)
  SS1<-SS[data_mat$acc==1,]
  list(mu=mu,DD=DD,theta=theta,disp_vary=(exp(SS1%*%beta)),
       alpha=alpha,beta=beta,D=sum(-DD/theta-(1-HH)*log(theta)),
       fisher_alpha=fisher_alpha_dp,
       fisher_beta=fisher_beta_dp
       )
}

EV<-function(data_mat){
  ## data_mat is a data frame with accident year (acc), 
  ## development year (dev), incremental claims (inc)
  data_mat$acc<-as.factor(data_mat$acc)
  data_mat$dev<-as.factor(data_mat$dev)
  data_mat$j<-as.numeric(data_mat$j)
  data_mat$inc<-as.numeric(data_mat$inc)
  data_mat$jj<-pmin(data_mat$j,max(data_mat$j)-1)
  odp<-glm(inc ~ acc + dev, family = quasipoisson, data=data_mat)
  alpha<-odp$coef
  N<-nrow(data_mat)
  p<-length(alpha)
  YY<-data_mat$inc
  TT<-model.matrix( ~ acc + dev , data=data_mat)
  eta<-TT%*%alpha
  mu<-exp(eta)
  DD<-abs(deviance_ddp(YY,mu))
  data_mat$DD<-DD
  disp_ev<-aggregate(data_mat$DD,by=list(data_mat$jj),mean)
  disp_ev$x<-disp_ev$x*(N/(N-p))
  disp_ev<-rbind(disp_ev,disp_ev[nrow(disp_ev),])
  disp_ev$Group.1[nrow(disp_ev)]<-max(data_mat$j)
  names(disp_ev)<-c("dev","theta")
  data_mat$theta<-NA
  for (i in 1:nrow(data_mat)){
    data_mat$theta[i]<-disp_ev$theta[disp_ev$dev==data_mat$dev[i]]
  }
  theta<-data_mat$theta
  W1<-diag(c(mu/theta))
  HH<-diag(W1^0.5%*%TT%*%solve(t(TT)%*%W1%*%(TT))%*%t(TT)%*%W1^0.5)
  fisher_alpha<-t(TT)%*%W1%*%(TT)
  list(mu=mu,DD=DD,theta=theta,disp_vary=disp_ev$theta,
       alpha=alpha,D=sum(-DD/theta-(1-HH)*log(theta)),
       fisher_alpha=fisher_alpha
       )
}

acc_pred<-function(I,disp,method,alpha,beta,fisher_alpha){
  pred_acc<-data.frame(acc=rep(NA,I*(I-1)/2+1),
                       dev=rep(NA,I*(I-1)/2+1),
                       j=rep(NA,I*(I-1)/2+1))
  pred_acc[1,]<-c(1,1,1)
  k=2
  for (i in 2:I){
    for (j in (I-i+2):I){
      pred_acc$acc[k]<-i
      pred_acc$dev[k]<-j
      k<-k+1
    }
  }
  pred_acc$j<-pred_acc$dev
  pred_acc$acc<-factor(pred_acc$acc)
  pred_acc$dev<-factor(pred_acc$dev)
  pred_TT<-model.matrix(~ acc + dev, data=pred_acc)
  if (method!="ev") {
    pred_SS<-model.matrix(as.formula(disp),data=pred_acc)
    pred_acc$theta<-exp(pred_SS%*%beta)
    }
  if (method=="ev") {
    for (i in 1:nrow(pred_acc)){
      pred_acc$theta[i]<-ev$theta[data_mat$dev==pred_acc$dev[i]][1]
    }
  }
  pred_acc$inc<-exp(pred_TT%*%alpha)
  pred_acc$pro_v<-pred_acc$theta*pred_acc$inc
  r_est<-aggregate(pred_acc$inc,by=list(pred_acc$acc),sum)$V1[-1]
  temp_mat<-pred_TT*c(pred_acc$inc);dim(temp_mat)
  xmu_mat<-matrix(NA,nrow=I,ncol=ncol(temp_mat))
  xmu_mat[1,]<-temp_mat[2,]
  for (i in 3:I){
    ind<-which(pred_acc$acc==i)
    xmu_mat[i-1,]<-apply(temp_mat[ind,],2,sum)
  }
  xmu_mat[I,]<-apply(temp_mat[-1,],2,sum)
  est_v<-rep(NA,I)
  for (i in 1:I){
    est_v[i]<-t(xmu_mat[i,])%*%solve(fisher_alpha)%*%c(xmu_mat[i,])
  }
  pro_v<-aggregate(pred_acc$pro_v,by=list(pred_acc$acc),sum)$V1[-1]
  pro_v<-c(pro_v,sum(pro_v))
  data.frame(mean=c(r_est,sum(r_est)),est_err=sqrt(est_v),
             pro_err=sqrt(pro_v),pre_err=sqrt(est_v+pro_v),
             cv=sqrt(est_v+pro_v)/c(r_est,sum(r_est))
             )
}

acc_pred2<-function(I,disp,method,alpha,beta,fisher_alpha,fisher_beta,seed){
  ## this function does not work well for non-gaussian distributed estimated
  ## parameters
  set.seed(seed)
  alpha<- t(rmvnorm(1, mean=alpha,sigma=solve(fisher_alpha)))
  beta<- t(rmvnorm(1, mean=beta,sigma=solve(fisher_beta)))
  pred_acc<-data.frame(acc=rep(NA,I*(I-1)/2+1),
                       dev=rep(NA,I*(I-1)/2+1),
                       j=rep(NA,I*(I-1)/2+1))
  pred_acc[1,]<-c(1,1,1)
  k=2
  for (i in 2:I){
    for (j in (I-i+2):I){
      pred_acc$acc[k]<-i
      pred_acc$dev[k]<-j
      k<-k+1
    }
  }
  pred_acc$j<-pred_acc$dev
  pred_acc$acc<-factor(pred_acc$acc)
  pred_acc$dev<-factor(pred_acc$dev)
  pred_TT<-model.matrix(~ acc + dev, data=pred_acc)
  if (method!="ev") {
    pred_SS<-model.matrix(as.formula(disp),data=pred_acc)
    pred_acc$theta<-exp(pred_SS%*%beta)
  }
  if (method=="ev") {
    for (i in 1:nrow(pred_acc)){
      pred_acc$theta[i]<-ev$theta[data_mat$dev==pred_acc$dev[i]][1]
    }
  }
  pred_acc$mu<-exp(pred_TT%*%alpha)
  set.seed(seed)
  pred_acc$inc<-pred_acc$theta*
    rpois(nrow(pred_acc),pred_acc$mu/pred_acc$theta)
  pred_acc[-1,c("acc","dev","inc","mu")]
}

BS<-function(data_mat,mu,theta,seed){
  ## this function bootstrap the residuals and return psedo data
  n<-nrow(data_mat)
  data_mat$mu<-mu
  data_mat$theta<-theta
  data_mat$scaled_pearson<-(data_mat$inc-data_mat$mu)/
    sqrt(data_mat$mu*data_mat$theta)
  set.seed(seed)
  data_mat$bs_pearson<-sample(data_mat$scaled_pearson,size = n,replace = T)
  data_mat$bs_inc<-(data_mat$mu+
    data_mat$bs_pearson*sqrt(data_mat$mu*data_mat$theta))
  
  ## this would underestimate reserves; 
  ## but it can avoid most errors in parameter estimation
  data_mat$bs_inc[data_mat$bs_inc<0]<-0.01 

    # while (data_mat$bs_inc[n]<0){
  #   bs_res<-sample(data_mat$scaled_pearson, size=1, replace = T)
  #   data_mat$bs_inc[n]<-(data_mat$mu[n]+
  #                          bs_res*sqrt(data_mat$mu[n]*data_mat$theta[n]))
  # }
  # while (data_mat$bs_inc[9]<0){
  #   bs_res<-sample(data_mat$scaled_pearson, size=1, replace = T)
  #   data_mat$bs_inc[9]<-(data_mat$mu[9]+
  #                          bs_res*sqrt(data_mat$mu[9]*data_mat$theta[9]))
  # }
  # while (data_mat$bs_inc[n-1]+data_mat$bs_inc[n-2]<0){
  #   bs_res<-sample(data_mat$scaled_pearson, size=2, replace = T)
  #   data_mat$bs_inc[c(n-2,n-1)]<-(data_mat$mu[c(n-2,n-1)]+
  #                          bs_res*sqrt(data_mat$mu[c(n-2,n-1)]*
  #                                        data_mat$theta[c(n-2,n-1)]))
  # }
  # while (data_mat$bs_inc[8]+data_mat$bs_inc[17]<0){
  #   bs_res<-sample(data_mat$scaled_pearson, size=2, replace = T)
  #   data_mat$bs_inc[c(8,17)]<-(data_mat$mu[c(8,17)]+
  #                                   bs_res*sqrt(data_mat$mu[c(8,17)]*
  #                                                 data_mat$theta[c(8,17)]))
  # }
  data.frame(acc=data_mat$acc,dev=data_mat$dev,
             inc=data_mat$bs_inc,j=data_mat$dev)
}

BS_pred<-function(I,disp,method,alpha,beta,seed){
  pred_acc<-data.frame(acc=rep(NA,I*(I-1)/2+1),
                       dev=rep(NA,I*(I-1)/2+1),
                       j=rep(NA,I*(I-1)/2+1))
  pred_acc[1,]<-c(1,1,1)
  k=2
  for (i in 2:I){
    for (j in (I-i+2):I){
      pred_acc$acc[k]<-i
      pred_acc$dev[k]<-j
      k<-k+1
    }
  }
  pred_acc$j<-pred_acc$dev
  pred_acc$acc<-factor(pred_acc$acc)
  pred_acc$dev<-factor(pred_acc$dev)
  pred_TT<-model.matrix(~ acc + dev, data=pred_acc)
  if (method!="ev") {
    pred_SS<-model.matrix(as.formula(disp),data=pred_acc)
    pred_acc$theta<-exp(pred_SS%*%beta)
  }
  if (method=="ev") {
    for (i in 1:nrow(pred_acc)){
      pred_acc$theta[i]<-beta[data_mat$dev==pred_acc$dev[i]][1]
    }
  }
  pred_acc$mu<-exp(pred_TT%*%alpha)
  pred_acc$inc<-NULL
  set.seed(seed)
  pred_acc$inc<-pred_acc$theta*
        rpois(nrow(pred_acc),pred_acc$mu/pred_acc$theta)
  pred_acc[-1,c("acc","dev","inc")]
}




