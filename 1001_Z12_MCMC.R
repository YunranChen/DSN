library(GPfit)
library(MASS)
library(dplyr)
library(purrr)
library(tidyr)
library(reshape2)
library("forcats")
library(mvtnorm)
library(rbenchmark)
library(microbenchmark)
library("coda")
library(MCMCpack)
rm(list=ls())

###read the params

# If you are running multiple runs in compute cluster:
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID', 0))
pars=scan(file = '1001_pars.txt',what = "",skip = TASKID, nlines = 1)
# pars: Y array, covariates X, Hstar,rho_mu,rho_ab,k,a1, a2.
#pars=c("0927Y_arr_42.RData","0927Z12t_top18_2006Q4_2016Q3.RData",5,0.5,0.5,0.1,2,2,"EU","Quarter")
#prepare

# Naive sampler for PG(1, z)
# (based on the finite approximation of infinite sum)
# can replace this using rpg function from the package BayesLogit
rpg_naive = function(z, n_terms = 100){
  g = rexp(n_terms, 1)
  out = 1 / (2*(pi^2)) * sum(g / ((1:n_terms - 1/2)^2 + z^2 / (4*(pi^2))))
  return(out)
}

## Read data

# read Y array
load(pars[1])
# read X
load(pars[2])
# number of nodes
V=dim(Y_arr)[1]
# number of time points
N=dim(Y_arr)[3]

Y_arr[,,N]=NA #Y_{40} a matrix of missing values


## Inference
H_star=pars[3]%>%as.numeric(.)
rho_x0=pars[4]%>%as.numeric(.)
rho_ab0=pars[5]%>%as.numeric(.)
k_u0=k_x0=k_beta01=k_beta02=k_ab0=pars[6]%>%as.numeric(.)
a1=pars[7]%>%as.numeric(.)
a2=pars[8]%>%as.numeric(.)

niter=5000*20
burnin=1000
set.seed(456)


# MCMC start here


now=proc.time()
logFile =paste0(TASKID,"1001log_file.txt")
cat(paste0(0,",",now[3]), file=logFile, append=FALSE, sep = "\n")


##prior
c_u0=corr_matrix(X = 1:N,beta = log10(k_u0),corr = list(type="exponential",power=2))
K_mu_inv=chol2inv(chol(c_u0+diag(rep(1e-8,N)))) #inverse
c_x0=corr_matrix(X = 1:N,beta = log10(k_x0),corr = list(type="exponential",power=2))
K_x0_inv=chol2inv(chol(c_x0+diag(rep(1e-8,N))))
x_x0=matrix(c(1,rho_x0,rho_x0,1),nrow=2)
x_x0_inv=chol2inv(chol(x_x0+diag(rep(1e-8,2))))
C_x0=x_x0%x%c_x0
K_x_inv=x_x0_inv%x%K_x0_inv

c_ab0=corr_matrix(X = 1:N,beta = log10(k_ab0),corr = list(type="exponential",power=2))
K_ab0_inv=chol2inv(chol(c_ab0+diag(rep(1e-8,N))))
x_ab0=matrix(c(1,rho_ab0,rho_ab0,1),nrow=2)
x_ab0_inv=chol2inv(chol(x_ab0+diag(rep(1e-8,2))))
C_ab0=x_ab0%x%c_ab0
K_ab_inv=x_ab0_inv%x%K_ab0_inv

c_beta01=corr_matrix(X = 1:N,beta = log10(k_beta01),corr = list(type="exponential",power=2))
K_beta_inv1=chol2inv(chol(c_beta01+diag(rep(1e-8,N))))
c_beta02=corr_matrix(X = 1:N,beta = log10(k_beta02),corr = list(type="exponential",power=2))
K_beta_inv2=chol2inv(chol(c_beta02+diag(rep(1e-8,N))))

mu0=mvrnorm(n=1,mu=rep(0,N),Sigma = c_u0)
beta01=mvrnorm(n=1,mu=rep(0,N),Sigma = c_beta01)
beta02=mvrnorm(n=1,mu=rep(0,N),Sigma = c_beta02)

ab0=mvrnorm(n=V,mu=rep(0,2*N),Sigma = C_ab0)
a0=ab0[,1:N]
b0=ab0[,(N+1):(2*N)]

vv=rep(1,H_star)
tao=accumulate(vv,prod) 
X_arr0=array(dim=c(V,H_star,2*N))

for(h in 1:H_star){
  X_arr0[,h,]=mvrnorm(n=V,mu=rep(0,2*N),Sigma = (1/tao[h])*C_x0)
}

Xs_arr0=X_arr0[,,1:N]
Xr_arr0=X_arr0[,,(N+1):(2*N)]
##Set Y_arr[,,N=40] as missing value
pi_arr0=map(1:N,~1/(1+exp(-(Xs_arr0[,,.x]%*%t(Xr_arr0[,,.x])+mu0[.x]+Z1[,,.x]*beta01[.x]+Z2[,,.x]*beta02[.x]+a0[,.x]+matrix(b0[,.x],nrow=V,ncol=V,byrow=T)))))
for (i in 1:V){
  for (j in 1:V){
    Y_arr[i,j,N]=sample(x = c(1,0),size = 1,prob = c(pi_arr0[[N]][i,j],1-pi_arr0[[N]][i,j]))
  }
}
diag(Y_arr[,,N])=NA

#Store the result: create empty bags

#W_cache=matrix(nrow=V*V*N,ncol=niter)
mu_cache=matrix(nrow=N,ncol=niter)
beta_cache1=matrix(nrow=N,ncol=niter)
beta_cache2=matrix(nrow=N,ncol=niter)
#X_cache=matrix(nrow=V*H_star*N,ncol=niter)
tao_cache=matrix(nrow=H_star,ncol=niter)
pi_arr_est_cache=matrix(nrow=V*V*N,ncol=niter)
a_cache=array(dim=c(V,N,niter))
b_cache=array(dim=c(V,N,niter))
Y_arrN_cache=array(dim=c(V,V,niter))
##posterior
for (iter in 1:niter){
  
  #sample W:using Xs_arr0,Xr_arr0,mu0,Z1,beta01,beta02
  S=map(1:N,~Xs_arr0[,,.x]%*%t(Xr_arr0[,,.x])) #S=Xt(X)
  W=array(data = NA,dim=c(V,V,N))
  for (t in 1:N){
    W[,,t]=map_dbl(S[[t]]+mu0[t]+Z1[,,t]*beta01[t]+Z2[,,t]*beta02[t]+a0[,t]+matrix(b0[,t],nrow=V,ncol=V,byrow=T),~rpg_naive(z = .x,n_terms = 100))%>%matrix(data = .,nrow = V)
    diag(W[,,t])=NA
  }
  
  #sample mu0:using W,S(Xs_arr0,Xr_arr0),Z1,beta01,beta02
  Sigma_mu0=(apply(X = W,MARGIN = 3,FUN = sum,na.rm=TRUE))%>%diag(.)
  Sigma_mu=chol2inv(chol(Sigma_mu0+K_mu_inv+diag(rep(1e-8,N))))
  mu_mu=Sigma_mu%*%map_dbl(1:N,~sum(Y_arr[,,.x]-0.5-W[,,.x]*(S[[.x]]+Z1[,,.x]*beta01[.x]+Z2[,,.x]*beta02[.x]+a0[,.x]+matrix(b0[,.x],nrow=V,ncol=V,byrow=T)),na.rm = TRUE)) #W-diag NA
  mu0=mvrnorm(1,mu_mu,Sigma_mu)
  mu_cache[,iter]=mu0
  
  #sample beta0:using Z1,W,Y_arr,S
  Sigma_beta01=(apply(X = Z1^2*W,MARGIN = 3,FUN = sum,na.rm=TRUE))%>%diag(.)
  Sigma_beta1=chol2inv(chol(Sigma_beta01+K_beta_inv1+diag(rep(1e-8,N))))
  beta_mu1=Sigma_beta1%*%map_dbl(1:N,~sum(Z1[,,.x]*(Y_arr[,,.x]-0.5-W[,,.x]*(mu0[.x]+S[[.x]]+Z2[,,.x]*beta02[.x]+a0[,.x]+matrix(b0[,.x],nrow=V,ncol=V,byrow=T))),na.rm = TRUE)) #W-diag NA
  beta01=mvrnorm(1,beta_mu1,Sigma_beta1)
  beta_cache1[,iter]=beta01
  
  Sigma_beta02=(apply(X = Z2^2*W,MARGIN = 3,FUN = sum,na.rm=TRUE))%>%diag(.)
  Sigma_beta2=chol2inv(chol(Sigma_beta02+K_beta_inv2+diag(rep(1e-8,N))))
  beta_mu2=Sigma_beta2%*%map_dbl(1:N,~sum(Z2[,,.x]*(Y_arr[,,.x]-0.5-W[,,.x]*(mu0[.x]+S[[.x]]+Z1[,,.x]*beta01[.x]+a0[,.x]+matrix(b0[,.x],nrow=V,ncol=V,byrow=T))),na.rm = TRUE)) #W-diag NA
  beta02=mvrnorm(1,beta_mu2,Sigma_beta2)
  beta_cache2[,iter]=beta02
  
  #sample X_arr0(Xs_arr0;Xr_arr0):using tao,W,Xs_arr0,Xr_arr0,mu0,beta0,Z1
  
  prior_sig_x=x_x0_inv%x%(diag(tao)%x%K_x0_inv)
  for (v in 1:V){
    X_tilta=matrix(0,nrow = 2*(V-1)*N,ncol = 2*H_star*N)
    w1=W[v,-v,]%>%t(.)%>%as.vector(.)
    w2=W[-v,v,]%>%t(.)%>%as.vector(.)
    w=c(w1,w2)
    Omega=diag(w)
    for (t in 1:N){
      X_tilta[seq(from=t,to=(V-1)*N,by=N),seq(from=t,to=H_star*N,by=N)]=Xr_arr0[-v,,t]
      X_tilta[seq(from=t+(V-1)*N,to=(V-1)*N*2,by=N),seq(from=t+H_star*N,to=H_star*N*2,by=N)]=Xs_arr0[-v,,t]
    }
    Sigma_x0=t(X_tilta)%*%Omega%*%X_tilta+prior_sig_x
    Sigma_x=chol2inv(chol(Sigma_x0+diag(rep(1e-8,N*H_star*2))))
    y1=Y_arr[v,-v,]%>%t(.)%>%as.vector(.)
    y2=Y_arr[-v,v,]%>%t(.)%>%as.vector(.)
    y=c(y1,y2)
    z11=Z1[v,-v,]%>%t(.)%>%as.vector(.)
    z21=Z1[-v,v,]%>%t(.)%>%as.vector(.)
    z_1=c(z11,z21)
    z12=Z2[v,-v,]%>%t(.)%>%as.vector(.)
    z22=Z2[-v,v,]%>%t(.)%>%as.vector(.)
    z_2=c(z12,z22)
    x_av=c(rep(a0[v,],(V-1)),a0[-v,]%>%t(.)%>%as.vector())
    x_bv=c(b0[-v,]%>%t(.)%>%as.vector(),rep(b0[v,],(V-1)))
    mu_x0=t(X_tilta)%*%(y-0.5-w*(rep(mu0,2*(V-1))+z_1*rep(beta01,2*(V-1))+z_2*rep(beta02,2*(V-1))+x_av+x_bv))
    mu_x=Sigma_x%*%mu_x0
    X_v=rmvnorm(n = 1,mean = mu_x,sigma = Sigma_x)
    Xs_arr0[v,,]=matrix(data = X_v[1:(H_star*N)],nrow=H_star,byrow = TRUE)
    Xr_arr0[v,,]=matrix(data = X_v[(H_star*N+1):(2*H_star*N)],nrow=H_star,byrow = TRUE)
  }
  
  #X_cache[,iter]=X_arr0%>%as.vector(),
  
  #sample tao,v:using Xr_arr0,Xs_arr0,vv
  xKx=map_dbl(.x = 1:H_star,.f = function(l){map_dbl(.x = 1:V,.f = ~ t(c(Xs_arr0[.x,l,],Xr_arr0[.x,l,]))%*%K_x_inv%*%c(Xs_arr0[.x,l,],Xr_arr0[.x,l,]))%>%sum(.)})
  
  tao_1=c(1,vv[-1])%>%accumulate(.,prod)
  rate1=1+0.5*sum(tao_1*xKx)
  v1=rgamma(n = 1,shape = a1+V*N*H_star,rate = rate1)
  rate_l=vector(length = H_star-1)
  for (h in 2:H_star){
    rate_l[h-1]=map_dbl(h:H_star,~prod(vv[1:.x][-h])*xKx[.x])%>%sum(.)*0.5+1
  }
  vl=map_dbl(2:H_star,~rgamma(n = 1,shape = a2+V*N*(H_star-.x+1),rate = rate_l[.x-1]))
  vv=c(v1,vl)
  tao=accumulate(vv,prod)
  
  tao_cache[,iter]=tao
  
  #update ab0
  for (v in 1:V){
    w_ab1=W[v,-v,]%>%apply(.,2,sum)
    w_ab2=W[-v,v,]%>%apply(.,2,sum)
    Sigma_abv0=diag(c(w_ab1,w_ab2))
    Sigma_abv=chol2inv(chol(Sigma_abv0+K_ab_inv+diag(rep(1e-8,2*N))))
    abmu1=map_dbl(1:N,~sum(Y_arr[v,-v,.x]-0.5-W[v,-v,.x]*(mu0[.x]+S[[.x]][v,-v]+Z1[v,-v,.x]*beta01[.x]+Z2[v,-v,.x]*beta02[.x]+b0[-v,.x]),na.rm = TRUE))
    abmu2=map_dbl(1:N,~sum(Y_arr[-v,v,.x]-0.5-W[-v,v,.x]*(mu0[.x]+S[[.x]][-v,v]+Z1[-v,v,.x]*beta01[.x]+Z2[-v,v,.x]*beta02[.x]+a0[-v,.x]),na.rm = TRUE))
    mu_abv=Sigma_abv%*%c(abmu1,abmu2)
    ab0[v,]=mvrnorm(n=1,mu=mu_abv,Sigma = Sigma_abv)
    a0[v,]=ab0[v,1:N]
    b0[v,]=ab0[v,(N+1):(2*N)]
    a_cache[v,,iter]=a0[v,]
    b_cache[v,,iter]=b0[v,]
  }
  
  #calculate pi!!!Finally!!!:using X_arr0,mu0,
  pi_arr0=map(1:N,~1/(1+exp(-(Xs_arr0[,,.x]%*%t(Xr_arr0[,,.x])+mu0[.x]+Z1[,,.x]*beta01[.x]+Z2[,,.x]*beta02[.x]+a0[,.x]+matrix(b0[,.x],nrow=V,ncol=V,byrow=T)))))
  pi_vec=pi_arr0%>%unlist(.)
  pi_arr_est_cache[,iter]=pi_vec
  for (i in 1:V){
    for (j in 1:V){
      Y_arr[i,j,N]=sample(x = c(1,0),size = 1,prob = c(pi_arr0[[N]][i,j],1-pi_arr0[[N]][i,j]))
    }
  }
  diag(Y_arr[,,N])=NA
  Y_arrN_cache[,,iter]=Y_arr[,,N]
  # track the number of iterations
if (iter%%100==0){
  cat(paste0(iter,",",proc.time()[3]), file=logFile, append=TRUE, sep = "\n")
}
  # save file every 5000 iterations
  if (iter%%5000==0){
    save.image(paste0("1001_TASKID",TASKID,"_code",pars[9],"_iter",iter,"_a1_",a1,"_Hstar",H_star,"_rho_x0",rho_x0,"_rho_ab0",rho_ab0,"_k0_",k_x0,".RData"),compress = "xz")
  }
}

# compute time
future=proc.time()
tictoc=future-now

# save files
save.image(paste0("1001_TASKID",TASKID,"_code",pars[9],"_iter",iter,"_a1_",a1,"_Hstar",H_star,"_rho_x0",rho_x0,"_rho_ab0",rho_ab0,"_k0_",k_x0,".RData"),compress = "xz")
