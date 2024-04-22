##fixed point iteration about unidimensional covariance
fix_point_iter_unidime<-function(V,t_max,w){
  #V denotes the list containing all subsets of covariance
  #t represents the max number of iteration
  #w donotes weight
  K<-length(V)
  detla<-as.list(rep(0,t_max))
  Z<-as.list(rep(0,t_max))
  dis<-rep(0,t_max)
  detla[[1]]<-1
  t=1
  while(t<t_max){
    B<-Reduce("+",lapply(V,function(x) w*(sqrt(detla[[t]])*x*sqrt(detla[[t]]))^(1/2)))
    detla[[t+1]]<-1/sqrt(detla[[t]])*(B)^2*(1/sqrt(detla[[t]]))
    dis[t]<-abs(detla[[t+1]]-detla[[t]])
    if(dis[t]<0.000001) break else t=t+1
  }
  return(list(t=t,detla=detla,dis=dis))
}

##fixed point iteration about multidimensional covariance
fix_point_iter_matrix<-function(V,t_max,w){
  #V denotes the list containing all subsets of covariance
  #t represents the max number of iteration
  #w donotes weight
  d<-nrow(V[[1]])
  #w<-1/K
  detla<-as.list(rep(0,t_max))
  dis<-rep(0,t_max)
  detla[[1]]<-diag(1,d)
  t=1
  while(t<t_max){
    sqrt_detla<-half(detla[[t]])
    B<-Reduce("+",lapply(V,function(x) w*half(detla[[t]]%*%x)))
    detla[[t+1]]<-solve(sqrt_detla)%*%B%*%t(B)%*%solve(sqrt_detla)
    dis[t]<-sum(diag(abs(detla[[t+1]]-detla[[t]])))
    if(dis[t]<0.000001) break else t=t+1
  }
  return(list(t=t,detla=detla,dis=dis[t]))
}


###################################Simulation###################################

##polya-gamma algorithm for subsets
library(truncnorm)
library(BayesLogit)
PG2PL_group_Par<-function(x){
  t1<-Sys.time()
  u<-y_subset[[x]]
  N<-nrow(u)
  J<-ncol(u)
  Rtheta<-matrix(0,N,ntime)
  Ra<-matrix(0,ntime,J)
  Rb<-matrix(0,ntime,J)
  R_xi<-as.list(rep(0,ntime))
  theta0<-matrix(0,nrow=N,ncol=1)
  a0<-matrix(1,nrow=1,ncol=J)
  b0<-matrix(0,nrow=1,ncol=J)
  
  M_theta0<-theta0%*%matrix(1,nrow=1,ncol=J)
  M_a0<-matrix(1,nrow=N,ncol=1)%*%a0
  M_b0<-matrix(1,nrow=N,ncol=1)%*%b0
  for(i in 1:ntime){
    temp<-M_a0*abs(M_theta0-M_b0)
    w0<-matrix(rpg(N*J,1,temp),N,J,byrow = F)
    
    ##theta, prior N(0,1)
    miu_theta<-0
    sig_theta<-1
    V_theta<-1/(rowSums(M_a0^2*w0)+1/sig_theta)
    m_theta<-V_theta*(rowSums(M_a0^2*w0*M_b0+M_a0*(u-1/2))+miu_theta/sig_theta)
    theta0<-rnorm(N,m_theta,sqrt(V_theta))
    Rtheta[,i]<-theta0
    M_theta0<-theta0%*%matrix(1,nrow=1,ncol=J)
    
    ##a, prior N(0,10)
    miu_a<-0
    sig_a<-10
    V_a<-1/(k*colSums((M_theta0-M_b0)^2*w0)+1/sig_a)
    m_a<-V_a*(k*colSums((M_theta0-M_b0)*(u-1/2))+miu_a/sig_a)
    a0<-rtruncnorm(J,a=0,b=Inf,mean=m_a,sd=sqrt(V_a))
    Ra[i,]<-a0
    M_a0<-matrix(1,nrow=N,ncol=1)%*%a0
    
    ##b, prior N(0,10)
    miu_b<-0
    sig_b<-10
    V_b<-1/(k*colSums(M_a0^2*w0)+1/sig_b)
    m_b<-V_b*(k*colSums(M_a0^2*w0*M_theta0-M_a0*(u-1/2))+miu_b/sig_b)
    b0<-rnorm(J,m_b,sqrt(V_b))
    Rb[i,]<-b0
    M_b0<-matrix(1,nrow=N,ncol=1)%*%b0
    R_xi[[i]]<-cbind(a0,b0)
  }
  t0<- Sys.time()
  time<-difftime(t0,t1,units="hours")
  ##return posterior draws
  return(list(Rtheta=Rtheta[,-(1:nburn)],R_xi=R_xi[-(1:nburn)],Ra=Ra[-(1:nburn),],Rb=Rb[-(1:nburn),],time=time))
}



##Data generation
set.seed(1818)
N<-5000
J<-20##item number
theta<-matrix(rnorm(N),ncol=1)
a<-matrix(rlnorm(J,0.3,0.2),nrow=1)
b<-matrix(rnorm(J),nrow=1)
#the matrix form of a, b and theta
M_a<-matrix(1,nrow=N,ncol=1)%*%a
M_b<-matrix(1,nrow=N,ncol=1)%*%b
M_theta<-theta%*%matrix(1,nrow=1,ncol=J)
p<-1/(1+exp(-(M_theta-M_b)*M_a))
set.seed(666)
u<-matrix(runif(N*J),nrow=N,ncol=J)
y<-ifelse(p>u,1,0)


###parallel process
library(parallel) 
cores <- detectCores() 
cl <- makeCluster(cores-2,type="PSOCK") #Initialize the number of cores to be used

set.seed(1818)
r_num<-sample(c(1:N),N)
k<-10##the number of subsets
s_num<-N/k ##subset size
##data partitioning
y_subset<-lapply(c(1:k),function(x) y[r_num[((x-1)*s_num+1):(x*s_num)],])
ntime<-10000
nburn<-5000
clusterExport(cl,c('y_subset','ntime','nburn','k'))
clusterEvalQ(cl, library(truncnorm))
clusterEvalQ(cl, library(BayesLogit))

##Starting parallel processing for each subset
results_group <- parLapply(cl=cl,1:k,PG2PL_group_Par) # apply的并行版本
stopCluster(cl) #Turn off parallel mode

##computing A^{1/2}
half<-function(X){
  Z<-eigen(X)
  U<-Z$vectors
  D<-Z$values
  sqrt_X<-U%*%diag(sqrt(D))%*%solve(U)
  return(sqrt_X)
}

#WASP Integration Operations
group_temp<-lapply(results_group,function(x){
  Rtheta<-apply(x$Rtheta,1,mean)
  Ra<-t(x$Ra)
  m_a<-apply(Ra,1,mean)
  V_a<-apply((Ra-m_a%*%t(rep(1,(ntime-nburn))))^2,1,mean)
  U_a<-(Ra-m_a%*%t(rep(1,(ntime-nburn))))/sqrt(V_a)
  
  Rb<-t(x$Rb)
  m_b<-apply(Rb,1,mean)
  V_b<-apply((Rb-m_b%*%t(rep(1,(ntime-nburn))))^2,1,mean)
  U_b<-(Rb-m_b%*%t(rep(1,(ntime-nburn))))/sqrt(V_b)
  
  Rxi<-x$R_xi
  m_xi<-Reduce("+",Rxi)/length(Rxi)
  V_xi1<-lapply(Rxi,function(r){
    temp1<-r-m_xi
    Sig<-as.list(rep(0,J))
    for (j in 1:J) {
      Sig[[j]]<-temp1[j,]%*%t(temp1[j,])
    }
    return(Sig)
  })
  V_xi<-as.list(rep(0,J))
  U_xi<-as.list(rep(0,J))
  for(j in 1:20){
    V_xi[[j]]<-Reduce("+",lapply(V_xi1,function(rr) rr[[j]]))/length(V_xi1)
    U_xi[[j]]<-lapply(Rxi,function(r) half(solve(V_xi[[j]]))%*%(r[j,]-m_xi[j,]))
  }
  return(list(Rtheta=Rtheta,m_a=m_a,m_b=m_b,m_xi=m_xi,V_a=V_a,U_a=U_a,V_b=V_b,U_b=U_b,V_xi=V_xi,U_xi=U_xi)) 
})



##fixed point iteration about the simulation of unidimensional covariance
#Take the discrimination parameter a as an example.  20 items represent 20 replications.
t_max<-10
w<-1/k
results<-as.list(rep(0,J))
for(j in 1:J){
  V_a<-unlist(lapply(group_temp,function(x) x$V_a[j]))
  results[[j]]<-fix_point_iter_unidime(V_a,t_max,w)
}
t<-unlist(lapply(results,function(x) x$t))
dis<-unlist(lapply(results,function(x) x$dis[2]))
mean(dis)


##fixed point iteration about the simulation of multidimensional covariance
V_xi2<-lapply(group_temp,function(x) x$V_xi)
results1<-as.list(rep(0,J))
for (j in 1:J) {
  V_xi_bar1<-lapply(V_xi2,function(x) x[[j]])
  results1[[j]]<-fix_point_iter_matrix(V_xi_bar1,t_max,w)
}
t<-unlist(lapply(results1,function(x) x$t))
t
mean(t)
dis<-unlist(lapply(results1,function(x) x$dis))
mean(dis)
