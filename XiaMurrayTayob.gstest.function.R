

######################################################################################
#Function of calculate the proposed test statistics with choices of safety bounds
#Function of calculate the restricted mean survival time (RMST) with follow-up windows
######################################################################################
#X is a vector of observed time
#delta is a vector of failure indicator, alive=0, dead=1
#E is a vector of entry time
#group is a vector of group indicator = 1 or 2
#s is a vector of analysis time
#Tau is a scalar of the length of follow-up windows
#safety.bound has choice of "OF", "Pocock", "JT"
#overall.alpha.safety.JT is only used when safety.bound="JT"


XiaMurrayTayob.gstest <- function(X, delta, E, group, s, Tau, 
                                  alpha.efficacy = 0.025, 
                                  alpha.safety = 0.025, 
                                  safety.bound = "JT", 
                                  overall.alpha.safety.JT = 0.2){
  
  require(survival)
  require(gdata)
  require(mvtnorm)
  require(MASS)
  
  #Check for errors of inputs
  if(sum(X<0)>0)
  {
    stop("Error in follow-up times input")
  }
  if(sum(delta<0 | delta>1 | floor(delta)!=delta)>0)
  {
    stop("Error in failure indicator input")
  }   
  if(Tau>max(X))
  {
    stop("Upper limit of integration greater than largest observed follow-up time")
  }
  if(Tau>min(s))
  {
    stop("Upper limit of integration greater than the first analysis time")
  }
  if(Tau<0)
  {
    stop("Invalid upper limit of integration")
  }
  
  
  #Get RMST with windows
  newRMST.summary1 <- RMST.windows(X[group==1],delta[group==1],E[group==1],s,Tau)
  newRMST.summary2 <- RMST.windows(X[group==2],delta[group==2],E[group==2],s,Tau)
  n1 <- newRMST.summary1$NumberEntered
  n2 <- newRMST.summary2$NumberEntered
  norm.const <- sqrt(n1*n2/(n1+n2))
  
  #Calculate test statistics T
  Tstat <- norm.const*(newRMST.summary1$RmeanST-newRMST.summary2$RmeanST)
  Cov.Tstat <- norm.const%*%t(norm.const)*(newRMST.summary1$cov+newRMST.summary2$cov)
  #Standardize
  Tstat.std <- Tstat/sqrt(diag(Cov.Tstat))
  Cor.Tstat <- cov2cor(Cov.Tstat)
  
  #Get Z
  Z.T <- rmvnorm(10000000,rep(0,length(s)),Cor.Tstat)
  
  #Identify critical values
  v <- s/max(s)
  #Efficacy-OF
  cum.alpha.effi <- 1-pnorm(qnorm(1-alpha.efficacy)/sqrt(v))
  alpha.each.effi <- c(cum.alpha.effi[1],(cum.alpha.effi[-1]-cum.alpha.effi[-5])/(1-cum.alpha.effi[-5]))
  
  #Safety
  if (safety.bound == "OF"){ cum.alpha.safe <- 1-pnorm(qnorm(1-alpha.safety)/sqrt(v)) }
  if (safety.bound == "Pocock"){ cum.alpha.safe <- alpha.safety*log(1+(exp(1)-1)*v) }
  if (safety.bound == "JT"){
    w <- log(alpha.safety/overall.alpha.safety.JT)/log(v[1])
    #Overall 1-sided incorrect stopping under H0 caused by trial safety boundary=overall.alpha.safety.JT
    #1-sided incorrect stopping under H0 caused by trial safety boundary at 1st look=alpha.safety
    #projected % information at 1st look=v[1]
    cum.alpha.safe <- overall.alpha.safety.JT*v^w
  }
  alpha.each.safe <- c(cum.alpha.safe[1],(cum.alpha.safe[-1]-cum.alpha.safe[-5])/(1-cum.alpha.safe[-5]))
  
  #Critical values
  critical.effi <- Find.critical.Z(Z.T, alpha.each.effi)
  critical.safe <- Find.critical.Z(Z.T, alpha.each.safe)
  stop.effi <- as.numeric(Tstat.std >= critical.effi)
  stop.safe <- as.numeric(Tstat.std <= -critical.safe)
  
  return(list(standarized.test.statistics = Tstat.std,
              correlation.of.statistics = Cor.Tstat,
              efficacy.critical.value = critical.effi,
              efficacy.stop = stop.effi,
              safety.critical.value = critical.safe,
              safety.stop = stop.safe))
}


#Function of calculate RMST with follow-up windows
RMST.windows=function(X, delta, E, s, Tau){
  
  #remove any subjects with missing data
  data <- data.frame(X,delta,E)
  data <- na.omit(data)
  
  #Setup function needed later
  sum_function <- function(X)
  {
    temp <- array(X,c(length(X),length(X)))
    lowerTriangle(temp,diag = F) <- 0
    apply(temp,2,sum)
  }
  
  #Get some parameters
  t <- seq(from=0, to=max(s), by=Tau/2) #Start of follow-up windows
  b <- length(t) #Number of follow-up windows
  s.final <- max(s)
  
  #Do the calculation for each si
  n_s <- array()
  Rmean <- array()
  williams_var <- array()
  Big.Z <- matrix(NA,nrow(data),length(s))
  Big.Include <- matrix(NA,nrow(data),length(s))
  for (si in 1:length(s)){
    #Get observed data
    data$in.s1 <- ifelse(data$E<s[si],1,0)
    data$C.s1 <- ifelse(data$E<s[si],s[si]-data$E,0)
    data$X.s1 <- ifelse(data$X<=data$C.s1,data$X,data$C.s1)
    data$delta.s1 <- ifelse(data$X<=data$C.s1,data$delta,0)
    
    data$C.s2 <- ifelse(data$E<s.final,s.final-data$E,0)
    data$X.s2 <- ifelse(data$X<=data$C.s2,data$X,data$C.s2)
    data$delta.s2 <- ifelse(data$X<=data$C.s2,data$delta,0)
    
    #Define some alias
    X1 <- data[data$in.s1==1,]$X.s1
    delta1 <- data[data$in.s1==1,]$delta.s1
    include1 <- data$in.s1
    X2 <- data$X.s2
    delta2 <- data$delta.s2
    
    n1 <- sum(include1)#Num of patients enrolled at si
    n2 <- length(X2)#Num of patients enrolled at s.final
    
    #generate X and delta for each window
    X_tj1 <- array(NA,c(n1,b))
    delta_tj1 <- array(NA,c(n1,b))
    X_tj2 <- array(NA,c(n2,b))
    delta_tj2 <- array(NA,c(n2,b))
    
    for(j in 1:b)
    {
      X_j1 <- X1-t[j]
      delta_j1 <- delta1
      delta_j1[X_j1<0] <- 0 #terminating events observed before time t[k]->patient censored at 0
      X_j1[X_j1<0] <- 0
      X_tj1[,j] <- X_j1
      delta_tj1[,j] <- delta_j1
      
      X_j2 <- X2-t[j]
      delta_j2 <- delta2
      delta_j2[X_j2<0] <- 0 #terminating events observed before time t[k]->patient censored at 0
      X_j2[X_j2<0] <- 0
      X_tj2[,j] <- X_j2
      delta_tj2[,j] <- delta_j2
    }
    
    #At the earlier analysis time, si
    #Get observed event time T1,T2.....TM from data collected to s.final (use as much data as possible)
    observed_events <- sort(X_tj2*delta_tj2)
    T <- unique(c(observed_events[observed_events<=Tau],Tau))#T0=0, T_M+1=Tau
    M <- length(T)-1
    T_array <- round(t(array(T[1:M],c(M,b))),6)
    
    #Get dN, Y for each patient combining all windows for si
    dN_i <- array(NA,c(n1,M))
    Y_i <- array(NA,c(n1,M))
    for(i in 1:n1)
    {
      temp3 <- array(round(X_tj1[i,],6),c(b,M))
      temp4 <- array(delta_tj1[i,],c(b,M))
      dN_i[i,] <- apply((temp3==T_array & temp4==1),2,sum) 
      Y_i[i,] <- apply(temp3>=T_array,2,sum)
    }
    dN <- apply(dN_i,2,sum)#dN, Y combining all patients across all windows
    Y <- apply(Y_i,2,sum)
    
    #Get the three components of Zi
    #First, get z_i_3 taking advantage of data collected to s.final
    #Get the numerator of Z_i_3
    dN_i_s1_hat <- array(NA,c(n1,M))
    Y_i_s1 <- array(NA,c(n1,M))
    for(i in 1:n1)
    {
      k <- which(include1==1)[i]
      temp5 <- array(round(X_tj2[k,],6),c(b,M))
      temp6 <- array(delta_tj2[k,],c(b,M))
      dN_ij_s2 <- ifelse(temp5==T_array & temp6==1,1,0)
      Y_ij_s2 <- ifelse(temp5>=T_array,1,0)
      
      temp3 <- array(round(X_tj1[i,],6),c(b,M))
      Y_ij_s1 <- ifelse(temp3>=T_array,1,0)
      
      dN_ij_s1_hat <- Y_ij_s1*dN_ij_s2/(Y_ij_s2)
      dN_ij_s1_hat[dN_ij_s1_hat=="NaN"] <- 0
      
      dN_i_s1_hat[i,] <- apply(dN_ij_s1_hat,2,sum)
      Y_i_s1[i,] <- apply(Y_ij_s1,2,sum)
    }
    
    #Get the denominator of Z_i_3
    surv_prob_event_s2=array(0,c(b,M))
    for (j in 1:b)
    {
      sur_time <- summary(survfit(Surv(X_tj2[,j], delta_tj2[,j]==1)~1))$time#get K-M estimator for event
      sur_prob <- summary(survfit(Surv(X_tj2[,j], delta_tj2[,j]==1)~1))$surv
      for (m in 1:M){
        rank <- max(which(sort(c(T[m],sur_time))==T[m]))
        if (rank==1) {surv_prob_event_s2[j,m] <- 1}
        if (rank>1) {surv_prob_event_s2[j,m] <- sur_prob[rank-1]}
      }
    }
    
    surv_prob_censor_s1 <- array(0,c(b,M))
    for (j in 1:b)
    {
      sur_time <- summary(survfit(Surv(X_tj1[,j], delta_tj1[,j]==0)~1))$time#get K-M estimator for censor
      sur_prob <- summary(survfit(Surv(X_tj1[,j], delta_tj1[,j]==0)~1))$surv
      for (m in 1:M){
        rank <- max(which(sort(c(T[m],sur_time))==T[m]))
        if (rank==1) {surv_prob_censor_s1[j,m] <- 1}
        if (rank>1) {surv_prob_censor_s1[j,m] <- sur_prob[rank-1]}
      }
    }
    
    #z_i_3
    denom <- apply(surv_prob_event_s2*surv_prob_censor_s1,2,sum)
    
    z_i_3 <- t(t(dN_i_s1_hat-t(dN/Y*t(Y_i_s1)))/denom)
    z_i_3[z_i_3=="NaN"] <- 0
    z_i_3_sum <- t(apply(z_i_3,1,sum_function))                                                                              
    
    #Get the first two
    time_int <- T[2:(M+1)]-T[1:M]
    time_int_array <- t(array(time_int,c(M,n1)))
    
    S_hat <- exp(-sum_function(dN/Y))
    S_hat[S_hat=="NaN"] <- 0
    S_hat_array <- t(array(S_hat,c(M,n1)))
    
    #Get product and sum, a vector of Zi with length n1 at si
    z_i_s1_Tau <- apply(time_int_array*S_hat_array*z_i_3_sum,1,sum)
    n_s[si] <- sum(include1)
    Rmean[si] <- sum(time_int*S_hat)
    williams_var[si] <- var(z_i_s1_Tau)/n1 #empirical variance estimate of estimate of overall tau restricted mean survival
    #Save Zi for covariance
    Big.Z[which(include1==1),si] <- z_i_s1_Tau
    Big.Include[,si] <- include1
  }
  
  #Get covariance
  empirical.cov <- matrix(NA,length(s),length(s))
  for (c1 in (1:length(s))){
    for (c2 in (1:length(s))){
      include <- Big.Include[,c1]
      if (sum(Big.Include[,c1])>sum(Big.Include[,c2])) {include <- Big.Include[,c2]}
      nd <- max(sum(Big.Include[,c1]),sum(Big.Include[,c2]))
      empirical.cov[c1,c2] <- cov(Big.Z[which(include==1),c1],Big.Z[which(include==1),c2])/nd
    }
  }
  
  list(NumberEntered=n_s,RmeanST=Rmean,Var=williams_var,cov=empirical.cov)
}


#Function of finding critical values
Find.critical.Z <- function(Z,alpha.each){
  Z <- abs(Z)
  critical.effi.Z <- array()
  cr=1
  repeat{
    if (cr==1){
      critical.effi.Z[cr] <- quantile(Z[,cr],1-2*alpha.each[cr])
      AcceptZ <- Z[Z[,cr]<=critical.effi.Z[cr],]
    }
    if (cr>1){
      critical.effi.Z[cr] <- quantile(AcceptZ[,cr],1-2*alpha.each[cr])
      AcceptZ <- AcceptZ[AcceptZ[,cr]<=critical.effi.Z[cr],]
    }
    cr=cr+1
    if (cr==length(Z[1,])+1){break}
  }
  return (critical.effi.Z)
}




#EXAMPLE
#Read data
data <- read.csv("example_data.csv")
#Run the test
XiaMurrayTayob.gstest(X=data$X, 
                      delta=data$delta, 
                      E=data$E, 
                      group=data$group, 
                      s=1:5,
                      Tau=1, 
                      safety.bound="JT")

