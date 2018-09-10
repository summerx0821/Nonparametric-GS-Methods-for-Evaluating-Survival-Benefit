
######################################################################################
#Function of calculate the closed form covariance of Zi(mu(s,tau)) in Web Appendix C
#We’ve used this closed form variance as a method to double-check that R code for
#our empirically calculated covariance is on target
######################################################################################

#Set up the scenario before running the function
A=NULL #Accrual time
s1=NULL #The 1st analysis time
s2=NULL #The 2nd analysis time
Tau=NULL #Time length of each window
lambda=NULL #hazard of events


#Define the function
get_closed_cov=function()
{
  #construct the windows
  t=seq(from=0, to=s2, by=Tau/2) 
  b=length(t)
  
  #Get the censoring "Survival" function
  minimum1=function(x){min(1,x)}
  survival_censor=function(u,s){  #half E~uniform(0, min(A,s1)), half E=0
    x=(s-u)/min(A,s)*ifelse(u<s,1,0)
    apply(x,c(1,2),minimum1)
  }
  
  #Create small interval for integration of u
  K=1000
  a_i=seq(from=0,to=Tau,length.out=K+1)
  d_i=c(a_i,Tau)-c(0,a_i) #length k+2
  b_i=(d_i[1:(K+1)]+d_i[2:(K+2)])/2 #length k+1, sum 1
  
  #Bottom-at risk process
  uv_plus_t=array(a_i,c(K+1,b))+t(array(t,c(b,K+1)))
  
  at_risk_u_t1j=exp(-lambda*uv_plus_t)*survival_censor(uv_plus_t,s1)
  at_risk_u_t1j[1,]=1 #By definition X(t_j)>=0 for all j
  at_risk_u=apply(at_risk_u_t1j,1,sum)
  
  at_risk_v_t2k=exp(-lambda*uv_plus_t)*survival_censor(uv_plus_t,s2)
  at_risk_v_t2k[1,]=1 #By definition X(t_k)>=0 for all k
  at_risk_v=apply(at_risk_v_t2k,1,sum)
  
  #Covariance on top
  cov_u_v=array(0,c(K+1,K+1))
  for(j in 1:b)
  {
    for(k in 1:b)
    {
      t1_j=t[j]
      t2_k=t[k]
      
      indicator1=array(as.numeric(array(uv_plus_t[,j],c(K+1,K+1))==t(array(uv_plus_t[,k],c(K+1,K+1)))),c(K+1,K+1))
      indicator2=array(as.numeric(array(uv_plus_t[,j],c(K+1,K+1))<=t(array(uv_plus_t[,k],c(K+1,K+1)))),c(K+1,K+1)) + array(as.numeric(array(a_i,c(K+1,K+1))==0),c(K+1,K+1))*t(array(as.numeric(array(a_i,c(K+1,K+1))<(t1_j-t2_k)),c(K+1,K+1)))
      indicator3=array(as.numeric(array(uv_plus_t[,j],c(K+1,K+1))>=t(array(uv_plus_t[,k],c(K+1,K+1)))),c(K+1,K+1)) + t(array(as.numeric(array(a_i,c(K+1,K+1))==0),c(K+1,K+1)))*array(as.numeric(array(a_i,c(K+1,K+1))<(t2_k-t1_j)),c(K+1,K+1))
      
      array1=array(lambda*at_risk_u_t1j[,j],c(K+1,K+1))*indicator1
      
      max_array=array(uv_plus_t[,j],c(K+1,K+1))*array(as.numeric(array(uv_plus_t[,j],c(K+1,K+1))>t(array(uv_plus_t[,k],c(K+1,K+1)))),c(K+1,K+1)) + t(array(uv_plus_t[,k],c(K+1,K+1)))*array(as.numeric(array(uv_plus_t[,j],c(K+1,K+1))<=t(array(uv_plus_t[,k],c(K+1,K+1)))),c(K+1,K+1))
      min_array=array(s1-uv_plus_t[,j],c(K+1,K+1))*array(as.numeric(array(s1-uv_plus_t[,j],c(K+1,K+1))<=t(array(s2-uv_plus_t[,k],c(K+1,K+1)))),c(K+1,K+1)) + t(array(s2-uv_plus_t[,k],c(K+1,K+1)))*array(as.numeric(array(s1-uv_plus_t[,j],c(K+1,K+1))>t(array(s2-uv_plus_t[,k],c(K+1,K+1)))),c(K+1,K+1))
      both_at_risk=exp(-lambda*max_array)*apply(min_array/min(A,s1)*ifelse(uv_plus_t[,j]<s1&uv_plus_t[,k]<s2,1,0),c(1,2),minimum1)
      
      both_at_risk[,1]=at_risk_u_t1j[,j]
      both_at_risk[1,]=at_risk_v_t2k[,k]
      
      array2=array(lambda^2,c(K+1,K+1))*both_at_risk*indicator2
      array3=t(array(lambda^2,c(K+1,K+1)))*both_at_risk*indicator3
      array4=array(lambda^2,c(K+1,K+1))*both_at_risk
      
      numerator=array(b_i,c(K+1,K+1))*array1-array(b_i,c(K+1,K+1))*(array2+array3-array4)*t(array(b_i,c(K+1,K+1)))
      denominator=array(at_risk_u,c(K+1,K+1))*t(array(at_risk_v,c(K+1,K+1)))
      cov_u_v=cov_u_v + numerator/denominator
    }
  }
  
  cov_u_v[cov_u_v=="NaN"]=0
  int_cov_u_v=array(NA,c(K+1,K+1))
  for(i in 1:(K+1))
  {
    for(j in 1:(K+1))
    {
      int_cov_u_v[i,j]=sum(cov_u_v[1:i,j])
    }
  }
  for(j in 2:(K+1))
  {
    int_cov_u_v[,j]=int_cov_u_v[,j-1]+int_cov_u_v[,j]
  }
  
  var_RM_array=array(exp(-lambda*a_i),c(K+1,K+1))*t(array(exp(-lambda*a_i),c(K+1,K+1)))*array(b_i,c(K+1,K+1))*int_cov_u_v*t(array(b_i,c(K+1,K+1)))
  sum(var_RM_array)
}





#EXAMPLE
#Set up the scenario
A=2 #Accrual time
s1=1 #The 1st analysis time
s2=3 #The 2nd analysis time
Tau=1 #Time length of each window
lambda=0.5

#Run the function
get_closed_cov()
#Result
#






