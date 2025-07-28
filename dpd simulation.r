# MULTINOMIAL REPRESENTATION OF INTERVAL-CENSORED UNDER GENERAL LIFETIMES
# Parameters
# ---------------------------------------
# model.distribution or distribution = {exponential, Weibull, gamma, lognormal}
# N : Number of devices
# n : sample
# tauvec : vector of times of stress change
# ITvec : vector of inspection times
# theta : parameter vector of general lifetime scale = exp(theta[1]+theta[2]x), common shape = theta[3]
# lambda: lifetime distribution parameters (different for each distribution) 
# --------------------------------------


simulate.sample <- function(theta, theta.cont, N, tauvec, stressvec, ITvec,  
                            model.distribution, seed){
  

  if (model.distribution == "lognormal"){
    lambda = c(theta[1]+ theta[2]*stressvec, theta[3]) #for logormal mui is linearly related
    lambda.cont =c(theta.cont[1]+theta.cont[2]*stressvec,theta.cont[3])
  }else{
    lambda = c(exp(theta[1]+ theta[2]*stressvec), theta[3])
    lambda.cont =c(exp(theta.cont[1]+theta.cont[2]*stressvec),theta.cont[3])
  }
 
  set.seed(seed)
  
  th = theoretical.probability(lambda, tauvec, stressvec, ITvec, model.distribution)
  th[th<1e-8]=0
  
  #contaminate j interval
  j = 3
  
  f <- paste0("model.distribution.", model.distribution)
  model.distribution = match.fun(f)
  
  th[j] = model.distribution(lambda, tauvec, stressvec, ITvec[j]) - model.distribution(lambda.cont, tauvec, stressvec, ITvec[j-1])
  th[j] = max(th[j],0)

  th = th/sum(th)
  n = drop(rmultinom(1, N, prob = th)) 
  
  return(n)
}



RMSE.function <- function(theta, theta.hat){return(mean(((theta-theta.hat)/theta)^2) )}

simulate <-function(theta, theta.cont.vec,  N, tauvec, stressvec, ITvec, model.distribution, initial = c(5, -0.02, 1), B=50){
  
  beta.list = c(0.2,0.4,0.6,0.8,1)
  nbetas = length(beta.list)+1
  res.mat = matrix(0, nrow =  length(theta.cont.vec), ncol = nbetas)
  
  for(cont in 1:length(theta.cont.vec)){
    count = 0
    
    for(b in 1:B){
      
      n = simulate.sample(theta, theta.cont.vec[[cont]], N, tauvec, stressvec, ITvec,  seed = b, model.distribution = model.distribution)
    
      estimators = estimate.function(N, tauvec, stressvec, ITvec, n, initial=initial, model.distribution = model.distribution)
      if(sum(sapply(estimators, function(x) sum(is.na(x))))<1){ #none of the estimators is nan
        res.mat[cont,-nbetas] = res.mat[cont,-nbetas]+ unlist(lapply(estimators, RMSE.function, theta = theta))
        count = count+1
      }
      res.mat[cont,nbetas] = res.mat[cont,nbetas]+ RMSE.function(theta = theta, optimum[[2]])
    }
    res.mat[cont,] = res.mat[cont,]/count
  }
  return(res.mat)
}



