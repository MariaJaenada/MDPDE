# MULTINOMIAL REPRESENTATION OF INTERVAL-CENSORED DEVICES UNDER GENERAL LIFETIMES
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

source("dpd estimation.R")
source("dpd simulation.R")

# Example of use
model.distribution <- "weibull"

N = 200
stressvec = c(30,40) 
ITvec= c(6,10,14,18,20,24,28,32,36,40,44,48,52) 
tauvec = c(18,52) 


theta = c(5.3, -0.05, 1.5) 
beta.list = c(0.2,0.4,0.6,0.8,1)
#outlying cell for i = 3

n <- simulate.sample(theta = theta, theta.cont= theta, 
                    N, tauvec, stressvec, ITvec, 
                    model.distribution, seed = 1234)
                    
estimators <- estimate.function(N, tauvec, stressvec, ITvec, n, 
                initial = c(1, 0, 1), beta.list = c(0.2,0.4,0.6,0.8,1), 
                model.distribution = model.distribution)
estimators
