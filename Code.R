#setting seed
set.seed(129)

#setting variables
K = 100
r = 0.02
sigma = 0.2
tau = 0.5
S0 = 102

#call price using Black-Scholes formula
d1 <- (log(S0/K) + (r + sigma^2/2) * T)/(sigma * sqrt(T))
d2 <- d1 - sigma * sqrt(T)
call_price <- S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
call_price

#Simple Monte Carlo Simulation

call_mc <- function(sims, tau, r, sigma, S0, K){
  
  #generating samples from the standard normal distribution
  Z <- rnorm(sims, mean=0, sd=1)
  #calculating W_T and S_T
  W_T <- Z*sqrt(tau)
  S_T <- S0*exp((r-(sigma^2)/2)*tau+sigma*W_T)
  #calculating call price for the sample
  simmed_calls <- exp(-r*tau)*pmax(S_T-K,0)
  #calculating simple monte carlo estimate and std. error
  mc_price <- mean(simmed_calls)
  mc_se <- sd(simmed_calls)/sqrt(sims)
  #returns the estimated call price and its standard error
  return(list(call_price = mc_price, standard_error = mc_se))
}

#the estimate
estimate_mc <- call_mc(sims = 1000000, tau = 0.5, r=0.02, sigma = 0.2, S0 = 102, K = 100)
estimate_mc

#Antithetic Sampling
call_as <- function(sims, tau, r, sigma, S0, K){
  
  #generating samples from the standard normal distribution
  Z <- rnorm(sims, mean=0, sd=1)
  #calculating W_T
  W_T <- Z*sqrt(tau)
  #generating antithetic samples S_T1 and S_T2
  S_T1 <- S0*exp((r-(sigma^2)/2)*tau+sigma*W_T)
  S_T2 <- S0*exp((r-(sigma^2)/2)*tau+sigma*(-W_T))
  #calculating call prices for both samples
  simmed_calls1 <- exp(-r*tau)*pmax(S_T1-K,0)
  simmed_calls2 <- exp(-r*tau)*pmax(S_T2-K,0)
  #obtaining the average of both samples
  avg_calls <- (simmed_calls1 + simmed_calls2)/2
  #calculating the estimated call price and std. error
  as_price <- mean(avg_calls)
  as_se <- sd(avg_calls)/sqrt(sims)
  #returns the estimated call price and its standard error
  return(list(call_price = as_price, standard_error = as_se))
}
set.seed(129)
estimate_as <- call_as(sims = 1000000, tau = 0.5, r=0.02, sigma = 0.2, S0 = 102, K = 100)
estimate_as

#antithetic sampling method with half of the simulations
set.seed(129)
estimate_as2 <- call_as(sims = 1000000/2, tau = 0.5, r=0.02, sigma = 0.2, S0 = 102, K = 100)
estimate_as2

call_is <- function(sims, tau, r, sigma, S0, K){
  
  #generating samples from the standard normal distribution
  Z <- rnorm(sims, mean=0, sd=1)
  #calculating W_T and S_T
  W_T <- Z*sqrt(tau)
  S_T <- S0*exp((r-(sigma^2)/2)*tau + sigma*W_T)
  #calculating call price
  simmed_calls <- (exp(-r*tau)*pmax(S_T-K,0))[S_T>K]
  #calculating the estimated call price and std. error
  is_price <- mean(simmed_calls*mean(S_T>K))
  is_se <- sd(simmed_calls*mean(S_T>K))/sqrt(sims)
  #returns the estimated call price and its standard error
  return(list(call_price = is_price, standard_error = is_se))
  
}

estimate_is <- call_is(sims = 1000000, tau = 0.5, r = 0.02, sigma = 0.2, S0 = 102, K = 100)
estimate_is
