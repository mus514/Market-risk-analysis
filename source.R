f_interpolate <- function(t1, r1, t2, r2, t) 
  {
  ### Calcule l'interpolation lineaire du taux d'interet entre t1 et t2
  #  INPUTS
  #   t1    : [scalar] temps minimum
  #   t2    : [scalar] temps maximum
  #   r1    : [scalar] interet au t1
  #   r2    : [scalar] interet au t2
  #   t    : [scalar]  temps target
  #  OUTPUTS
  #   prix   : [scalar] taux d'interet 
  
  r <- r1 + (r2 - r1)*(t - t1)/(t2 - t1)
  r
}


f_BS <- function(s, k, sigma, r, t, type = "Call")
{
  ### Calcule le prix d'option Call et Put avec la formule Black-Schules
  #  INPUTS
  #   S    : [scalar] prix du S-J
  #   K    : [scalar] Strike pirce
  #   sigma: [scalar] volatility
  #   r    : [scalar] taux interet
  #   t    : [scalar] maturity
  #   type : [string] type de l'option "Call" ou "Put"
  #  OUTPUTS
  #   p   : [scalar] prix de l'option

  d1 <- (log(s/k) + (r + sigma^2/2)*t)/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)

  if(type == "Call")
    prix <- s*pnorm(d1, mean = 0, sd = 1) - k*exp(-r*t)*pnorm(d2, mean = 0, sd = 1)

  else
    prix <- -s*pnorm(-d1, mean = 0, sd = 0) + k*exp(-r*t)*pnorm(-d2, mean = 0, sd = 0)

  prix
}




f_sim <- function(s, n_sim, n_day, mu, sigma, t)
{
  ### Fonction qui simule le mouvement Brownien 
  #  INPUTS
  #   s     : [scalar] prix du S-J au t = 0
  #   mu    : [scalar] rendement moyen
  #   sigma : [scalar] volatility
  #   n_day : [scalar] nombre de jours de simulation
  #   n_sim : [scalar] nombre de simulation
  #   t     : [scalar] interval de simulation
  #  OUTPUTS
  #   output   : [Matrix] prix simulés
  
  output <- matrix(NA, nrow = n_sim, ncol = n_day+1)
  output[,1] <- rep(s, n_sim)
  
  for(i in 1:n_sim)
  {
    for(j in 2:(n_day+1))
    {
      output[i,j] <- output[i, j-1]*exp((mu - 0.5*sigma^2)*t + sigma*sqrt(t)*rnorm(1))
    }
  }
  
  output
}




f_simBV <- function(s, vix, n_sim, n_day, mu, sigma, phi, t)
{
  ### Fonction qui simule le mouvement Brownien bivarié 
  #  INPUTS
  #   s     : [scalar] prix du S-J au t = 0
  #   vix   : [scalar] indice du vix t = 0
  #   mu    : [vector] rendement moyen du s-j et du vix
  #   sigma : [vector] volatilite du s-j et vix
  #   phi   : [scalar] correlation entre s-j et vix 
  #   n_day : [scalar] nombre de jours de simulation
  #   n_sim : [scalar] nombre de simulation
  #   t     : [scalar] interval de simulation
  #  OUTPUTS
  #   output: [Matrix] vix et s-j au jours n_day
  
  output1 <- matrix(NA, nrow = n_sim, ncol = n_day+1)
  output2 <- matrix(NA, nrow = n_sim, ncol = n_day+1)
  
  output1[,1] <- rep(s, n_sim)
  output2[,1] <- rep(vix, n_sim)
  
  m <- matrix(c(1, phi, phi, 1), nrow = 2)
  
  for(i in 1:n_sim)
  {
    for(j in 2:(n_day+1))
    {
      x <- mvtnorm::rmvnorm(n = 1, mean = c(0, 0), sigma = m)
      output1[i,j] <- output1[i, j-1]*exp((mu[1] - 0.5*sigma[1]^2)*t + sigma[1]*sqrt(t)*x[1])
      output2[i,j] <- output2[i, j-1]*exp((mu[2] - 0.5*sigma[2]^2)*t + sigma[2]*sqrt(t)*x[2])
    }
  }
  
  output <- list(output1, output2)
  output
}



f_simCP <- function(s, vix, n_day, n_sim, theta1, theta2, fit, t)
{
  ### Fonction qui simule le mouvement Brownien bivarié 
  #  INPUTS
  #   s     : [scalar] prix du S-J au t = 0
  #   vix   : [scalar] indice du vix t = 0
  #   theta1: [vector] moyenne et volatilite du s-j
  #   theta2: [vector] moyenne et volatilite du vix
  #   fit   : [function] fonction de copule utilisee 
  #   n_day : [scalar] nombre de jours de simulation
  #   n_sim : [scalar] nombre de simulation
  #   t     : [scalar] interval de simulation
  #  OUTPUTS
  #   output: [Matrix] simulation de vix et s-j pour les n_day
  
  output1 <- matrix(NA, nrow = n_sim, ncol = n_day+1)
  output2 <- matrix(NA, nrow = n_sim, ncol = n_day+1)
  
  output1[,1] <- rep(s, n_sim)
  output2[,1] <- rep(vix, n_sim)
  
  for(i in 1:n_sim)
  {
    for(j in 2:(n_day+1))
    {
      x <- rCopula(1, fit@copula)
      q1 <- qt(x[1], 10)
      q2 <- qt(x[2], 5)
      
      output1[i,j] <- output1[i, j-1]*exp((theta1[1] - 0.5*theta1[2]^2)*t + theta1[2]*sqrt(t)*q1)
      output2[i,j] <- output2[i, j-1]*exp((theta2[1] - 0.5*theta2[2]^2)*t + theta2[2]*sqrt(t)*q2)
    }
  }
  
  output <- list(output1, output2)
  output
}


