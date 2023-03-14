# Initialization and defining functions
setwd("D:/Academia/SFU/Fall 2022/ACMA830/Project_1")
library(readr)
S0 <- 100
V0 <- 0.04
kappa <- 3
theta <- 0.04
eta <- sigma <- 0.3
rho <- -0.6
r <- mu <- 0.05
tau <- 1
Nsteps <- 252
dt <- tau/Nsteps
obs_rho <- -0.306167243

# Function to simulate Heston model
HestonSimulation <- function(S0, V0, kappa, theta, rho, sigma, r, Nsteps, tau){
  set.seed(3)
  V <- numeric() 
  S <- numeric()
  returns <- numeric()
  dt <- tau/Nsteps 
  dw1 <- sqrt(dt) * rnorm(Nsteps) 
  dw2 <- rho*dw1+sqrt((1-rho^2)*dt)*rnorm(Nsteps) 
  S[1] <- log(S0) 
  V[1] <- V0 
  returns[1] <- 0
  for(i in 1:Nsteps){
    V[i+1] <- V[i] + kappa*(theta-V[i])*dt + sigma*sqrt(V[i]) * dw1[i]
    if( V[i+1] < 0 ){ V[i+1] <- 0 } 
    S[i+1] = S[i] + (r-0.5 * V[i]) * dt + sqrt(V[i]) * dw2[i]
    returns[i+1] <- (exp(S[i+1])-exp(S[i]))/exp(S[i])
  }
  output <- data.frame(exp(S),S,returns,V)
  colnames(output) <- c("StockPrice","LogPrice","Returns","Volatility")
  return(output)
}

# function to generate stock prices from volatility values
StockPredictor <- function(UKF_Vol,S0,V0,r,rho,tau,Nsteps){
  set.seed(3)
  V <- c(V0,V0,UKF_Vol$X1)
  S <- numeric()
  returns <- numeric()
  dt <- tau/Nsteps 
  dw1 <- sqrt(dt) * rnorm(Nsteps) 
  dw2 <- rho*dw1+sqrt((1-rho^2)*dt)*rnorm(Nsteps) 
  S[1] <- log(S0)
  returns[1] <- 0
  for(i in 1:Nsteps){
    if(V[i+1] < 0){V[i+1] <- 0} 
    S[i+1] = S[i] + (r-0.5 * V[i]) * dt + sqrt(V[i]) * dw2[i]
    returns[i+1] <- (exp(S[i+1])-exp(S[i]))/exp(S[i])
  }
  output <- data.frame(exp(S),S,returns,V)
  colnames(output) <- c("StockPrice","LogPrice","Returns","UKF_Volatility")
  return(output)
}

# Predictions using filtered and smoother vol for simulated data

path1 <- HestonSimulation(S0,V0,kappa,theta,rho,sigma,r,Nsteps,tau)
write.csv(path1,"SimulatedHestonModel.csv")

system("matlab -nodisplay -r \"run('UKS_sim.m'); exit\"")

UKF_Vol_1 <- read_csv("UKF_Vol_1.csv", col_names = FALSE, show_col_types = F)
UKF_Vol_2 <- read_csv("UKF_Vol_2.csv", col_names = FALSE, show_col_types = F)

pred1 <- StockPredictor(UKF_Vol_1,S0,V0,r,rho,tau,Nsteps)
pred2 <- StockPredictor(UKF_Vol_2,S0,V0,r,rho,tau,Nsteps)

# Graphs for comparison

par(mfrow=c(3,2))
plot(path1$Returns,type="l", main="Simulated Stock Returns", xlab = "Days", ylab="Returns", col="blue")
plot(path1$Volatility,type = "l", main="Simulated Volatility", xlab = "Days", ylab="Volatility", col="blue")
plot(pred1$Returns,type="l", main="Filtered Stock Returns", xlab = "Days", ylab="Returns", col="red")
lines(path1$Returns, col="blue")
plot(pred1$UKF_Volatility,type="l", main="Filtered Volatility", xlab = "Days", ylab="Volatility", col="red")
lines(path1$Volatility, col="blue")
plot(pred2$Returns,type="l", main="Smoothed Stock Returns", xlab = "Days", ylab="Returns", col="red")
lines(path1$Returns, col="blue")
plot(pred2$UKF_Volatility,type="l", main="Smoothed Volatility", xlab = "Days", ylab="Returns", col="red")
lines(path1$Volatility, col="blue")

# Predictions using filtered and smoothed vol on real data

SP500 <- read_csv("S&P500.csv",col_types=cols(Date=col_date(format="%Y-%m-%d")))

system("matlab -nodisplay -r \"run('UKS_SP500.m'); exit\"")

UKF_Vol_3 <- read_csv("UKF_Vol_3.csv", col_names = FALSE, show_col_types = F)
UKF_Vol_4 <- read_csv("UKF_Vol_4.csv", col_names = FALSE, show_col_types = F)

pred3 <- StockPredictor(UKF_Vol_3,SP500$StockPrice[1],SP500$Volatility[1],0.053,obs_rho,1,252)
pred4 <- StockPredictor(UKF_Vol_4,SP500$StockPrice[1],SP500$Volatility[1],0.053,obs_rho,1,252)

# Graphs for Comparison

par(mfrow=c(3,2))
plot(SP500$Returns,type="l", main="Realized Stock Returns", xlab = "Days", ylab="Returns", col="blue")
plot(SP500$Volatility,type = "l", main="B.S. Implied Volatility", xlab = "Days", ylab="Volatility", col="blue")
plot(pred3$Returns,type="l", main="Filtered Stock Returns", xlab = "Days", ylab="Returns", col="red")
lines(SP500$Returns, col="blue")
plot(pred3$UKF_Volatility,type="l", main="Filtered Volatility", xlab = "Days", ylab="Volatility", col="red")
lines(SP500$Volatility, col="blue")
plot(pred4$Returns,type="l", main="Smoothed Stock Returns", xlab = "Days", ylab="Returns", col="red")
lines(SP500$Returns, col="blue")
plot(pred4$UKF_Volatility,type="l", main="Smoothed Volatility", xlab = "Days", ylab="Returns", col="red")
lines(SP500$Volatility, col="blue")

# Trying different values of parameters on simulated data

path2 <- HestonSimulation(S0,V0,0.5,theta,rho,sigma,r,Nsteps,tau)
write.csv(path2,"SimulatedHestonModel.csv")
system("matlab -nodisplay -r \"run('UKS_sim_parameters.m'); exit\"")

UKF_Vol_5 <- read_csv("UKF_Vol_5.csv", col_names = FALSE, show_col_types = F)
UKF_Vol_6 <- read_csv("UKF_Vol_6.csv", col_names = FALSE, show_col_types = F)
pred5 <- StockPredictor(UKF_Vol_5,S0,V0,mu,rho,tau,Nsteps)
pred6 <- StockPredictor(UKF_Vol_6,S0,V0,mu,rho,tau,Nsteps)



path3 <- HestonSimulation(S0,V0,10,theta,rho,sigma,r,Nsteps,tau)
system("matlab -nodisplay -r \"run('UKS_sim_parameters.m'); exit\"")

UKF_Vol_7 <- read_csv("UKF_Vol_5.csv", col_names = FALSE, show_col_types = F)
UKF_Vol_8 <- read_csv("UKF_Vol_6.csv", col_names = FALSE, show_col_types = F)
pred7 <- StockPredictor(UKF_Vol_7,S0,V0,mu,rho,tau,Nsteps)
pred8 <- StockPredictor(UKF_Vol_8,S0,V0,mu,rho,tau,Nsteps)

# Comparison graphs for different parameters

par(mfrow=c(3,2))
plot(path2$Returns,type="l", main="Simulated Returns (k=0.5)", xlab = "Days", ylab="Returns", col="blue")
plot(path3$Returns,type="l", main="Simulated Returns (k=10)", xlab = "Days", ylab="Returns", col="blue")
plot(path2$Volatility,type = "l", main="Filtered Volatility", xlab = "Days", ylab="Volatility", col="blue")
lines(pred5$UKF_Volatility, col="red")
plot(path3$Volatility,type = "l", main="Filtered Volatility", xlab = "Days", ylab="Volatility", col="blue")
lines(pred7$UKF_Volatility, col="red")
plot(path2$Volatility,type = "l", main="Smoothed Volatility", xlab = "Days", ylab="Volatility", col="blue")
lines(pred6$UKF_Volatility, col="red")
plot(path3$Volatility,type = "l", main="Smoothed Volatility", xlab = "Days", ylab="Volatility", col="blue")
lines(pred8$UKF_Volatility, col="red")

# Real Data Error Graphs

par(mfrow=c(3,2))
plot(SP500$Returns,type="l", main="Realized Stock Returns", xlab = "Days", ylab="Returns", col="blue")
plot(SP500$Volatility,type = "l", main="B.S. Implied Volatility", xlab = "Days", ylab="Volatility", col="blue")
plot(pred3$UKF_Volatility,type="l", main="Filtered Volatility", xlab = "Days", ylab="Volatility", col="red")
lines(SP500$Volatility, col="blue")
plot(pred4$UKF_Volatility,type="l", main="Smoothed Volatility", xlab = "Days", ylab="Returns", col="red")
lines(SP500$Volatility, col="blue")
plot(pred3$UKF_Volatility-SP500$Volatility, type="l", main = "Difference: Filtered", xlab = "Days", ylab = "Error", col="blue")
plot(pred4$UKF_Volatility-SP500$Volatility, type="l", main = "Difference: Smoothed", xlab = "Days", ylab = "Error", col="blue")
