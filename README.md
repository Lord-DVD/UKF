# Unscanted Kalman Filter and Smoother for Volatility Extraction
  
This project tries to apply a Gaussian smoothing algorithm to scaled uscented trasnformation using a Kalman filter to extract volatility from the Heston Stock Price model. The simulation and real data study shows that both stock prices and option prices are needed for accurately capturing volatility dynamics.  
  
In the state-space model framework, Bayesian optimal smoothing, also known as belief inference, refers to statistical methodology that can be used to infer state estimate using all information, which is available not only in the past and at the current, but also in the future. Optimal smoothing is closely related to optimal filtering, which makes inference based on information only available in the past and at the current time.  
  
At the core of financial econometrics is volatility estimation. Volatility pervades almost everywhere in financial markets. For example, it is used in option pricing, in portfolio allocation to control and manage risks, and in computation of risk adjusted returns for comparison of relative performance of various financial investments. Time-varying/stochastic volatility is well documented in empirical studies. There are mainly two modeling approaches for volatility. One is the class of ARCH/GARCH models (Engle, 1982; Bollerslev, 1986), where conditional volatility is a deterministic function of past volatility and return innovations, and the other is the stochastic volatility models (Shephard, 2005), which assume that volatility is unobservable and is driven by a different random process. In the past thirty years, the diffusion process has become a common tool used to model dynamics of financial data. The diffusion stochastic volatility models (Hull and White, 1987; Heston, 1993) may be the most popular ones both in academia and in practice because of their flexibility in pricing derivatives and risk management.  
  
However, in the case of stochastic volatility models, there are two main problems that persist while performing any statistical analysis.  
1. Variances of stochastic variables are state dependent.  
2. The pricing formula for derivatives is not linear.  
  
These problems combined mean that most of the conventional methods such as the standard Kalman filter and smoother are rarely applicable to these models. Even the non-linear Kalman filters, such as the extended Kalman filter and the unscented Kalman filter do not use all the information available in the past and upto current time. Another artifact to note while using these methods is that when the system is highly non-linear and high dimensional (consist of more sources of variability), the extended Kalman filter tends to perform rather poorly.  
  
Here, I have tried to use the unscented transformation approach for non-linear Gaussian system. We first use an unscented Kalman filter to approximate latent states in the model for a discretized version of the Heston stochastic volatility model, and then we smooth the data by assuming Gaussian distribution in the latent variables without making any assumption about the observation states.  
  
This project implements a simulation study using the Heston stochastic volatility model. I found that the unscented Kalman filter and smoother both perform almost identically when only stock price are supplied as observations, but according to the original paper by Junye Li (2013), the smoother improves significantly when option prices are also taken into account.  
  
I have also applied the same algorithm to the real data on S&P 500 index, and have compared errors when only applying a filter with filter and smoother.  
  
The rest of the report discusses about the scaled unscented transformation, the unscented Kalman filter as a direct application of scaled unscented transformation and its effectiveness on simulated data as well as real world data.  

You can find the link to the full report [here](https://github.com/Lord-DVD/UKF/blob/main/DesaiVatsal_ACMA830_Project1_Report.pdf)
