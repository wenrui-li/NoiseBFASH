library(nimble)
library(Matrix)
library(statmod)
library(pgdraw)
library(MASS)
source("MCMC.R")
 
# load data
load("../Data/data.RData")

# set parameters 
H <- 3
L <- 9
nu1 <- 1
nu2 <- 1
a_rho <- b_rho <- a_phi <- b_phi <- 1
a_omega <- 2
b_omega <- 1  
sigma2_m <- 1
n_iter <-  1e4
burnin <- 1e3
k <- 200
a_q <- b_q <- rep(1,H)
uj1 <- 0.5
u00_mean_prior <- rep(0,H)
u10_mean_prior <- rep(0,H)
u00_star_mean_prior <- rep(0,H)
u10_star_mean_prior <- rep(0,H) 
sigma_u00 <- sigma_u10 <- sigma_u00_star <- sigma_u10_star <- rep(1,H) 
r <- rep(5,H)
r_star <- rep(5,H)
r_total <- r + r_star 
u1_mean_prior <- u0_mean_prior <- lapply(1:H,function(h) matrix(0,r_total[h]-2,1))
Sigma0 <- Sigma1 <- lapply(1:H,function(h) diag(r_total[h]-2))


# run MCMC
fit <- MCMC(data$X,data$A_tilde_list,A_star_list,L,nu1,nu2,a_rho,b_rho,a_omega,b_omega,sigma2_m,a_phi,b_phi,
                        a_q,b_q,u00_mean_prior,sigma_u00,u10_mean_prior,sigma_u10,
                        u00_star_mean_prior,sigma_u00_star,u10_star_mean_prior,sigma_u10_star,
                        r,r_star,uj1,u0_mean_prior,u1_mean_prior,Sigma0,Sigma1, n_iter,k,burnin) 

# estimated factor loading matrix W
W <- apply(fit$W[,,(burnin+1):n_iter],c(1,2),mean)
# estimated factor matrix Z
Z <- apply(fit$Z[,,(burnin+1):n_iter],c(1,2),mean)
# DIC
DIC <- fit$DIC
