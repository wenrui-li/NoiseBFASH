# MCMC algorihtm

# X: features, p by n matrix
# A_tilde: list of the existing knowledge graph, list of p_h by p_h matrix, h=1,...,H
# A_star: list of the graph estimated by graphical lasso, list of p_h by p_h matrix, h=1,...,H
# L: the number of factors, integer
# nu1, nu2 : parameters of the prior distribution of theta, numeric
# a_rho, b_rho: parameters of the prior distribution of $rho_j^{(h)}$, numeric
# a_omega, b_omega: parameters of the prior distribution of Omega, numeric
# sigma2_m: parameters of the prior distribution of $m_j^{(h)}$, numeric
# a_phi,b_phi: parameters of the prior distribution of $\phi_l^{(h)}$, numeric
# a_q, b_q: parameters of the prior distribution of the true graph, numeric
# u00_mean_prior, sigma_u00: parameters of the prior distribution of $\tilde s_0^{(h)}$, h=1,...,H, list of numeric values
# u10_mean_prior, sigma_u10: parameters of the prior distribution of $\tilde s_1^{(h)}$, h=1,...,H, list of numeric values
# u00_star_mean_prior, sigma_u00_star: parameters of the prior distribution of $s*_0^{(h)}$, h=1,...,H, list of numeric values
# u10_star_mean_prior, sigma_u10_star: parameters of the prior distribution of $s*_1^{(h)}$, h=1,...,H, list of numeric values
# r: length of latent scales $\tilde u_{0j}^{(h)}$ and $\tilde u_{1j}^{(h)}$, h=1,...,H, list of integers
# r_star: length of latent scales $u*_{0j}^{(h)}$ and $u*_{1j}^{(h)}$, h=1,...,H, list of integers
# uj1: the first element of every latent scale, numeric
# u0_mean_prior: mean of the prior distribution of (\tilde u_{0j}^{(h)}[2:r[h]],u*_{0j}^{(h)}[2:r_star[h]]), h=1,...,H, list of vector of length r[h]+r_star[h]-2
# u1_mean_prior: mean of the prior distribution of (\tilde u_{1j}^{(h)}[2:r[h]],u*_{1j}^{(h)}[2:r_star[h]]), h=1,...,H, list of vector of length r[h]+r_star[h]-2
# Sigma0: variance of the prior distribution of (\tilde u_{0j}^{(h)}[2:r[h]],u*_{0j}^{(h)}[2:r_star[h]]),  h=1,...,H, list of r[h]+r_star[h]-2 by r[h]+r_star[h]-2 matrix
# Sigma1: variance of the prior distribution of (\tilde u_{1j}^{(h)}[2:r[h]],u*_{1j}^{(h)}[2:r_star[h]]),  h=1,...,H, list of r[h]+r_star[h]-2 by r[h]+r_star[h]-2 matrix
# n_iter: number of iterations, integer
# k: interval length for computing the acceptance rate of Metropolis/MH algorithm, integer. (set k < n_iter if we want to update proposal variance)
# burnin: number of burn-in, integer

MCMC <- function(X,A_tilde,A_star,L,nu1,nu2,a_rho,b_rho,a_omega,b_omega,sigma2_m,a_phi,b_phi,
                 a_q,b_q,u00_mean_prior,sigma_u00,u10_mean_prior,sigma_u10,
                 u00_star_mean_prior,sigma_u00_star,u10_star_mean_prior,sigma_u10_star,
                 r,r_star,uj1,u0_mean_prior,u1_mean_prior,Sigma0,Sigma1, n_iter,k,burnin) 
{
  n = ncol(X)
  p = nrow(X)
  H = length(A_tilde)
  # numbers of features p_h, p_vec[h]
  p_vec = sqrt(lengths(A_tilde))
  # index for modality h, p_index[[h]]
  p_group = cut(1:p,c(0,cumsum(p_vec)), labels = 1:H)
  p_index = split(1:p, p_group)
  # get h, p_mod[j]
  p_mod = as.numeric(p_group)
  # get j for each modality
  j_vec = unlist(sapply(1:H, function(h) 1:p_vec[h]))
  
  W = array(0,c(p,L,n_iter))
  M = matrix(0, p,n_iter)
  RHO = matrix(1, p,n_iter)
  Z = array(0,c(L,n,n_iter))
  PHI = array(1,c(H,L,n_iter))
  THETA = array(0,c(p,L,n_iter))
  TAU = array(1,c(p,L,n_iter))
  psi = matrix(0, p,n)
  loglikelihood_x= NULL
  
  m = M[,1]
  w = W[,,1]
  z = Z[,,1]
  rho = RHO[,1]
  phi = PHI[,,1]
  theta = THETA[,,1]
  tau = TAU[,,1]
  theta_sd_prop = matrix(1,p,L)
  acs_theta = array(0,c(p,L,n_iter))
  r_theta = array(0,c(p,L,floor(n_iter/k)))
  r_opt = 0.44
  c1 = 0.8
  
  allone_col = matrix(1,n,1)
  allone_row = matrix(1,1,n)
  allzero_L= rep(0,L)
  
  rho_shape_post = a_rho+n/2
  phi_shape_post = a_phi + 3*p_vec/2
  
  a = A_tilde
  a_ind = list()
  for(h in 1:H){
    a_ind = append(a_ind,lapply(1:p_vec[h], function(pp) which(a[[h]][,pp]==1))) 
  }
  
  gamma_a_omega = gamma(a_omega)
  Q = matrix(a_q/(a_q+b_q), H,n_iter)
  q = Q[,1]
  Kappa_tilde = lapply(1:H, function(h) A_tilde[[h]]-0.5+diag(0.5,p_vec[h],p_vec[h]))
  Kappa_star = lapply(1:H, function(h) A_star[[h]]-0.5+diag(0.5,p_vec[h],p_vec[h]))
  
  U00 = matrix(u00_mean_prior,H,n_iter)
  U10 = matrix(u10_mean_prior,H,n_iter)
  u00 = U00[,1]
  u10 = U10[,1]
  U00_star = matrix(u00_star_mean_prior,H,n_iter)
  U10_star = matrix(u10_star_mean_prior,H,n_iter)
  u00_star = U00_star[,1]
  u10_star = U10_star[,1]
  
  r_total = r + r_star 
  index_r_tilde = lapply(1:H, function(h) 1:(r[h]-1))
  index_r_star = lapply(1:H, function(h) r[h]:(r_total[h]-2))
  
  u0 = lapply(1:H, function(h)  matrix(c(uj1,u0_mean_prior[[h]][index_r_tilde[[h]]]),p_vec[h],r[h],byrow = TRUE ) )
  u1 = lapply(1:H, function(h)  matrix(c(uj1,u1_mean_prior[[h]][index_r_tilde[[h]]]),p_vec[h],r[h],byrow = TRUE ) )
  u0_star = lapply(1:H, function(h)  matrix(c(uj1,u0_mean_prior[[h]][index_r_star[[h]]]),p_vec[h],r_star[h],byrow = TRUE ) )
  u1_star = lapply(1:H, function(h)  matrix(c(uj1,u1_mean_prior[[h]][index_r_star[[h]]]),p_vec[h],r_star[h],byrow = TRUE ) )
  
  Xi0_tilde =  Xi1_tilde = Xi0_star =  Xi1_star = lapply(1:H, function(h)  matrix(1,p_vec[h],p_vec[h]))
  
  upper_index = lapply(1:H, function(h) upper.tri(a[[h]]))
  a_sum = rep(0,H)
  
  tau_inv_u00_tilde = .5/sigma_u00^2
  tau_inv_u10_tilde = .5/sigma_u10^2
  tau_inv_u00_star = .5/sigma_u00_star^2
  tau_inv_u10_star = .5/sigma_u10_star^2
  
  Sigma0_11 = lapply(1:H, function(h) Sigma0[[h]][index_r_tilde[[h]],index_r_tilde[[h]]])
  Sigma0_12 = lapply(1:H, function(h) Sigma0[[h]][index_r_tilde[[h]],index_r_star[[h]]])
  Sigma0_21 = lapply(1:H, function(h) Sigma0[[h]][index_r_star[[h]],index_r_tilde[[h]]])
  Sigma0_22 = lapply(1:H, function(h) Sigma0[[h]][index_r_star[[h]],index_r_star[[h]]])
  Sigma1_11 = lapply(1:H, function(h) Sigma1[[h]][index_r_tilde[[h]],index_r_tilde[[h]]])
  Sigma1_12 = lapply(1:H, function(h) Sigma1[[h]][index_r_tilde[[h]],index_r_star[[h]]])
  Sigma1_21 = lapply(1:H, function(h) Sigma1[[h]][index_r_star[[h]],index_r_tilde[[h]]])
  Sigma1_22 = lapply(1:H, function(h) Sigma1[[h]][index_r_star[[h]],index_r_star[[h]]])
  
  u0_mean_prior_1 = lapply(1:H, function(h) u0_mean_prior[[h]][index_r_tilde[[h]]])
  u1_mean_prior_1 = lapply(1:H, function(h) u1_mean_prior[[h]][index_r_tilde[[h]]])
  u0_mean_prior_2 = lapply(1:H, function(h) u0_mean_prior[[h]][index_r_star[[h]]])
  u1_mean_prior_2 = lapply(1:H, function(h) u1_mean_prior[[h]][index_r_star[[h]]])  
  
  
  Sigma0_part1 = lapply(1:H, function(h) Sigma0_12[[h]]%*%solve(Sigma0_22[[h]]))
  Sigma0_part2 = lapply(1:H, function(h) Sigma0_21[[h]]%*%solve(Sigma0_11[[h]]))
  Sigma1_part1 = lapply(1:H, function(h) Sigma1_12[[h]]%*%solve(Sigma1_22[[h]]))
  Sigma1_part2 = lapply(1:H, function(h) Sigma1_21[[h]]%*%solve(Sigma1_11[[h]]))
  
  Sigma0_tilde =  lapply(1:H, function(h) Sigma0_11[[h]]-Sigma0_part1[[h]]%*%Sigma0_21[[h]])
  Sigma0_tilde_inv = lapply(1:H, function(h) solve(Sigma0_tilde[[h]]))
  Sigma1_tilde =  lapply(1:H, function(h) Sigma1_11[[h]]-Sigma1_part1[[h]]%*%Sigma1_21[[h]])
  Sigma1_tilde_inv = lapply(1:H, function(h) solve(Sigma1_tilde[[h]]))
  
  Sigma0_star = lapply(1:H, function(h) Sigma0_22[[h]]-Sigma0_part2[[h]]%*%Sigma0_12[[h]])
  Sigma0_star_inv = lapply(1:H, function(h) solve(Sigma0_star[[h]]))
  Sigma1_star = lapply(1:H, function(h) Sigma1_22[[h]]-Sigma1_part2[[h]]%*%Sigma1_12[[h]])
  Sigma1_star_inv = lapply(1:H, function(h) solve(Sigma1_star[[h]]))
  
  for (i in 1:n_iter) {
    # update m
    m_var_post = 1/(1/sigma2_m+n*rho)
    m_mean_post = m_var_post*rho*rowSums(X-w%*%z)
    m = rnorm(p,m_mean_post,sqrt(m_var_post))
    M[,i] = m 
    X_m = X-m
    
    # update rho
    rho_rate_post=rowSums((X_m-w%*%z)^2)/2+b_rho
    rho = rgamma(p,shape = rho_shape_post,rate = rho_rate_post)
    RHO[,i] = rho
    
    # update theta
    for(j in 1:p) {
      h_j = p_mod[j]
      theta_j = theta[j,]
      index_k = a_ind[[j]]
      lengh_k = length(index_k)
      
      if(lengh_k>0) {
        theta_k = matrix(theta[index_k,],lengh_k,L)
        sum_jk = colSums((t(theta_k)-theta_j)^2)/2/nu2+b_omega
        for(l in 1:L) {
          theta_j_prop = theta_j
          theta_jl = theta_j[l]
          theta_jl_prop = rnorm(1,theta_jl,theta_sd_prop[j,l])
          theta_j_prop[l] = theta_jl_prop
          sum_jk_diff = (theta_jl^2-theta_jl_prop^2)/2/nu2+(theta_jl_prop-theta_jl)/nu2*theta_k[,l]
          sum_jk_prop = sum_jk - sum_jk_diff
          lhr = 2*theta_jl_prop - 2*theta_jl + 
            (exp(2*theta_jl)-exp(2*theta_jl_prop))*tau[j,l]*phi[h_j,l]/2 +
            - (theta_jl_prop-nu1)^2/2/nu2 + (theta_jl-nu1)^2/2/nu2+
            -a_omega*sum(log(sum_jk_prop))+a_omega*sum(log(sum_jk))
          if( log(runif(1)) <lhr ) {theta_j[l] = theta_jl_prop; acs_theta[j,l,i] = acs_theta[j,l,i] + 1;sum_jk = sum_jk_prop}
        }
      } else {
        for(l in 1:L) {
          theta_j_prop = theta_j
          theta_jl = theta_j[l]
          theta_jl_prop = rnorm(1,theta_jl,theta_sd_prop[j,l])
          theta_j_prop[l] = theta_jl_prop 
          lhr = 2*theta_jl_prop - 2*theta_jl + 
            (exp(2*theta_jl)-exp(2*theta_jl_prop))*tau[j,l]*phi[h_j,l]/2 +
            - (theta_jl_prop-nu1)^2/2/nu2 + (theta_jl-nu1)^2/2/nu2
          if( log(runif(1)) <lhr ) {theta_j[l] = theta_jl_prop; acs_theta[j,l,i] = acs_theta[j,l,i] + 1}
        }
      }
      
      theta[j,] = theta_j
    } 
    THETA[,,i] = theta
    
    # update z
    B_z_inv = crossprod(w,w*rho) 
    diag(B_z_inv) = diag(B_z_inv) + 1 
    B_z_inv_chol = chol(B_z_inv)
    B_z = chol2inv(B_z_inv_chol)
    z_mean_post = tcrossprod(B_z,w)%*%(X_m*rho)
    z = z_mean_post+replicate(n,rmnorm_chol(1,allzero_L,B_z_inv_chol, prec_param = TRUE) )
    Z[,,i] = z
    
    # update phi
    phi_rate_tmp = tau*exp(2*theta)+w^2/tau
    for(h in 1:H){
      phi[h,] = rgamma(L,shape = phi_shape_post[h], rate = b_phi+colSums(phi_rate_tmp[p_index[[h]],])/2) 
    }
    PHI[,,i] = phi
    
    # update w_j & tau_j
    zz = tcrossprod(z,z) 
    for(j in 1:p){
      phi_j = phi[p_mod[j],]
      B_j_inv = rho[j]*zz
      diag(B_j_inv) = diag(B_j_inv) + phi_j/tau[j,]
      B_j_inv_chol = chol(B_j_inv)
      B_j = chol2inv(B_j_inv_chol)
      w_j_mean_post = B_j%*%z%*%X_m[j,]*rho[j] 
      w_j = w_j_mean_post+rmnorm_chol(1,allzero_L,B_j_inv_chol, prec_param = TRUE)  
      w[j,] = w_j 
      lambda_j = exp(theta[j,])
      tau[j,] = 1/rinvgauss(L,lambda_j/abs(w_j),lambda_j^2*phi_j)
    }
    W[,,i]=w
    TAU[,,i]=tau
    
    for(h in 1:H){
      # update a
      p_h = p_vec[h]
      a_tmp = Xi0_tilde_tmp=Xi1_tilde_tmp=Xi0_star_tmp=Xi1_star_tmp=matrix(0, p_h, p_h)
      upper_index_h = upper_index[[h]]
      uu1_part = tcrossprod(u1[[h]],u1[[h]])
      uu1 = u10[h]+ uu1_part
      uu0_part = tcrossprod(u0[[h]],u0[[h]])
      uu0 = u00[h]+ uu0_part 
      uu1_star_part = tcrossprod(u1_star[[h]],u1_star[[h]])
      uu1_star = u10_star[h]+uu1_star_part
      uu0_star_part =  tcrossprod(u0_star[[h]],u0_star[[h]])
      uu0_star = u00_star[h]+uu0_star_part
      theta_diff = matrix(rowSums(sapply(1:L, function(l) (outer(theta[p_index[[h]],l],theta[p_index[[h]],l],"-"))^2) ),p_h,p_h)
      prob_tmp = q[h]/(1-q[h])*gamma_a_omega/(b_omega+(theta_diff)/nu2/2)^(a_omega)*exp(Kappa_tilde[[h]]*(uu1-uu0)-Xi1_tilde[[h]]*uu1*uu1/2+Xi0_tilde[[h]]*uu0*uu0/2+
                                                                                          Kappa_star[[h]]*(uu1_star-uu0_star)-Xi1_star[[h]]*uu1_star*uu1_star/2+Xi0_star[[h]]*uu0_star*uu0_star/2 )
      prob_tmp = prob_tmp/(1+prob_tmp)
      a_tmp[upper_index_h] = rbinom(p_h*(p_h-1)/2,1,prob_tmp[upper_index_h])
      
      a_tmp_op = 1-a_tmp
      a[[h]] = a_tmp +t(a_tmp)
      a_sum[h] =sum(a[[h]])/2
      
      # update Xi0_tilde, Xi1_tilde, Xi0_star,Xi1_star
      Xi0_tilde_para = uu0-a_tmp*uu0
      Xi0_tilde_tmp[upper_index_h] = pgdraw(1, Xi0_tilde_para[upper_index_h])
      Xi1_tilde_para = a_tmp*uu1
      Xi1_tilde_tmp[upper_index_h] = pgdraw(1, Xi1_tilde_para[upper_index_h])
      Xi0_star_para = uu0_star-a_tmp*uu0_star
      Xi0_star_tmp[upper_index_h] = pgdraw(1, Xi0_star_para[upper_index_h])
      Xi1_star_para = a_tmp*uu1_star
      Xi1_star_tmp[upper_index_h] = pgdraw(1, Xi1_star_para[upper_index_h])
      
      Xi0_tilde[[h]] = Xi0_tilde_tmp+t(Xi0_tilde_tmp)
      Xi1_tilde[[h]] = Xi1_tilde_tmp+t(Xi1_tilde_tmp)
      Xi0_star[[h]] = Xi0_star_tmp+t(Xi0_star_tmp)
      Xi1_star[[h]] = Xi1_star_tmp+t(Xi1_star_tmp)
      
      # update u00, u10, u00_star, u10_star
      Xi0_tilde_sum0 = sum(Xi0_tilde_tmp*a_tmp_op) 
      Xi1_tilde_sum1 = sum(Xi1_tilde_tmp*a_tmp)
      Xi0_star_sum0 = sum(Xi0_star_tmp*a_tmp_op) 
      Xi1_star_sum1 = sum(Xi1_star_tmp*a_tmp)
      
      Kappa_tilde_sum1 = sum(Kappa_tilde[[h]]*a_tmp) 
      Kappa_tilde_sum0 = sum(Kappa_tilde[[h]])/2- Kappa_tilde_sum1 
      Kappa_star_sum1 = sum(Kappa_star[[h]]*a_tmp) 
      Kappa_star_sum0 = sum(Kappa_star[[h]])/2- Kappa_star_sum1 
      
      u00[h] = rnorm(1,(u00_mean_prior[h]*tau_inv_u00_tilde[h]+Kappa_tilde_sum0-sum(a_tmp_op*Xi0_tilde_tmp*uu0_part))/(tau_inv_u00_tilde[h]+Xi0_tilde_sum0),1/sqrt(tau_inv_u00_tilde[h]+Xi0_tilde_sum0))
      u10[h] = rnorm(1,(u10_mean_prior[h]*tau_inv_u10_tilde[h]+Kappa_tilde_sum1-sum(a_tmp*Xi1_tilde_tmp*uu1_part))/(tau_inv_u10_tilde[h]+Xi1_tilde_sum1),1/sqrt(tau_inv_u10_tilde[h]+Xi1_tilde_sum1))
      u00_star[h] = rnorm(1,(u00_star_mean_prior[h]*tau_inv_u00_star[h]+Kappa_star_sum0-sum(a_tmp_op*Xi0_star_tmp*uu0_star_part))/(tau_inv_u00_star[h]+Xi0_star_sum0),1/sqrt(tau_inv_u00_star[h]+Xi0_star_sum0))
      u10_star[h] = rnorm(1,(u10_star_mean_prior[h]*tau_inv_u10_star[h]+Kappa_star_sum1-sum(a_tmp*Xi1_star_tmp*uu1_star_part))/(tau_inv_u10_star[h]+Xi1_star_sum1),1/sqrt(tau_inv_u10_star[h]+Xi1_star_sum1))
    }
    
    # non-zero & zero in a
    a_ind =  a_ind0 =list()
    for(h in 1:H){
      a_ind = append(a_ind,lapply(1:p_vec[h], function(pp) which(a[[h]][,pp]==1))) 
      a_ind0 = append(a_ind0,lapply(1:p_vec[h], function(pp) setdiff(which(a[[h]][,pp]==0),pp))) 
    }
    
    # update q
    a_q_pos = a_q+ a_sum
    b_q_pos = b_q+ p_vec*(p_vec-1)/2- a_sum
    q = rbeta(H,a_q_pos,b_q_pos)
    
    # update tau_inv_u00_tilde, tau_inv_u10_tilde, tau_inv_u00_star, tau_inv_u10_star
    tau_inv_u00_tilde = rinvgauss(1/sigma_u00/abs(u00-u00_mean_prior),1/sigma_u00^2)
    tau_inv_u10_tilde = rinvgauss(1/sigma_u10/abs(u10-u10_mean_prior),1/sigma_u10^2)
    tau_inv_u00_star = rinvgauss(1/sigma_u00_star/abs(u00_star-u00_star_mean_prior),1/sigma_u00_star^2)
    tau_inv_u10_star = rinvgauss(1/sigma_u10_star/abs(u10_star-u10_star_mean_prior),1/sigma_u10_star^2)
    
    # update u0 
    for(j in 1:p) {
      h_j = p_mod[j]
      j_j = j_vec[j]
      index_k = a_ind0[[j]]
      lengh_k = length(index_k)
      
      if(lengh_k>0) {
        U0j_tilde = matrix(u0[[h_j]][index_k,-1],nrow=lengh_k)
        xi0j_tilde = Xi0_tilde[[h_j]][index_k,j_j]
        kappa0j_tilde = Kappa_tilde[[h_j]][index_k,j_j]
        U0j_xi0j_tilde=U0j_tilde*xi0j_tilde
        O0j_tilde = kappa0j_tilde/xi0j_tilde+u00[h_j]
        V0j_tilde = crossprod(U0j_xi0j_tilde,U0j_tilde)+Sigma0_tilde_inv[[h_j]]
        V0j_tilde_inv=solve(V0j_tilde)
        u0[[h_j]][j_j,-1] = mvrnorm(1,V0j_tilde_inv%*%(crossprod(U0j_xi0j_tilde,O0j_tilde)+Sigma0_tilde_inv[[h_j]]%*%(u0_mean_prior_1[[h_j]]+Sigma0_part1[[h_j]]%*%(u0_star[[h_j]][j_j,-1]-u0_mean_prior_2[[h_j]]) ) ),V0j_tilde_inv)
      } else {
        u0[[h_j]][j_j,-1] = mvrnorm(1,u0_mean_prior_1[[h_j]]+Sigma0_part1[[h_j]]%*%(u0_star[[h_j]][j_j,-1]-u0_mean_prior_2[[h_j]]),Sigma0_tilde[[h_j]])
      }
    }
    
    # update u0_star
    for(j in 1:p) {
      h_j = p_mod[j]
      j_j = j_vec[j]
      index_k = a_ind0[[j]]
      lengh_k = length(index_k)
      
      if(lengh_k>0) {
        U0j_star = matrix(u0_star[[h_j]][index_k,-1],nrow= lengh_k )
        xi0j_star = Xi0_star[[h_j]][index_k,j_j]
        kappa0j_star = Kappa_star[[h_j]][index_k,j_j]
        U0j_xi0j_star=U0j_star*xi0j_star
        O0j_star = kappa0j_star/xi0j_star+u00_star[h_j]
        V0j_star = crossprod(U0j_xi0j_star,U0j_star) +Sigma0_star_inv[[h_j]]
        V0j_star_inv=solve(V0j_star)
        u0_star[[h_j]][j_j,-1] = mvrnorm(1,V0j_star_inv%*%(crossprod(U0j_xi0j_star,O0j_star)+Sigma0_star_inv[[h_j]]%*%(u0_mean_prior_2[[h_j]]+Sigma0_part2[[h_j]]%*%(u0[[h_j]][j_j,-1]-u0_mean_prior_1[[h_j]]) ) ),V0j_star_inv)
      } else {
        u0_star[[h_j]][j_j,-1] = mvrnorm(1,u0_mean_prior_2[[h_j]]+Sigma0_part2[[h_j]]%*%(u0[[h_j]][j_j,-1]-u0_mean_prior_1[[h_j]]),Sigma0_star[[h_j]])
      }
    }
    
    # update u1 
    for(j in 1:p) {
      h_j = p_mod[j]
      j_j = j_vec[j]
      index_k = a_ind[[j]]
      lengh_k = length(index_k)
      
      if(lengh_k>0) {
        U1j_tilde = matrix(u1[[h_j]][index_k,-1],nrow = lengh_k)
        xi1j_tilde = Xi1_tilde[[h_j]][index_k,j_j] 
        kappa1j_tilde = Kappa_tilde[[h_j]][index_k,j_j] 
        U1j_xi1j_tilde=U1j_tilde*xi1j_tilde
        O1j_tilde = kappa1j_tilde/xi1j_tilde+u10[h_j]
        V1j_tilde = crossprod(U1j_xi1j_tilde,U1j_tilde) +Sigma1_tilde_inv[[h_j]]
        V1j_tilde_inv=solve(V1j_tilde)
        u1[[h_j]][j_j,-1] = mvrnorm(1,V1j_tilde_inv%*%(crossprod(U1j_xi1j_tilde,O1j_tilde)+Sigma1_tilde_inv[[h_j]]%*%(u1_mean_prior_1[[h_j]]+Sigma1_part1[[h_j]]%*%(u1_star[[h_j]][j_j,-1]-u1_mean_prior_2[[h_j]]) ) ),V1j_tilde_inv)
      } else {
        u1[[h_j]][j_j,-1] = mvrnorm(1,u1_mean_prior_1[[h_j]]+Sigma1_part1[[h_j]]%*%(u1_star[[h_j]][j_j,-1]-u1_mean_prior_2[[h_j]]),Sigma1_tilde[[h_j]])
      }
    }
    
    # update u1_star
    for(j in 1:p) {
      h_j = p_mod[j]
      j_j = j_vec[j]
      index_k = a_ind[[j]]
      lengh_k = length(index_k)
      
      if(lengh_k>0) {
        U1j_star = matrix(u1_star[[h_j]][index_k,-1],nrow= lengh_k )
        xi1j_star = Xi1_star[[h_j]][index_k,j_j]
        kappa1j_star = Kappa_star[[h_j]][index_k,j_j]
        U1j_xi1j_star=U1j_star*xi1j_star
        O1j_star = kappa1j_star/xi1j_star+u10_star[h_j]
        V1j_star = crossprod(U1j_xi1j_star,U1j_star) +Sigma1_star_inv[[h_j]]
        V1j_star_inv=solve(V1j_star)
        u1_star[[h_j]][j_j,-1] = mvrnorm(1,V1j_star_inv%*%(crossprod(U1j_xi1j_star,O1j_star)+Sigma1_star_inv[[h_j]]%*%(u1_mean_prior_2[[h_j]]+Sigma1_part2[[h_j]]%*%(u1[[h_j]][j_j,-1]-u1_mean_prior_1[[h_j]]) ) ),V1j_star_inv)
      } else {
        u1_star[[h_j]][j_j,-1] = mvrnorm(1,u1_mean_prior_2[[h_j]]+Sigma1_part2[[h_j]]%*%(u1[[h_j]][j_j,-1]-u1_mean_prior_1[[h_j]]),Sigma1_star[[h_j]])
      }
    }
    
    
    if(i>burnin) {
      # compute log likelihood
      psi_i = w%*%z+m
      loglikelihood_x = c(loglikelihood_x,n*sum(log(rho))/2-sum((X-psi_i)^2*rho)/2-n*p*log(2*pi)/2) 
      psi = psi + psi_i 
    }
    # update proposal variance
    if(i%%k==0){
      s = i%/%k
      gamma1 = 1/(s+3)^c1
      gamma2 = 10 * gamma1
      r_theta[,,s] = apply(acs_theta[,,(s*k-k+1):(s*k)],c(1,2), mean)
      adaptFactor = exp(gamma2 * (r_theta[,,s] - r_opt))
      theta_sd_prop = theta_sd_prop*adaptFactor
    }
    
  }
  
  psi_mean = psi/(n_iter-burnin)
  rho_mean = rowMeans(RHO[,(burnin+1):n_iter])
  # compute DIC
  DIC =  2*(n*sum(log(rho_mean))/2-sum((X-psi_mean)^2*rho_mean)/2 -n*p*log(2*pi)/2- 2*mean(loglikelihood_x))
  return(list(Z=Z,W=W,DIC=DIC,psi_mean=psi_mean,M=M))
} 
