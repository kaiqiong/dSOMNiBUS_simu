

source("../../functions/BSMethEM_fix_scale_as_init_intelli_fixed_correct_pi_y_use_mean_negative_phiS_remove_return_phiS_truncated_median_directly_With_truncation.R")
#source("../../../functions/BSMethEM_fix_scale_as_init.R")
load("../../functions/BANK1betas.RData")
load("../../functions/BANK1data.RData")
# Load the population parameters
load("../../functions/Settings_Final_MD_one_direction_new_smooth.RData")


sss <- commandArgs(TRUE)

sss <- as.numeric(sss)

M = 100 # MC size
total_index <- 1:M
sep_index <- split(total_index, ceiling(total_index/1))


phi_no <- 2 ; sig_no <-2

phi_v <- c(1,3)
sig_v <- c(0, 1, 3, 9)

phi <- rep(phi_v[phi_no] ,length(pos))
sigma.ee <- sig_v[sig_no]

index_now = expand.grid(1:ncol(BETA.0), 1:length(sep_index))



index_now[sss, 1]

index_now[sss, 2]


sep_index[[index_now[sss, 2]]]

beta.0 <- BETA.0[,index_now[sss,1]]; beta.1 <- BETA.1[,index_now[sss,1]]

my.samp <- 300


BSMethSim_bbinom <-function(n, posit, theta.0, beta, phi, random.eff = F, mu.e=0,
                            sigma.ee=1,p0 = 0.003, p1 = 0.9, X, Z,binom.link="logit"){
  if( !is.matrix((Z)) ) message ("covariate Z is not a matrix")
  #  if( !is.matrix(beta) ) message ("the covariate effect parameter beta is not a matrix")
  
  if( !(nrow(X)==nrow(Z) & nrow(X) == n) ) message("Both X and Z should have n rows")
  if( !(ncol(X)==length(theta.0) & ncol(X) ==nrow (beta) & ncol(X)==length(posit) )) message ("The columns of X should be the same as length of beta theta.0 and posit; They all equals to the number of CpGs")
  if( ncol(beta)!= ncol(Z)) message("beta and Z should have the same dimentions")
  
  # the random effect term
  if(random.eff == T){
    my.e <- rnorm(n, mean=mu.e, sd = sqrt(sigma.ee))
  }else{
    my.e <- rep(mu.e, n)
  }
  
  my.theta <- t(sapply(1:n, function(i){
    theta.0 + rowSums(sapply(1:ncol(Z), function(j){Z[i,j] * beta[,j]})) + my.e[i]
  }))
  
  
  # Transform my.theta to my.pi for each (i, j)
  
  my.pi <- t(sapply(1:nrow(my.theta), function(i){
    #exp(my.theta[i,])/(1+exp(my.theta[i,]))
    binomial(link=binom.link)$linkinv(my.theta[i,])
  }))
  #  Generate S-ij based on the my.pi and my.Z
  #---------------------------------------------------#
  my.S <- my.pi
  for ( i in 1:nrow(my.S)){
    for ( j in 1:ncol(my.S)){
      #my.S[i,j] <- rbinom (1, size = X[i,j], prob= my.pi[i,j])
      my.S[i,j] <-VGAM::rbetabinom(1, size = X[i,j], prob = my.pi[i,j], rho = (phi[j]-1)/(X[i,j]-1) )
    }
  }
  #---------------------------------------------------#
  # Generate Y-ij based on the S-ij and the error rate (1-p1) and p0
  #---------------------------------------------------#
  my.Y <- my.S
  for ( i in 1:nrow(my.Y)){
    for ( j in 1:ncol(my.Y)){
      my.Y[i,j] <- sum(rbinom(my.S[i,j], size =1, prob=p1)) +
        sum(rbinom(X[i,j]-my.S[i,j], size = 1, prob=p0))
    }
  }
  out = list(S = my.S, Y = my.Y, theta = my.theta, pi = my.pi)
}




add_read_depth =4

n.knots = 5; my.p1 = 0.9



covs_use <- c( "cell_type") 

# store results from the other two methods
pos.index <- matrix(NA, nrow=length(pos), ncol =  length(sep_index[[index_now[sss,2]]]))



Beta.est <- array(NA, c(length(pos),  length(sep_index[[index_now[sss,2]]]), 
                        length(covs_use)+1),dimnames = list(NULL, NULL, c("Intercept",covs_use)))



sigma00 <- rep(NA,  length(sep_index[[index_now[sss,2]]]))


sum.est <- array(NA, c( length(covs_use)+1, length(sep_index[[index_now[sss,2]]]) ,3, 3), dimnames = list(c("Intercept",covs_use), NULL, c("EDF","F", "pvalue"),c("SOMNiBUS_Fletcher", "SOMNiBUS_REML", "gam") ))



Beta.se <-array(NA, c(length(pos),  length(sep_index[[index_now[sss,2]]]), length(covs_use)+1, 3),dimnames = list(NULL, NULL, c("Intercept",covs_use), c("SOMNiBUS_Fletcher", "SOMNiBUS_REML", "gam") ))


phi.est <- data.frame(matrix(NA, nrow =  length(sep_index[[index_now[sss,2]]]), ncol = 3))
colnames(phi.est) <- c("SOMNiBUS_Fletcher", "SOMNiBUS_REML", "gam") 
lambda.est <- matrix(NA, nrow =  length(sep_index[[index_now[sss,2]]]), ncol = length(covs_use)+2)


time.0 <- Sys.time()
for(mm in sep_index[[index_now[sss,2]]]){
  set.seed(3432421+mm)
  
  #------------------- Use bootstrap to build the Z matrix-------------------------------------------#
  Z <- data.frame(matrix(NA, nrow= my.samp, ncol = 1)); colnames(Z) <- c( "cell_type")
  
  Z$cell_type <- sample(c(0,1), size = my.samp, replace = T ) # simulate Z from binomial distribution
  Z <-as.matrix(Z);rownames(Z)<- NULL
  samp.Z <- Z
  #Z$disease <- ifelse(Z$disease=="RA",1,0)
  #Z$cell_type <- ifelse(Z$cell_type=="TCELL",1,0)
  
  
  # Use bootstrap to build the read-depth matrix 
  my.X <- matrix(sample(0:1, my.samp*length(pos), replace = T), nrow = my.samp, ncol = length(pos))
  
  #my.X <- matrix(sample(as.vector(dat.use.total), size = nrow(Z)*length(pos), replace = T) ,
  #               nrow = nrow(Z), ncol = length(pos))
  #plot(pos, totalMat[,1])
  
  #plot(dat.use$Position,dat.use$Total_Counts, pch = 19)
  ff = smooth.spline(pos, apply(totalMat, 1, median), nknots = 10)
  
  spacial_shape <- round(predict(ff, pos)$y)
  for ( i in 1:my.samp){
    my.X[i,] <- my.X[i,] + spacial_shape 
    #+ round(10*beta.0) + round(beta.1 * Z[i,1] + beta.2 * Z[i, 2])
  } 
  
  my.X <-  (my.X + add_read_depth)
  
  
  #plot(pos, apply(my.X, 2, median))
  #--- Simulate the data ---#
  
  
  
  if(sigma.ee == 0){
    sim.dat<-BSMethSim_bbinom(n= my.samp, posit = pos, theta.0 =as.matrix(beta.0, ncol=1), beta= cbind(beta.1), phi=phi, 
                              X = my.X, Z =Z,p0 = 0.003, p1 = 0.9,random.eff = F)
  }else{
    
    sim.dat<-BSMethSim_bbinom(n= my.samp, posit = pos, theta.0 =as.matrix(beta.0, ncol=1), beta= cbind(beta.1), phi=phi, 
                              X = my.X, Z =Z,p0 = 0.003, p1 = 0.9,random.eff = T, sigma.ee = sigma.ee)
  }
  
  
  
  X <-my.X; Y <- sim.dat$Y 
  samp.size <- nrow(Y); my.p <- ncol(Y)
  my.span.dat <- data.frame(Y=as.vector(t(Y)), X=as.vector(t(X)), Posit = rep(pos, samp.size),
                            ID = rep(1:samp.size, each=my.p),
                            sapply(1:ncol(Z), function(i){rep(Z[,i], each = my.p)}))
  colnames(my.span.dat)[-(1:4)] <- colnames(Z)
  my.span.dat <- my.span.dat[my.span.dat$X>0,]
  
  # my.span.dat<- data.frame(my.span.dat, null = sample(c(0,1), size = nrow(my.span.dat), replace = T))
  if(ncol(Z)==1){
    Z <- matrix(my.span.dat[,-(1:4)], nrow = nrow(my.span.dat))
  }else{
    Z <- my.span.dat[,-(1:4)]
  }
  
  # -- Fit the model use the proposed method ---#
  
  
  
  colnames(my.span.dat)[1:3] <-  c("Meth_Counts", "Total_Counts", "Position")
  
  
  out <- BSMethEM(data=my.span.dat, n.k = rep(n.knots,ncol(Z)+1 ),epsilon = 10^(-6)/300,p0 = 0.003,
                  p1 = 0.9,maxStep = 500, method="REML", Quasi = T, RanEff = T, reml.scale = F)
  
  imm <- match(mm, sep_index[[index_now[sss,2]]])
  
  # -- Extract the functional parameter estimates ---#
  
  Beta.est[,imm,] <-  out$Beta.out
  # beta.0.est[,mm,jj] <- out$Beta.out[,"Intercept"]
  #  beta.1.est[,mm,jj] <- out$Beta.out[,2]
  
  pos.index[,imm] <- out$uni.pos
  
  Beta.se [,imm,, 1] <- out$SE.out
  Beta.se [,imm,, 2] <- out$SE.out.REML.scale
  Beta.se [,imm,, 3] <- out$SE.out/sqrt(out$phi_fletcher) * sqrt(out$phi_gam)
  
  
  
  # GamObj[[mm]]<-out$FinalGamObj
  
  sum.est[,imm, ,1]<- out$reg.out[1:(length(covs_use)+1),]
  sum.est[,imm, ,2]<- out$reg.out.reml.scale [1:(length(covs_use)+1),]
  sum.est[,imm, ,3]<- out$reg.out.gam[1:(length(covs_use)+1),-3]
  
  pos.index[,imm] <- out$uni.pos
  
  # sum.est.gam[,mm,] <-  mgcv::summary.gam(out$GamObj)$s.table[1:(length(covs_use)+2),2:4]
  
  lambda.est[imm,] <- out$lambda[1:(length(covs_use)+2)] 
  
  phi.est$SOMNiBUS_Fletcher[imm] <- out$phi_fletcher
  phi.est$SOMNiBUS_REML [imm] <- out$phi_reml
  phi.est$gam[imm] <- out$phi_gam
  sigma00[imm] <- out$sigma00
  
  # phi.est.gam[mm]<-out$phi_gam
  
  print(mm)
  print(Sys.time()-time.0)

}

save.image(paste0("Results/S", index_now[sss,1],"Batch", index_now[sss,2],  "Samp", my.samp, "Simu", M,".RData")) 



