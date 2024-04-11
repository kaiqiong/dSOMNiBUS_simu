
sss <- commandArgs(TRUE)

source("../../functions/BSMethEM_fix_scale_as_init_intelli_fixed_correct_pi_y_use_mean_negative_phiS_remove_return_phiS_truncated_median_directly_With_truncation.R")

# real data results

load("../../functions/BANK1data.RData")
#RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0,])
#out_quasi = BSMethEM(data=RAdat.f, n.k = rep(5,3), p0 = 0.003, p1 = 0.9, Quasi = T, RanEff = F)

#out_quasi_re = BSMethEM(data=RAdat.f, n.k = rep(5,3), p0 = 0.003, p1 = 0.9, Quasi = T, RanEff = T)
#res_som <- out$reg[-1,]



library(magrittr)

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

load("../../functions/BANK1betas.RData")
beta.0 <- beta.0/4 + 2
beta_all <- cbind(beta.0, beta.1, beta.2, rep(0, length(beta.0)))

beta_all[,1] <- beta_all[,1]-3.5


phi <- rep(3,length(pos))


my.samp <- 100

#-----------------------------------------------------------------------#
n.k = rep(5, 4)#my.p1 = 0.9


add_read_depth =4
time.0 <- Sys.time()


covs_use <- c("disease", "cell_type", "NullZ") 
#--------


M = 1000 # MC size
pos.index <- matrix(NA, nrow=length(pos), ncol = M)



Beta.est <- array(NA, c(length(pos), M, length(covs_use)+1),dimnames = list(NULL, NULL, c("Intercept",covs_use)))

sum.est <- array(NA, c( length(covs_use)+1,M ,3), dimnames = list(c("Intercept",covs_use), NULL, c("EDF","F", "pvalue")))

Beta.se <- Beta.est

phi.est <- rep(NA, M)
phi.est.gam <- phi.est
sum.est.gam <- sum.est

sigma00 <- rep(NA, M)




total_index <- 1:M
sep_index <- split(total_index, ceiling(total_index/2))

#file_names <- list.files(pattern = "Res_all_Samp100Simu1000.RData")
#fail_ids <- which(1:1000 %in% as.numeric(gsub(".*?([0-9]+).*", "\\1", file_names))==FALSE)

for ( mm in sep_index[[sss]]){


  set.seed(3432421+mm)
  
  
  Z <- data.frame(matrix(NA, nrow= my.samp, ncol = 2)); colnames(Z) <- c("disease", "cell_type")
  
  #Z$disease <- sample(my.pheno$disease, size = my.samp, replace = T)
  #Z$cell_type<- sample(my.pheno$cell_type, size = my.samp, replace = T)
  #table(my.pheno$disease)
  
  Z$disease <- sample(c(0,1), size = my.samp, replace = T, 
                      prob = c(sum(my.pheno$disease =="control"),sum(my.pheno$disease =="RA") )/nrow(my.pheno) )
  Z$cell_type <- sample(c(0,1), size = my.samp, replace = T, 
                        prob = c(sum(my.pheno$cell_type =="MONO"),sum(my.pheno$cell_type =="TCELL") )/nrow(my.pheno) )
  Z$NullZ <-  sample(c(0,1), size = nrow(Z), replace = T)
  #Z$disease <- ifelse(Z$disease=="RA",1,0)
  #Z$cell_type <- ifelse(Z$cell_type=="TCELL",1,0)
  
  Z <-as.matrix(Z);rownames(Z)<- NULL
  
  samp.Z <- Z
  #-----------------------------
  # Build a read-depth matrix which sort of preserve the dependence structure in read-depths
  
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
  
  sim.dat<-BSMethSim_bbinom(n= my.samp, posit = pos, theta.0 =beta_all[,1], beta= beta_all[,-1], phi=phi, 
                            X = my.X, Z =Z,p0 = 0.003, p1 = 0.9,random.eff = T,  sigma.ee=3)
  
  
  #--- Organize the data before EM-smooth ---#
  X <-my.X; Y <- sim.dat$Y 
  samp.size <- nrow(Y); my.p <- ncol(Y)
  
  dat.use <- data.frame(Meth_Counts=as.vector(t(Y)), 
                        Total_Counts=as.vector(t(X)), 
                        Position = rep(pos, samp.size),
                        ID = rep(1:samp.size, each=my.p))
  
  
  covs_use <- colnames(samp.Z)
  for( j in 1:length(covs_use)){
    dat.use <- data.frame(dat.use, rep(samp.Z[,j], each = my.p))
  }
  colnames(dat.use)[-c(1:4)] <- covs_use
  
  dat.use <- dat.use[dat.use$Total_Counts>0,]
  
  #my.span.dat<- data.frame(my.span.dat, null = sample(c(0,1), size = nrow(my.span.dat), replace = T))
  
  #Z <- dat.use[,-c(1:4)]
  
  # -- Fit the model use the proposed method ---#
  

  
  out <- BSMethEM(data=dat.use, n.k = n.k,epsilon = 10^(-6)/300,p0 = 0.003,
                  p1 = 0.9,maxStep = 500, method="REML", Quasi = T, RanEff = T, reml.scale = F)
  
  
  print(Sys.time()-time.0)
  
 
  save(Y, out, file =paste0("Results/Res_",mm,"V0Res_all_Samp",my.samp, "Simu", M, ".RData") )
}


#save.image(file = paste0("V2Res_all_Samp",my.samp, "Simu", M, ".RData"))




