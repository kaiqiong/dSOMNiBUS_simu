#' @title A smoothed-EM algorithm to estimate covariate effects and test regional association in Bisulfite Sequencing-derived methylation data
#'
#' @description This function fits a binomial regression model where the outcome - methylated reads- are contaminated by known error rates \code{p0} and \code{p1} and the covariate effects are smoothly varying across genomic positions. The functional parameters of the smooth covariate effects are first represented by a linear combination of a bunch of restricted cubic splines (with dimention \code{n.k}), and a smoothness penalization term which depends on the smoothing parameters \code{lambdas} is also added to control smoothness.
#' @description The estimation is performed by an iterated EM algorithm. Each M step constitutes an outer Newton's iteration to estimate smoothing parameters \code{lambdas} and an inner P-IRLS iteration to estimate spline coefficients \code{alpha} for the covariate effects. Currently, the computation in the M step depends the implementation of \code{gam()} in package \code{mgcv}.
#' @param data a data frame with rows as individual CpGs appeared in all the samples. The first 4 columns should contain the information of `Meth_Counts` (methylated counts), `Total_Counts` (read depths), `Position` (Genomic position for the CpG site) and `ID` (sample ID). The covariate information, such as disease status or cell type composition are listed in column 5 and onwards.
#' @param n.k a vector of basis dimensions for the intercept and individual covariates. \code{n.k} specifies an upper limit of the degrees of each functional parameters.
#' @param p0 the probability of observing a methylated read when the underlying true status is unmethylated.
#' @param p1 the probability of observing a methylated read when the underlying true status is methylated.
#' @param covs a vector of covariate names. The covariates with names in \code{covs} will be included in the model and their covariate effects will be estimated.
#' @param Quasi whether a Quasi-likelihood estimation approach will be used
#' @param epsilon numeric; stopping criterion for the closeness of estimates of spline coefficients from two consecutive iterations.
#' @param epsilon.lambda numeric; stopping criterion for the closeness of estimates of smoothing parameter \code{lambda} from two consecutive iterations.
#' @param maxStep the algorithm will step if the iteration steps exceed \code{maxStep}
#' @param detail indicate whether print the number of iterations
#' @param binom.link the link function used in the binomial regression model; the default is the logit link
#' @param method the method used to estimate the smoothing parameters. The default is the "REML" method which is generally better than prediction based criterion \code{GCV.cp}
#' @param RanEff whether sample-level random effects are added or not
#' @param reml.scale whether a REML-based scale (dispersion) estimator is used. The default is Fletcher-based estimator
#' @param scale nagative values mean scale paramter should be estimated; if a positive value is provided, a fixed scale will be used.
#' @return This function return a \code{list} including objects:
#' \itemize{
#' \item \code{est}: estimates of the spline basis coefficients alphas
#' \item \code{lambda}: estimates of the smoothing parameters for each functional paramters
#' \item \code{est.pi}: predicted methylation lelves for each row in the input \code{data}
#' \item \code{ite.points}: estimates of \code{est}, \code{lambda} at each EM iteration
#' \item \code{cov1}: estimated variance-covariance matrix of the basis coefficients alphs
#' \item \code{reg.out}: output from the reginal zero effect test
#' \item \code{chi.sq}: chi-square statistics for each covariates (including intercept) w.r.t the zero regional effect tests
#' \item \code{pvalue}: the p value from the regional tests for each covariate
#' \item \code{pvalue.log}: the log of \code{pvalue}
#' \item \code{SE.out} a matrix of the estimated pointwise Standard Errors (SE); number of rows are the number of unique CpG sites in the input data and the number of columns equal to the total number of covariates fitted in the model (the first one is the intercept)
#' \item \code{SE.pos} the genomic postions for each row of CpG sites in the matrix \code{SE.out}
#' \item \code{Beta.out} a matrix of the estimated covariate effects beta(t), here t denots the genomic positions.
#' \item \code{ncovs} number of functional paramters in the model (including the intercept)
#' }
#' @author  Kaiqiong Zhao
#' @seealso  \link[mgcv]{gam}
#' @examples #------------------------------------------------------------#
#' data(RAdat)
#' head(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0,])
#' out <- BSMethEM(data=RAdat.f, n.k = rep(5,3), p0 = 0.003307034, p1 = 0.9,
#' epsilon = 10^(-6), epsilon.lambda = 10^(-3), maxStep = 200, detail=FALSE)
#' @importFrom mgcv gam
#' @importFrom mgcv predict.gam
#' @importFrom mgcv model.matrix.gam
#' @importFrom mgcv s
#' @importFrom Matrix bdiag
#' @importFrom stats as.formula binomial pchisq rbinom rnorm
#' @export
BSMethEM = function (data, n.k, p0 = 0.003, p1 = 0.9, Quasi = TRUE, epsilon = 10^(-6),  epsilon.lambda = 10^(-3), maxStep = 200,  detail=FALSE, binom.link = "logit",method="REML", covs = NULL, RanEff = T, reml.scale=F, scale = -2){

  n.k <<-n.k # an error of 'n.k is not found' would appear if without this golable environment assignment; so I save n.k in a parent scope
  data <- data.frame(data)

  if(is.factor(data$Position) ) {
    #message("The Position in the data set should be numeric other than a factor")
    data$Position <- as.numeric(as.character(data$Position))
  }

  if( any(!c("Meth_Counts", "Total_Counts", "Position")%in% colnames(data))) stop('Please make sure object "data" have columns named as "Meth_Counts", "Total_Counts" and "Position" ')


  colnames(data)[match(c("Meth_Counts", "Total_Counts", "Position") ,colnames(data)) ] <- c("Y", "X", "Posit")

  if(is.null(covs)){
    Z <- as.matrix(data[,-(1:4)], ncol = ncol(data)-4)
    colnames(Z) <- colnames(data)[-c(1:4)]
  }else{
    id <- match(covs, colnames(data))
    if(any(is.na(id))){
      stop(paste(covs[is.na(id)], " is not found in the input data frame"))
    }else{
      Z <- as.matrix(data[,id], ncol = length(id))
      colnames(Z) <- covs
    }
  }

  if(length(n.k)!= (ncol(Z) +1)) stop('The length of n.k should equal to the number of covariates plus 1 (for the intercept)')
  if(any(data$X == 0)) stop('The rows with Total_Counts equal to 0 should be deleted beforehand')
  if(any(is.na(Z))) stop('The covariate information should not have missing values')
  if(any(!is.numeric(Z))) stop('Please transform the covariate information into numeric values, eg. use dummy variables for the categorical covariates')

  # The smoothing formula corresponding to the Z
  formula.z.part <- sapply(1:(ncol(Z)), function(i){
    paste0("s(Posit, k = n.k[",i+1,"], fx=F, bs=\"cr\", by = Z[,", i, "])"  ) } )
  my.covar.fm <- paste(c("s(Posit, k=n.k[1], fx=F, bs=\"cr\")", formula.z.part), collapse="+")
  if(RanEff){
    my.covar.fm <- paste0(my.covar.fm, "+ s(ID, bs = \"re\")" )
    data$ID <- as.factor(data$ID)
  }
  # Fit gam for the initial value
  if(Quasi){
    gam.int <- mgcv::gam(as.formula( paste0("Y/X ~", my.covar.fm)), family =quasibinomial(link = binom.link),weights=X, data = data, method = method, scale = scale)
  }else{
    gam.int <- mgcv::gam(as.formula( paste0("Y/X ~", my.covar.fm)), family =binomial(link = binom.link),weights=X, data = data, method = method, scale = scale)
  }
  # Estimates
  old.pi.ij <- gam.int$fitted.values;old.par <-gam.int$coefficients;lambda <- gam.int$sp
  #phi_fletcher = summary(gam.int)$dispersion
  p_res <- residuals(gam.int, type ='pearson')
  d_res <- residuals(gam.int, type ="deviance")
  edf.out <- gam.int$edf
  edf1.out <- gam.int$edf1

  N <- length(unique(data$ID))
  # Note: this phi_fletcher can be also self-calculated
  if(Quasi & scale<=0){   # calculate the estimate of phi if Quasi = T and scale is unknown
    my_s <- (1-2*old.pi.ij)/(data$X*old.pi.ij*(1-old.pi.ij))*(data$Y-data$X*old.pi.ij)#*sqrt(data$X)

    # Note: the estimator implemented in the mgcv calculated my_s with an additional multiplier sqrt(data$X)
    # But from the paper there shouldn't be this one
    #phi_p = sum( p_res^2)/(length(data$Y) - sum(edf.out))

    # Feb 21, 2020  use edf1 to calculate the pearson's dispersion estimate not the edf; because I will use the asympototic chi-square dist. of phi.est
    # March 3, 2020, the dispersion parameter estimated in gam use the edf.out instead of edf1.out

   # my_s_gam <- my_s*sqrt(data$X)

    phi_p = sum( p_res^2)/(length(data$Y) - sum(edf.out))

    phi_fletcher <- phi_p/(1+mean(my_s))
   # phi_gfletcher <- phi_p/(1+mean(my_s_gam))  # This gives the exact results as mgcv::gam -dispersion

    if(reml.scale){phi_fletcher <- gam.int$reml.scale}
    #phi_fletcher <- out$phi_fletcher
    #phi_fletcher <- summary(gam.int)$dispersion
  }

  if(!Quasi){
    phi_fletcher <- 1
  }

  if(scale>0 ){
    phi_fletcher = scale
  }



  if(p0==0 & p1==1){

    #phi_S = phi_fletcher; sd_phi_S = 0
    out <- list( pi.ij = gam.int$fitted.values, par = gam.int$coefficients,
                 lambda = gam.int$sp, edf1 = gam.int$edf1, pearson_res = p_res, deviance_res=d_res,
                 edf=gam.int$edf, phi_fletcher= phi_fletcher, GamObj = gam.int )
    new.par<-out$par; new.lambda <- out$lambda;new.pi.ij <- out$pi.ij

    Est.points <-  c(new.par, new.lambda, phi_fletcher)
    new.p0 = 0; new.p1 = 1

  }else{
    
    # transform the phi_fletcher (for y) to the phi_fletcher (for S) by a transformation -- Oct, 15, 2020
    #piy <- old.pi.ij
    
    #piS <- (old.pi.ij-p0)/(p1-p0)
      #old.pi.ij*p1 + (1-old.pi.ij)*p0
    #phi_S <- (phi_fletcher-1)/(p0-p1)^2*(old.pi.ij*(1-old.pi.ij))/(piS*(1-piS)) + 1
    
    
    #piy_used <- old.pi.ij[which((old.pi.ij>p0+0.05) & (old.pi.ij<p1-0.05)) ]
    # we don't need use all the piy
    #piy_used <- old.pi.ij[-which((old.pi.ij<p0)|old.pi.ij>p1)]
    #piy_used <- old.pi.ij[-which((abs(old.pi.ij-p0)<0.001)|(abs(old.pi.ij-p1)<0.001))]
    #piy_used <- old.pi.ij
    #phi_S <- (phi_fletcher-1)*(piy_used*(1-piy_used))/((piy_used-p0)*(p1-piy_used))+1
    
    #see = seq(0.1, 0.8, 0.01)
    
    print("phi when ignoring the error")
    print(phi_fletcher)
    
    #-----------------
    # option 1 - use mean
    #----------------
    mean_part <- mean(((old.pi.ij-p0)*(p1-old.pi.ij))/(old.pi.ij*(1-old.pi.ij)) )
    print("mean of part")
    print(mean_part)
    
    phi_S_mean <- (phi_fletcher-1)/mean_part+1
    
    print("phi for S after transformation use mean")
    print(phi_S_mean)
    
    #-----------------
    # option 2 - use median
    #----------------
    median_part <- median(((old.pi.ij-p0)*(p1-old.pi.ij))/(old.pi.ij*(1-old.pi.ij)) )
    print("median of part")
    print(median_part)
    phi_S_median <- (phi_fletcher-1)/median_part+1
    
    print("phi for S after transformation use median")
    print(phi_S_median)
    print("the SD of pi_y")
    print(sd(old.pi.ij))
    
    #-----------------
    # option 3 - use mean -- after removing the ijs such that old.pi.ij-p0<0 or p1-old.pi.ij<0
    #----------------
    
    piy_prop0 = mean(old.pi.ij-p0<0)
    piy_prop1 = mean(p1-old.pi.ij<0)


    print("How many piY less than p0")
    print(mean(old.pi.ij-p0<0))
    
  
    print("How many piY greater than p1")
    print(mean(p1-old.pi.ij<0))
    
    old.pi.ij_truc = old.pi.ij[(old.pi.ij>p0) & (old.pi.ij<p1)] 
    mean_part_truc = mean(((old.pi.ij_truc-p0)*(p1-old.pi.ij_truc))/(old.pi.ij_truc*(1-old.pi.ij_truc)) )
    median_part_truc = median(((old.pi.ij_truc-p0)*(p1-old.pi.ij_truc))/(old.pi.ij_truc*(1-old.pi.ij_truc)) )

 phi_S_mean_truc <- (phi_fletcher-1)/mean_part_truc +1
 phi_S_median_truc <- (phi_fletcher-1)/median_part_truc +1
    #save(phi_S, phi_fletcher, file = "see.RData")
    print("phi_S_mean_truc")
    print(phi_S_mean_truc)

    print("phi_S_median_truc")
    print(phi_S_median_truc)

  phi_S = phi_S_mean_truc
  #phi_S = phi_S_mean

  
    
    if(phi_S >0){
      phi_fletcher <- phi_S
    }
    
    print("The final dispersion parameter")
    print(phi_fletcher)

    #---------
    old.p0 <- p0; old.p1 <- p1
    out <-  BSMethEMUpdate (data, old.pi.ij, p0 = old.p0, p1 = old.p1, n.k=n.k, binom.link = binom.link, method = method, Z = Z, my.covar.fm = my.covar.fm, Quasi = Quasi, scale = phi_fletcher ,RanEff = RanEff, N = N)
    new.par<-out$par; new.lambda <- out$lambda;new.pi.ij <- out$pi.ij
    new.phi <- out$phi_fletcher
    new.p1 <- out$p1_est
    new.p0 <- out$p0_est
    i <- 1; Est.points <- rbind(c(old.par, lambda, phi_fletcher, old.p0, old.p1, NA, NA, NA, NA, NA),
                                c(new.par, new.lambda, new.phi, new.p0, new.p1, sqrt(sum((new.pi.ij - old.pi.ij)^2)),  
                                  sqrt(sum((c(old.p0, old.p1) - c(new.p0, new.p1) )^2)), out$Q1, out$Qp1, out$Qp0 ))


  
    # Do the iteration
    # The stopping criterion is that estimator a and b are close enough
    # I exclude the criterio that lambda-a and lambda-b are close
    # while ( sum((old.par - new.par)^2) > (epsilon/15 * length(new.par)) & i < maxStep & sum((lambda-new.lambda)^2) > epsilon.lambda){

    while ( sqrt(sum((new.pi.ij - old.pi.ij)^2)) > epsilon & i < maxStep ){
      i <- i +1;
      old.par <- new.par
      old.pi.ij <- new.pi.ij
      old.p0 <- new.p0
      old.p1 <- new.p1
      out <-  BSMethEMUpdate (data, old.pi.ij, p0 = old.p0, p1 = old.p1 ,binom.link = binom.link, method = method, Z = Z, my.covar.fm = my.covar.fm, Quasi = Quasi, scale = phi_fletcher,RanEff = RanEff, N = N)
      new.par<-out$par; new.lambda <- out$lambda;new.pi.ij <- out$pi.ij;new.phi <- out$phi_fletcher
      new.p1 <- out$p1_est
      new.p0 <- out$p0_est
      Est.points <- rbind(Est.points, c(new.par,new.lambda, new.phi, new.p0, new.p1, sqrt(sum((new.pi.ij - old.pi.ij)^2)),
                                        sqrt(sum((c(old.p0, old.p1) - c(new.p0, new.p1) )^2)), out$Q1, out$Qp1, out$Qp0 
                                        ))
      if(detail){
        print(paste0("iteration", i))
      }
    }
    # Update


  }
  # Effective degrees of freedom:  edf1 -- good for chisquare test and p value calculation tr(2A - A^2)
  edf1.out <- out$edf1
  # Effective degree of freedom: edf --trace of the hat matrix
  edf.out <- out$edf
  # the residuals degrees of freedom
  resi_df <- nrow(data) - sum(edf.out)
  # Pearson Residuals
  p_res <- out$pearson_res
  # Estimated dispersion paramters (Fletcher adjusted)
  #phi_fletcher <- out$phi_fletcher

  phi_fletcher <- out$phi_fletcher

  GamObj <- out$GamObj



  #--------------------------------------------
  # Calculate SE of alpha
  #-------------------------------------------
  # the part to calculate the SE of alpha
  my.design.matrix <-  mgcv::model.matrix.gam(GamObj) # the model matrix for the GamObj and the FinalGamObj is the same. the difference is only on the outcomes
  Y <- data$Y ; X <- data$X
  pred.pi <- new.pi.ij
  #N <- length(unique(data$ID))
  #---------------- this part is inside  the trycatch

  phi_reml <- GamObj$reml.scale
 # phi_gam_default <- GamObj$scale
  #-------------------------------------------------------------------------
  # from the variance-covariance of alpha, to report
  # 1. var(beta_p(t))
  # 2. Report the chi-square statistic for each covariate
  # 3. Report the p-value for each covariate
  #-------------------------------------------------------------------------

  # ------ 3: estimate of beta(t) --- #
  uni.pos <- unique(data$Posit); uni.id <- match(uni.pos, data$Posit)
  BZ <- my.design.matrix[uni.id, 1:n.k[1]]
  BZ.beta = lapply(1:ncol(Z), function(i){mgcv::smooth.construct(mgcv::s(Posit, k = n.k[i+1],
                                                                         fx = F, bs = "cr") ,
                                                                 data = data[uni.id,],
                                                                 knots = NULL)$X })

  # Use PredictMat to get BZ.beta





  cum_s <-cumsum(n.k)
  alpha.sep <- lapply(1:ncol(Z), function(i){new.par[ (cum_s[i]+1):cum_s[i+1]]}); alpha.0 <- new.par[1:n.k[1]]

  Beta.out<-  cbind( BZ %*% alpha.0, sapply(1:ncol(Z), function(i){BZ.beta[[i]] %*% alpha.sep[[i]]}))

  colnames(Beta.out) <- c("Intercept", colnames(Z))

  # Simply check to see
  # tryCatch({
  #----------------------------------------------------------------
  # calculate var_cov (for alpha & beta)  from the Hessian matrix
  #----------------------------------------------------------------

  H = Hessian(w_ij =pred.pi * (1-pred.pi)/phi_fletcher,        new.par, new.lambda, X, Y, my.design.matrix, GamObj, Z,pred.pi, new.p0, new.p1, disp_est = phi_fletcher, RanEff = RanEff, N=N)

  var.cov.alpha <- solve(-H)
  var.alpha.0 <- var.cov.alpha[1:n.k[1], 1:n.k[1]]
  var.alpha.sep <- lapply(1:ncol(Z), function(i){var.cov.alpha[ (cum_s[i]+1):cum_s[i+1],(cum_s[i]+1):cum_s[i+1]]})  # SE of the effect of Zs [beta.1(t), beta.2(t), beta.3(t) ...]
  #var.cov.alpha  <- matlib::Ginv(-H)
  # solve and Ginv give very similar results, but MASS is very different from the other two

  #var.cov.alpha <- MASS::ginv(-H)
  # SE.out --- pointwise standard deviation


  # A more efficient way to calculate SE rowsum(A*B) is faster than diag(A %*% t(B))
  # sqrt(pmax(0,rowSums((BZ.beta_now[[2]]%*%var.cov.alpha[7:12,7:12,drop=FALSE])*BZ.beta_now[[2]])))
  SE.out <- cbind(sqrt(pmax(0, rowSums( ( BZ %*% var.alpha.0) * BZ))),
                  sapply(1:ncol(Z), function(i){sqrt(pmax(0,rowSums((BZ.beta[[i]] %*% var.alpha.sep[[i]])*BZ.beta[[i]])))
                  }) )

  rownames(SE.out) <- uni.pos; colnames(SE.out) <-  c("Intercept", colnames(Z))
  SE.pos <- uni.pos


  #SE.out <- vector(mode = "list", length = (ncol(Z)+1))
  #names(SE.out) <- c("Intercept", colnames(Z))


  #var.beta.0 <-  BZ %*% var.alpha.0 %*% t(BZ)
  #var.beta <- lapply(1:ncol(Z), function(i){BZ.beta[[i]] %*% var.alpha.sep[[i]] %*% t(BZ.beta[[i]])})

  #SE.out <- cbind(sqrt(diag(var.beta.0)), sapply(1:ncol(Z), function(i){sqrt(diag(var.beta[[i]]))}))

  SE.out.REML.scael <- SE.out/sqrt(phi_fletcher) * sqrt(phi_reml)  #phi_reml
  #SE.out.gFletcher <- SE.out/sqrt(phi_fletcher) * sqrt(phi_gfletcher)

  #----------------------------------------------------------------
  # calculate the region-based statistic from the testStats function in mgcv
  #----------------------------------------------------------------

  # A more efficient way to extract design matrix. use a random sample of rows of the data to reduce the computational cost
  if(RanEff) {
    re.test=T
  }else{
    re.test= F
  }

  if (!is.null(GamObj$R))  X_d <- GamObj$R else {
    sub.samp <- max(1000,2*length(GamObj$coefficients))
    if (nrow(GamObj$model)>sub.samp) { ## subsample to get X for p-values calc.
      seed <- try(get(".Random.seed",envir=.GlobalEnv),silent=TRUE) ## store RNG seed
      if (inherits(seed,"try-error")) {
        runif(1)
        seed <- get(".Random.seed",envir=.GlobalEnv)
      }
      kind <- RNGkind(NULL)
      RNGkind("default","default")
      set.seed(11) ## ensure repeatability
      ind <- sample(1:nrow(GamObj$model),sub.samp,replace=FALSE)  ## sample these rows from X
      X_d <- predict(GamObj,GamObj$model[ind,],type="lpmatrix")
      RNGkind(kind[1],kind[2])
      assign(".Random.seed",seed,envir=.GlobalEnv) ## RNG behaves as if it had not been used
    } else { ## don't need to subsample
      X_d <- model.matrix(GamObj)
    }
    X_d  <- X_d [!is.na(rowSums(X_d )),] ## exclude NA's (possible under na.exclude)
  } ## end if (m>0)

  #p : estimated paramters --- alpha.0, alpha.seq
  #Xt: the design matrix --- my.design.matrix
  #V: estimated variance  matrix
  ## on entry `rank' should be an edf estimate
  ## 0. Default using the fractionally truncated pinv.
  ## 1. Round down to k if k<= rank < k+0.05, otherwise up.
  ## res.df is residual dof used to estimate scale. <=0 implies
  ## fixed scale.

             s.table <- BSMethEM_summary(GamObj, var.cov.alpha, new.par, edf.out, edf1.out, X_d, resi_df, Quasi, scale, RanEff,re.test,Z)
  s.table.REML.scale <- BSMethEM_summary(GamObj, var.cov.alpha/phi_fletcher*phi_reml, new.par, edf.out, edf1.out, X_d, resi_df, Quasi, scale, RanEff,re.test,Z)
  #var_out = list(cov1 = var.cov.alpha, reg.out = reg.out,  SE.out = SE.out, uni.pos = SE.pos,  pvalue = pvalue , ncovs = ncol(Z)+1)
  #Est_out = list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij, ite.points = Est.points,
  #               Beta.out = Beta.out, phi_fletcher = phi_fletcher)


  if(RanEff){
    sigma00 <- GamObj$reml.scale/GamObj$sp["s(ID)"]
  }else{
    sigma00 <- NA
  }

  reg.out.gam = summary(GamObj)$s.table

  if(p0==0 & p1==1){
    return(out=list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij,
                    Beta.out = Beta.out,
                    phi_fletcher = phi_fletcher,
                    phi_reml = phi_reml,
                    phi_gam = GamObj$scale,
                    reg.out = s.table,
                    reg.out.reml.scale = s.table.REML.scale,
                    cov1 = var.cov.alpha,
                    reg.out.gam= reg.out.gam,
                    SE.out = SE.out,
                    SE.out.REML.scale = SE.out.REML.scael,
                    uni.pos = SE.pos,
                    ncovs = ncol(Z)+1 , ite.points = Est.points, sigma00 = sigma00))
  }else{
  
  return(out=list(est = new.par, lambda = new.lambda, est.pi = new.pi.ij,
                  Beta.out = Beta.out,
                  phi_fletcher = phi_fletcher,
                  phi_reml = phi_reml,
                  phi_gam = GamObj$scale,
                  reg.out = s.table,
                  reg.out.reml.scale = s.table.REML.scale,
                  cov1 = var.cov.alpha,
                  reg.out.gam= reg.out.gam,
                  SE.out = SE.out,
                  SE.out.REML.scale = SE.out.REML.scael,
                  uni.pos = SE.pos,
                  ncovs = ncol(Z)+1 , ite.points = Est.points, sigma00 = sigma00,
                  new.p0=new.p0, new.p1 = new.p1,
                  piy_prop0= piy_prop0,
                  piy_prop1= piy_prop1,
                  phi_S_mean = phi_S_mean,
                  phi_S_median=phi_S_median,
                  phi_S_mean_truc=phi_S_mean_truc,
                  phi_S_median_truc=phi_S_median_truc))
  }

}
Hessian <- function(w_ij, new.par, new.lambda, X, Y, my.design.matrix, gam.int, Z,pred.pi, p0, p1, disp_est, RanEff, N){

  # Q1: the second partial derivative w.r.t alpha^2
  # Q2: the second derivative w.r.t alpha & alpha_star
  res <- outer( 1:length(new.par), 1:length(new.par), Vectorize(function(l,m)  sum(-X * w_ij * my.design.matrix[, m]*my.design.matrix[,l] )) )
  smoth.mat <-lapply(as.list(1:(ncol(Z)+1)), function(i){gam.int$smooth[[i]]$S[[1]] * new.lambda[i]})  # extract the penalty matrix

  smoth.mat[[length(smoth.mat) + 1]] <- 0  # assume the lambda for the constant of the intercept is 0 -- no penalization
  if(RanEff){
    #new.lambda[ncol(Z)+2]
    #smoth.mat[[length(smoth.mat) + 1]] <- diag(N) /(re_sd^2) # NOTE that the smoothing matrix for the random effects is not identity remember to multiple lambda (1/sigma_0^2)
    smoth.mat[[length(smoth.mat) + 1]] <- diag(N) *new.lambda[ncol(Z)+2] # !!!! Otherwise, we get very wide CI
    span.penal.matrix <- as.matrix(Matrix::bdiag( smoth.mat[c(length(smoth.mat)-1, (1:(length(smoth.mat)-2)), length(smoth.mat))] ))
  }else{

    span.penal.matrix <- as.matrix(Matrix::bdiag( smoth.mat[c(length(smoth.mat), (1:(length(smoth.mat)-1)))] ))
  }

  Q1_with_lambda <- res - span.penal.matrix/disp_est
  Q1_no_lambda <- res

  Q2 <- outer(1:length(new.par), 1:length(new.par), Vectorize(function(l,m){
    term1 <- Y*p1*p0/(p1*pred.pi + p0 * (1-pred.pi))^2 + (X-Y)*(1-p1)*(1-p0)/((1-p1)*pred.pi + (1-p0) * (1-pred.pi))^2
    sum(term1 * w_ij  * my.design.matrix[,m] * my.design.matrix[,l])
  }))

  return(Q1_with_lambda + Q2)

}
BSMethEM_summary <- function(GamObj, var.cov.alpha, new.par, edf.out, edf1.out, X_d, resi_df, Quasi, scale, RanEff,re.test, Z){
  ii <- 0
  m <- length(GamObj$smooth)
  df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, m)
  for ( i in 1:m){
    start <- GamObj$smooth[[i]]$first.para; stop <- GamObj$smooth[[i]]$last.para

    V <- var.cov.alpha[start:stop,start:stop,drop=FALSE] ## Bayesian

    p <- new.par[start:stop]  # params for smooth
    edfi <- sum(edf.out[start:stop]) # edf for this smooth
    ## extract alternative edf estimate for this smooth, if possible...
    edf1i <-  sum(edf1.out[start:stop])
    Xt <- X_d [,start:stop,drop=FALSE]
    fx <- if (inherits(GamObj$smooth[[i]],"tensor.smooth")&&
              !is.null(GamObj$smooth[[i]]$fx)) all(GamObj$smooth[[i]]$fx) else GamObj$smooth[[i]]$fixed
    if (!fx&&GamObj$smooth[[i]]$null.space.dim==0&&!is.null(GamObj$R)) { ## random effect or fully penalized term
      res <- if (re.test) mgcv:::reTest(GamObj,i) else NULL     # Test the mth smooth for equality to zero  (m is not the RE term)
      ## and accounting for all random effects in model
    } else { ## Inverted Nychka interval statistics

      if (Quasi) rdf <-  resi_df else rdf <- -1
      res <- testStat(p,Xt,V,min(ncol(Xt),edf1i),type=0,res.df = rdf)
    }

    if (!is.null(res)) {
      ii <- ii + 1
      df[ii] <- res$rank
      chi.sq[ii] <- res$stat
      s.pv[ii] <- res$pval
      edf1[ii] <- edf1i
      edf[ii] <- edfi
      names(chi.sq)[ii]<- GamObj$smooth[[i]]$label
    }

    if (ii==0) df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, 0) else {
      df <- df[1:ii];chi.sq <- chi.sq[1:ii];edf1 <- edf1[1:ii]
      edf <- edf[1:ii];s.pv <- s.pv[1:ii]
    }
    if (!Quasi || scale>=0) {
      s.table <- cbind(edf, df, chi.sq, s.pv)
      dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "Chi.sq", "p-value"))
    } else {
      s.table <- cbind(edf, df, chi.sq/df, s.pv)
      dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "F", "p-value"))
    }
  }
  s.table <- s.table[, -1]
  if(RanEff){ rownames(s.table) <- c("Intercept", colnames(Z), "ID")
  }else{
    rownames(s.table) <- c("Intercept", colnames(Z))}
  colnames(s.table)[1] <- "EDF"
  s.table
}

## testStat function is directly copied from the mgcv package
## Implements Wood (2013) Biometrika 100(1), 221-228
## The type argument specifies the type of truncation to use.


# res.dif = -1 when no dispersion

# res.dif = residual.df = residual.df<-length(object$y)-sum(object$edf) when est.dispersion = TRUE

testStat <- function(p,X,V,rank=NULL,type=0,res.df= -1) {
  ## Implements Wood (2013) Biometrika 100(1), 221-228
  ## The type argument specifies the type of truncation to use.
  ## on entry `rank' should be an edf estimate
  ## 0. Default using the fractionally truncated pinv.
  ## 1. Round down to k if k<= rank < k+0.05, otherwise up.
  ## res.df is residual dof used to estimate scale. <=0 implies
  ## fixed scale.

  qrx <- qr(X,tol=0)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot,drop=FALSE]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)
  ## remove possible ambiguity from statistic...
  siv <- sign(ed$vectors[1,]);siv[siv==0] <- 1
  ed$vectors <- sweep(ed$vectors,2,siv,"*")

  k <- max(0,floor(rank))
  nu <- abs(rank - k)     ## fractional part of supplied edf
  if (type==1) { ## round up is more than .05 above lower
    if (rank > k + .05||k==0) k <- k + 1
    nu <- 0;rank <- k
  }

  if (nu>0) k1 <- k+1 else k1 <- k

  ## check that actual rank is not below supplied rank+1
  r.est <- sum(ed$values > max(ed$values)*.Machine$double.eps^.9)
  if (r.est<k1) {k1 <- k <- r.est;nu <- 0;rank <- r.est}

  ## Get the eigenvectors...
  # vec <- qr.qy(qrx,rbind(ed$vectors,matrix(0,nrow(X)-ncol(X),ncol(X))))
  vec <- ed$vectors
  if (k1<ncol(vec)) vec <- vec[,1:k1,drop=FALSE]

  ## deal with the fractional part of the pinv...
  if (nu>0&&k>0) {
    if (k>1) vec[,1:(k-1)] <- t(t(vec[,1:(k-1)])/sqrt(ed$val[1:(k-1)]))
    b12 <- .5*nu*(1-nu)
    if (b12<0) b12 <- 0
    b12 <- sqrt(b12)
    B <- matrix(c(1,b12,b12,nu),2,2)
    ev <- diag(ed$values[k:k1]^-.5,nrow=k1-k+1)
    B <- ev%*%B%*%ev
    eb <- eigen(B,symmetric=TRUE)
    rB <- eb$vectors%*%diag(sqrt(eb$values))%*%t(eb$vectors)
    vec1 <- vec
    vec1[,k:k1] <- t(rB%*%diag(c(-1,1))%*%t(vec[,k:k1]))
    vec[,k:k1] <- t(rB%*%t(vec[,k:k1]))
  } else {
    vec1 <- vec <- if (k==0) t(t(vec)*sqrt(1/ed$val[1])) else
      t(t(vec)/sqrt(ed$val[1:k]))
    if (k==1) rank <- 1
  }
  ## there is an ambiguity in the choise of test statistic, leading to slight
  ## differences in the p-value computation depending on which of 2 alternatives
  ## is arbitrarily selected. Following allows both to be computed and p-values
  ## averaged (can't average test stat as dist then unknown)
  d <- t(vec)%*%(R%*%p)
  d <- sum(d^2)
  d1 <- t(vec1)%*%(R%*%p)
  d1 <- sum(d1^2)
  ##d <- d1 ## uncomment to avoid averaging

  rank1 <- rank ## rank for lower tail pval computation below

  ## note that for <1 edf then d is not weighted by EDF, and instead is
  ## simply refered to a chi-squared 1

  if (nu>0) { ## mixture of chi^2 ref dist
    if (k1==1) rank1 <- val <- 1 else {
      val <- rep(1,k1) ##ed$val[1:k1]
      rp <- nu+1
      val[k] <- (rp + sqrt(rp*(2-rp)))/2
      val[k1] <- (rp - val[k])
    }

    if (res.df <= 0) pval <- (mgcv:::liu2(d,val) + mgcv:::liu2(d1,val))/2 else ##  pval <- davies(d,val)$Qq else
      pval <- (mgcv:::simf(d,val,res.df) + mgcv:::simf(d1,val,res.df))/2
  } else { pval <- 2 }
  ## integer case still needs computing, also liu/pearson approx only good in
  ## upper tail. In lower tail, 2 moment approximation is better (Can check this
  ## by simply plotting the whole interesting range as a contour plot!)
  if (pval > .5) {
    if (res.df <= 0) pval <- (pchisq(d,df=rank1,lower.tail=FALSE)+pchisq(d1,df=rank1,lower.tail=FALSE))/2 else
      pval <- (pf(d/rank1,rank1,res.df,lower.tail=FALSE)+pf(d1/rank1,rank1,res.df,lower.tail=FALSE))/2
  }
  list(stat=d,pval=min(1,pval),rank=rank)
} ## end of testStat


plot_BSMethEM = function (BEM.obj, mfrow = NULL, same.range= F){
  ncovs = ncol(BEM.obj$Beta.out)
  
  if(is.null(mfrow)){
    par(mfrow = c(1, ncovs),mar = c(4,2.6,3,1))
  }else{
    par(mfrow = mfrow, mar = c(4,2.6,3,1))
  }
  covs.names = rownames(BEM.obj$reg.out)
  ll <- BEM.obj$Beta.out - 1.96* BEM.obj$SE.out ; hh <- BEM.obj$Beta.out + 1.96* BEM.obj$SE.out
  
  yylim <- matrix(NA, nrow = ncovs, ncol = 2)
  if(same.range){
    yylim = matrix(rep(c(min(ll, na.rm = T), max(hh, na.rm =T)), ncovs), ncol = 2, byrow = T)
  }else{
    yylim<-  t(sapply(1:ncovs, function(i){ yylim[i,] = c(ifelse(min(ll[,i], na.rm =T)>0, 0, min(ll[,i], na.rm =T)), max(hh[,i], na.rm =T))}))
  }
  
  for ( ii in 1:ncovs){
    
    plot(BEM.obj$uni.pos[order(BEM.obj$uni.pos)], BEM.obj$Beta.out[order(BEM.obj$uni.pos),ii],col="red", xaxt ="n",
         type="l",xlab="Genomic Position", ylab=" ", main = covs.names[ii] , lwd = 2, ylim = yylim[ii, ])
    
    if(ii !=1){
      mtext(paste("Pval = ", formatC(BEM.obj$reg.out[ii,"p-value"],
                                     digits=2, format="e")), side = 3, line = -1.2, cex = 0.8, font = 0.8)
    }
    
    axis(side = 1, at = BEM.obj$uni.pos[order(BEM.obj$uni.pos)],  labels=F, lwd=0.5, lwd.ticks = 0.5, tck=0.03)
    axis(side=1, at = seq(round(min(BEM.obj$uni.pos)), round(max(BEM.obj$uni.pos)),length.out = 10 ) , tck= -0.02)
    lines(BEM.obj$uni.pos[order(BEM.obj$uni.pos)], ll[order(BEM.obj$uni.pos), ii], lty = 2, col ="red")
    lines(BEM.obj$uni.pos[order(BEM.obj$uni.pos)], hh[order(BEM.obj$uni.pos), ii], lty = 2, col ="red")
    abline(h=0, lty=2)
  }
}



BSMethEMUpdate <- function (data, pi.ij, p0 , p1 , n.k,binom.link, method, Z, my.covar.fm, Quasi, scale,RanEff, N ){
  if( !(nrow(data)==length(pi.ij)) ) message("The row of data should be compatible with the length of initial value pi.ij")
  # The E-step
  # Calculate the "posterior" probability
  eta.1 <-  p1 *pi.ij /(p1*pi.ij + p0*(1-pi.ij)) # posterior probability given an observed methylated rates, what is the probability that the reads are truely methylated
  eta.0 <-  (1-p1) * pi.ij /((1-p1) * pi.ij + (1-p0)*(1-pi.ij))
  
  Y <- data$Y ; X <- data$X
  E.S <- Y * eta.1 + (X-Y) * eta.0
  
  #p0_est <- mean((Y*(1-eta.1))/(X-E.S))
  #p1_est <- mean((Y *eta.1)/E.S)
  
 # p0_est <- sum(Y*(1-eta.1))/(sum(X-E.S))
 # p1_est <- sum(Y *eta.1)/sum(E.S)
 # p0_est <- p0
 # p0_est <- sum(Y*(1-eta.1))/(sum(X-E.S))
  
  
  p0_est = p0
  p1_est = p1
  
  
  if(Quasi){
    gam.int.see <- suppressWarnings( mgcv::gam(as.formula( paste0("E.S/X ~", my.covar.fm)), family = quasibinomial(link=binom.link),weights=X,
                                               data = data, method=method, scale = scale))
  }else{
    gam.int.see <- suppressWarnings( mgcv::gam(as.formula( paste0("E.S/X ~", my.covar.fm)), family = binomial(link=binom.link),weights=X,
                                               data = data, method=method))
  }
  

  p_res <- residuals(gam.int.see, type ='pearson')
  d_res <- residuals(gam.int.see, type ='deviance')
  if(Quasi & scale<=0){
  my_s <- (1-2*gam.int.see$fitted.values)/(data$X*gam.int.see$fitted.values*(1-gam.int.see$fitted.values))*(E.S-data$X*gam.int.see$fitted.values)#*sqrt(data$X)
  phi_p = sum( p_res^2)/(length(data$Y) - sum(gam.int.see$edf))
  
  phi_fletcher <- phi_p/(1+mean(my_s))
  }
  
  if(scale >0) {
    phi_fletcher = scale
  }
  if(!Quasi){
    phi_fletcher = 1
  }

  
  # The values of Q function
  
  new.lambda = gam.int.see$sp
  pi.ij = gam.int.see$fitted.values
  par = gam.int.see$coefficients
  smoth.mat <-lapply(as.list(1:(ncol(Z)+1)), function(i){gam.int.see$smooth[[i]]$S[[1]] * new.lambda[i]})  # extract the penalty matrix
  
  smoth.mat[[length(smoth.mat) + 1]] <- 0  # assume the lambda for the constant of the intercept is 0 -- no penalization
  if(RanEff){
    #new.lambda[ncol(Z)+2]
    #smoth.mat[[length(smoth.mat) + 1]] <- diag(N) /(re_sd^2) # NOTE that the smoothing matrix for the random effects is not identity remember to multiple lambda (1/sigma_0^2)
    smoth.mat[[length(smoth.mat) + 1]] <- diag(N) *new.lambda[ncol(Z)+2] # !!!! Otherwise, we get very wide CI
    span.penal.matrix <- as.matrix(Matrix::bdiag( smoth.mat[c(length(smoth.mat)-1, (1:(length(smoth.mat)-2)), length(smoth.mat))] ))
  }else{
    
    span.penal.matrix <- as.matrix(Matrix::bdiag( smoth.mat[c(length(smoth.mat), (1:(length(smoth.mat)-1)))] ))
  }
  
  Q1 = sum(E.S*log(pi.ij) +(X-E.S)*log(1-pi.ij)) - 0.5 * t(par) %*% span.penal.matrix %*% par
  Qp1 = sum( log(p1_est) * Y *eta.1 + log(1-p1_est) * (X-Y)*eta.0 )
  Qp0 = sum( log(p0_est) * Y *(1-eta.1) + log(1-p0_est) * (X-Y)*(1-eta.0) )
  
 # phi_fletcher = summary(gam.int.see)$dispersion # this is actually the fixed scale paramters in the input
  out <- list( pi.ij = gam.int.see$fitted.values, par = gam.int.see$coefficients,
               lambda = gam.int.see$sp, edf1 = gam.int.see$edf1, pearson_res = p_res, 
               edf=gam.int.see$edf, phi_fletcher= phi_fletcher, GamObj = gam.int.see , E.S = E.S,
               p0_est=p0_est,
               p1_est=p1_est, Q1 = Q1, Qp1 = Qp1, Qp0 = Qp0)
  return(out)
}

