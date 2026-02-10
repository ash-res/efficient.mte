#' Efficient Semiparametric Estimation of Marginal Treatment Effects
#'
#' @description A method for heterogeneous treatment effects analysis in a semiparametric MTE model
#'
#' @param A a vector of a binary treatment variable taking values 0 and 1
#' @param X a vector (if there is only a single covariate) or a matrix (if multiple covariates) of covariate(s), where the k-th column gives the k-th covariate
#' @param Y a vector of an outcome variable
#' @param Z a vector of an instrumental variable
#' @param ps a vector of propensity scores obtained from fitting either a parametric or nonparametric regression of A on (X,Z)
#' @param S (optional) a positive integer indicating the order of polynomial of MTE model (default is S=1, i.e. a linear model)
#' @param alpha (optional) (1-alpha)*100% is the nominal coverage of confidence intervals reported in the table of results; default is 0.05
#' @param interaction (optional) a logical; if \code{TRUE} then an MTE model allowing for X and V interaction effects is estimated. The default is \code{FALSE}
#' @param trim (optional) a non-negative constant between 0 and 1 indicating whether observations with extreme propensity scores less than \code{1-trim} and greater than \code{trim} should be excluded. The default is \code{trim=0} so that no observations are excluded
#' @param pc.pval (optional) p-value threshold to select PC-based instruments for each exposure. Note that a very strict p-value threshold may leave very few instruments (perhaps none), so the choice should be made with the sample size in mind. The default choice is 5e-8
#' @param clump.thres (optional) R^2 clumping threshold for variant correlations; if unspecified the default value is R^2 = 0.95
#' @param exposure.names (optional) a J x 1 character vector of the names of exposures included in the model
#'
#' @details Efficient semiparametric MTE estimation - for details of the method please see the paper "Efficient semiparametric estimation of marginal treatment effects with genetic instrumental variables"
#'
#' @return Output is a list containing:
#' \itemize{
#'  \item \code{gamma_lr} \cr
#'  a vector of efficient MTE parameter \code{gamma} estimates
#'  \item \code{gamma_ols} \cr
#'  a vector of conventional (non-robust) MTE parameter \code{gamma} estimates by a simple linear regression
#'  \item \code{se_gamma_lr} \cr
#'  a vector of standard errors relating to \code{gamma_lr}
#'  \item \code{var_gamma} \cr
#'  an estimate of the asymptotic variance-covariance matrix relating to \code{gamma_lr}
#'  \item \code{se_gamma_ols} \cr
#'  a vector of standard errors relating to \code{gamma_ols}
#'  \item \code{ate_eff} \cr
#'  an efficient estimate of the average treatment effect estimand
#'  \item \code{ate_est} \cr
#'  a conventional (non-efficient, non-robust) estimate of the average treatment effect estimand
#'  \item \code{ate_se} \cr
#'  a standard error relating to \code{ate_eff}
#'  \item \code{att_eff} \cr
#'  an efficient estimate of the average treatment effect on the treated estimand
#'  \item \code{att_est} \cr
#'  a conventional (non-efficient, non-robust) estimate of the average treatment effect on the treated estimand
#'  \item \code{att_se} \cr
#'  a standard error relating to \code{att_eff}
#'  \item \code{atu_eff} \cr
#'  an efficient estimate of the average treatment effect on the untreated estimand
#'  \item \code{atu_est} \cr
#'  a conventional (non-efficient, non-robust) estimate of the average treatment effect on the untreated estimand
#'  \item \code{atu_se} \cr
#'  a standard error relating to \code{atu_eff}
#'  \item \code{asg_eff} \cr
#'  an efficient estimate of the average selection on the gain estimand
#'  \item \code{asg_est} \cr
#'  a conventional (non-efficient, non-robust) estimate of the average selection on the gain estimand
#'  \item \code{asg_se} \cr
#'  a standard error relating to \code{asg_eff}
#'  \item \code{plot.mte} \cr
#'  a plot of the efficient estimate of the MTE function averaged over covariates (reported only if there is a single covariate)
#' }
#'
#'
#' @author Ashish Patel, Francis J. DiTraglia, and Stephen Burgess
#'
#' @export

eff.mte <- function(A,X,Y,Z,ps,S=1,alpha=0.05,interaction=FALSE,trim=0){
  if(interaction==FALSE){
    X <- as.matrix(X)
    sel.rm <- which(ps<trim|ps>(1-trim))
    if(length(sel.rm)>0){Y=Y[-sel.rm];A=A[-sel.rm];X=X[-sel.rm,];Z=Z[-sel.rm];ps=ps[-sel.rm]}
    
    # first, lets estimate parameters based on the local IV approach regression equation
    n <- length(A)
    tS <- matrix(NA,nrow=n,ncol=S); tS1 <- matrix(NA,nrow=n,ncol=S)
    for (s in 1:S){tS[,s] <- ps^(s+1); tS1[,s] <- (ps^s)*(s+1)}
    r <- cbind(1,X,ps,X*ps,tS)
    R <- cbind(0,matrix(0,ncol=ncol(X),nrow=nrow(X)),rep(1,length=length(ps)),X,tS1)
    g <- function(gamma){r*(Y-as.vector(r%*%gamma))}
    phi <- function(gamma){-r*as.vector((R*(A-ps))%*%gamma)}
    psi <- function(gamma){g(gamma)+phi(gamma)}
    Omega <- (t(r)%*%r)/n
    Gamma <- t(r)%*%(R*(A-ps))/n
    varg <- function(gamma){t(g(gamma))%*%g(gamma)/n}
    varpsi <- function(gamma){t(psi(gamma))%*%psi(gamma)/n}
    Qg <- function(gamma){as.numeric(t(colMeans(g(gamma)))%*%solve(varg(gamma))%*%colMeans(g(gamma)))}
    DQg <- function(gamma){-as.vector(2*t(Omega)%*%solve(varg(gamma))%*%colMeans(g(gamma)))}
    # gamma_ols <- nlminb(rep(0,ncol(r)),objective=Qg,gradient=DQg)$par is the same as below
    ols.reg <- lm(Y~r-1)
    gamma_ols <- as.vector(ols.reg$coefficients)
    se_gamma_ols <- sqrt(diag((solve(t(Omega))%*%varpsi(gamma_ols)%*%solve(Omega))/n))
    se_gamma_ols_naive <- sqrt(diag((solve(t(Omega))%*%varg(gamma_ols)%*%solve(Omega))/n))
    gamma_lr <- as.vector(solve(Omega+Gamma)%*%(t(r)%*%Y/n))
    se_gamma_lr <- sqrt(diag((solve(t(Omega+Gamma))%*%varpsi(gamma_lr)%*%solve(Omega+Gamma))/n))
    var_gamma <- ((solve(t(Omega+Gamma))%*%varpsi(gamma_lr)%*%solve(Omega+Gamma))/n)
    
    # ATE estimand and point estimates
    r_ate <- cbind(0, matrix(0,ncol=ncol(X),nrow=nrow(X)), 1, X, matrix(1,nrow=n,ncol=S))
    ate_est <- as.numeric(t(colMeans(r_ate))%*%gamma_ols)
    ate_eff <- as.numeric(t(colMeans(r_ate))%*%gamma_ols)+as.numeric(t(colMeans(r_ate))%*%solve(Omega)%*%colMeans(psi(gamma_ols)))
    phi_ate <- as.vector(r_ate%*%gamma_ols) - ate_eff + as.vector(psi(gamma_ols)%*%(solve(Omega)%*%colMeans(r_ate)))
    avar_ate <- (t(phi_ate)%*%phi_ate/n)
    ate_se <- as.numeric(sqrt(avar_ate/n)) # see Hines et al.: 1/n times the sample variance of influence function is the asymp var of estimand estimator 
    
    # ATT estimand and point estimates
    r_att <- cbind(0, matrix(0,ncol=ncol(X),nrow=nrow(X)), ps, X*ps, tS)/mean(A==1)
    R_att <- cbind(0,matrix(0,ncol=ncol(X),nrow=nrow(X)),1,X,tS1)/mean(A==1)
    att_est <- as.numeric(t(colMeans(r_att))%*%gamma_ols)
    att_eff <- (as.numeric(t(colMeans(r_att))%*%gamma_ols))+(as.numeric(t(colMeans(R_att*(A-ps)))%*%gamma_ols))+(as.numeric(t(colMeans(r_att))%*%solve(Omega)%*%colMeans(psi(gamma_ols))))
    phi_att <- -(A*att_eff/mean(A==1))+(as.vector(r_att%*%gamma_ols))+(as.vector((R_att*(A-ps))%*%gamma_ols))+(as.vector(psi(gamma_ols)%*%(solve(Omega)%*%colMeans(r_att))))
    avar_att <- (t(phi_att)%*%phi_att/n)
    att_se <- as.numeric(sqrt(avar_att/n))
    
    # ATU estimand and point estimates
    r_atu <- cbind(0, matrix(0,ncol=ncol(X),nrow=nrow(X)), 1-ps, X*(1-ps),1-tS)/mean(A==0)
    R_atu <- cbind(0,matrix(0,ncol=ncol(X),nrow=nrow(X)),-1,-X,-tS1)/mean(A==0)
    atu_est <- as.numeric(t(colMeans(r_atu))%*%gamma_ols)
    atu_eff <- (as.numeric(t(colMeans(r_atu))%*%gamma_ols))+(as.numeric(t(colMeans(R_atu*(A-ps)))%*%gamma_ols))+(as.numeric(t(colMeans(r_atu))%*%solve(Omega)%*%colMeans(psi(gamma_ols))))
    phi_atu <- ((A-1)*atu_eff/mean(A==0))+(as.vector(r_atu%*%gamma_ols))+(as.vector((R_atu*(A-ps))%*%gamma_ols))+(as.vector(psi(gamma_ols)%*%(solve(Omega)%*%colMeans(r_atu))))
    avar_atu <- (t(phi_atu)%*%phi_atu/n)
    atu_se <- as.numeric(sqrt(avar_atu/n))
    
    # ASG estimand and point estimates
    asg_est <- att_est-atu_est
    asg_eff <- att_eff-atu_eff
    avar_asg <- (t(phi_att-phi_atu)%*%(phi_att-phi_atu)/n)
    asg_se <- as.numeric(sqrt(avar_asg/n))
    
    # Show table
    decimals <- function(number, places) format(round(number, places), nsmall = places)
    del_par <- c("ATE","ATT","ATU","ASG")
    del_est <- c(ate_eff,att_eff,atu_eff,asg_eff)
    del_se <- c(ate_se,att_se,atu_se,asg_se)
    del_lower <- c(ate_eff-qnorm(1-alpha/2)*ate_se,att_eff-qnorm(1-alpha/2)*att_se,atu_eff-qnorm(1-alpha/2)*atu_se,asg_eff-qnorm(1-alpha/2)*asg_se)
    del_upper <- c(ate_eff+qnorm(1-alpha/2)*ate_se,att_eff+qnorm(1-alpha/2)*att_se,atu_eff+qnorm(1-alpha/2)*atu_se,asg_eff+qnorm(1-alpha/2)*asg_se)
    res.table <- data.frame(del_par,decimals(del_est,3),decimals(del_se,3),decimals(del_lower,3),decimals(del_upper,3))
    colnames(res.table) <- c("param","estimate","std err","95% CI lower","95% CI upper")
    cat("\nMTE model: X and V have separable effects\n")
    cat("\n---------------Efficient Target Parameter Estimates---------------\n")
    print(res.table,row.names=F)
    cat("------------------------------------------------------------------\n")
    
    if(ncol(X)>1){
      return(list("gamma_lr"=gamma_lr,"gamma_ols"=gamma_ols,"se_gamma_lr"=se_gamma_lr,"se_gamma_ols"=se_gamma_ols,"ate_est"=ate_est,"ate_eff"=ate_eff,"ate_se"=ate_se,"att_est"=att_est,"att_eff"=att_eff,"att_se"=att_se,"atu_est"=atu_est,"atu_eff"=atu_eff,"atu_se"=atu_se,"asg_est"=asg_est,"asg_eff"=asg_eff,"asg_se"=asg_se,"var_gamma"=var_gamma))
    }
    
    # MTE plot if X is a single binary variable
    if(ncol(X)==1){
      mte <- function(v0){
        mteS <- vector(,length=1)
        for (s in 1:S){mteS[s] <- (v0^s)*(s+1)}
        r_mte <- function(x0){c(0,0,1,x0,mteS)}
        suppX <- unique(X[,1])
        mte_eff <- sum(sapply(suppX,function(x){mean(X[,1]==x)*as.numeric(t(r_mte(x))%*%gamma_lr)}))
        mte_var <- sum(sapply(suppX,function(x){mean(X[,1]==x)*as.numeric(t(r_mte(x))%*%var_gamma%*%r_mte(x))}))
        mte_lower <- mte_eff - qnorm(1-alpha/2)*sqrt(mte_var)
        mte_upper <- mte_eff + qnorm(1-alpha/2)*sqrt(mte_var)
        return(list("est"=mean(mte_eff),"lower"=mean(mte_lower),"upper"=mean(mte_upper)))
      }
      mte.sq <- seq(0,1,0.1)
      mte <- sapply(mte.sq,mte)
      mte.est <- function(v){mte[,v]$est}; mte.est <- sapply(1:length(mte.sq),mte.est)
      mte.lower <- function(v){mte[,v]$lower}; mte.lower <- sapply(1:length(mte.sq),mte.lower)
      mte.upper <- function(v){mte[,v]$upper}; mte.upper <- sapply(1:length(mte.sq),mte.upper)
      
      V.breaks <- seq(0,1,0.2)
      dat.mte <- data.frame(V=rep(mte.sq,3), method=factor(rep(c("point",paste0((1-alpha)*100,"% CI (lower)"),paste0((1-alpha)*100,"% CI (upper)")), each=(length(mte.sq)))), est=c(mte.est,mte.lower,mte.upper))
      plot.mte <- ggplot2::ggplot(data=dat.mte, ggplot2::aes(x=V, y=est, group=method)) + 
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), colour="grey60", linetype = "dotted") + 
        ggplot2::theme_bw() + 
        ggplot2::scale_x_continuous(breaks=V.breaks) +  
        ggplot2::geom_line(ggplot2::aes(color = method, linetype = method)) + 
        ggplot2::scale_linetype_manual(values=c(2,2,1)) + 
        ggplot2::scale_color_manual(values=c('grey20','grey20','dodgerblue')) + 
        ggplot2::xlab("unobserved resistance to take treatment: V") + 
        ggplot2::ylab("MTE(V)") + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                       plot.subtitle = ggplot2::element_text(hjust = 0.5)) + 
        ggplot2::theme(panel.border = ggplot2::element_blank(), 
                       panel.grid.major = ggplot2::element_blank(), 
                       panel.grid.minor = ggplot2::element_blank(), 
                       axis.line = ggplot2::element_line(colour = "black")) + 
        ggplot2::theme(axis.text = ggplot2::element_text(size=12), 
                       legend.text = ggplot2::element_text(size=12), 
                       legend.title = ggplot2::element_text(size=14))
    return(list("gamma_lr"=gamma_lr,"gamma_ols"=gamma_ols,"se_gamma_lr"=se_gamma_lr,"se_gamma_ols"=se_gamma_ols,"ate_est"=ate_est,"ate_eff"=ate_eff,"ate_se"=ate_se,"att_est"=att_est,"att_eff"=att_eff,"att_se"=att_se,"atu_est"=atu_est,"atu_eff"=atu_eff,"atu_se"=atu_se,"asg_est"=asg_est,"asg_eff"=asg_eff,"asg_se"=asg_se,"var_gamma"=var_gamma,"plot"=plot.mte))
    }
  }
  if(interaction==TRUE){
    X <- as.matrix(X)
    sel.rm <- which(ps<trim|ps>(1-trim))
    if(length(sel.rm)>0){Y=Y[-sel.rm];A=A[-sel.rm];X=X[-sel.rm,];Z=Z[-sel.rm];ps=ps[-sel.rm]}
    
    # first, lets estimate parameters based on the local IV approach regression equation
    n <- length(A)
    tS <- matrix(NA,nrow=n,ncol=S); tS1 <- matrix(NA,nrow=n,ncol=S); tSX <- list(); tSX1 <- list()
    for (s in 1:S){tS[,s] <- ps^(s+1); tS1[,s] <- (ps^s)*(s+1)}
    for (k in 1:ncol(X)){tSX[[k]] <- tS*X[,k]; tSX1[[k]] <- tS1*X[,k]}
    tSX <- do.call(cbind,tSX); tSX1 <- do.call(cbind,tSX1)
    r <- cbind(1,X,ps,X*ps,tS,tSX)
    R <- cbind(0,matrix(0,ncol=ncol(X),nrow=nrow(X)),rep(1,length=length(ps)),X,tS1,tSX1)
    g <- function(gamma){r*(Y-as.vector(r%*%gamma))}
    phi <- function(gamma){-r*as.vector((R*(A-ps))%*%gamma)}
    psi <- function(gamma){g(gamma)+phi(gamma)}
    Omega <- (t(r)%*%r)/n
    Gamma <- t(r)%*%(R*(A-ps))/n
    varg <- function(gamma){t(g(gamma))%*%g(gamma)/n}
    varpsi <- function(gamma){t(psi(gamma))%*%psi(gamma)/n}
    Qg <- function(gamma){as.numeric(t(colMeans(g(gamma)))%*%solve(varg(gamma))%*%colMeans(g(gamma)))}
    DQg <- function(gamma){-as.vector(2*t(Omega)%*%solve(varg(gamma))%*%colMeans(g(gamma)))}
    # gamma_ols <- nlminb(rep(0,ncol(r)),objective=Qg,gradient=DQg)$par is the same as below
    ols.reg <- lm(Y~r-1)
    gamma_ols <- as.vector(ols.reg$coefficients)
    se_gamma_ols <- sqrt(diag((solve(t(Omega))%*%varpsi(gamma_ols)%*%solve(Omega))/n))
    se_gamma_ols_naive <- sqrt(diag((solve(t(Omega))%*%varg(gamma_ols)%*%solve(Omega))/n))
    gamma_lr <- as.vector(solve(Omega+Gamma)%*%(t(r)%*%Y/n))
    se_gamma_lr <- sqrt(diag((solve(t(Omega+Gamma))%*%varpsi(gamma_lr)%*%solve(Omega+Gamma))/n))
    var_gamma <- ((solve(t(Omega+Gamma))%*%varpsi(gamma_lr)%*%solve(Omega+Gamma))/n)
    
    # ATE estimand and point estimates
    r_ate <- cbind(0, matrix(0,ncol=ncol(X),nrow=nrow(X)), 1, X, matrix(1,nrow=n,ncol=S),matrix(rep(X,S),ncol=(ncol(X)*S)))
    ate_est <- as.numeric(t(colMeans(r_ate))%*%gamma_ols)
    ate_eff <- as.numeric(t(colMeans(r_ate))%*%gamma_ols)+as.numeric(t(colMeans(r_ate))%*%solve(Omega)%*%colMeans(psi(gamma_ols)))
    phi_ate <- as.vector(r_ate%*%gamma_ols) - ate_eff + as.vector(psi(gamma_ols)%*%(solve(Omega)%*%colMeans(r_ate)))
    avar_ate <- (t(phi_ate)%*%phi_ate/n)
    ate_se <- as.numeric(sqrt(avar_ate/n)) # see Hines et al.: 1/n times the sample variance of influence function is the asymp var of estimand estimator 
    
    # ATT estimand and point estimates
    r_att <- cbind(0, matrix(0,ncol=ncol(X),nrow=nrow(X)), ps, X*ps, tS,tSX)/mean(A==1)
    R_att <- cbind(0,matrix(0,ncol=ncol(X),nrow=nrow(X)),1,X,tS1,tSX1)/mean(A==1)
    att_est <- as.numeric(t(colMeans(r_att))%*%gamma_ols)
    att_eff <- (as.numeric(t(colMeans(r_att))%*%gamma_ols))+(as.numeric(t(colMeans(R_att*(A-ps)))%*%gamma_ols))+(as.numeric(t(colMeans(r_att))%*%solve(Omega)%*%colMeans(psi(gamma_ols))))
    phi_att <- -(A*att_eff/mean(A==1))+(as.vector(r_att%*%gamma_ols))+(as.vector((R_att*(A-ps))%*%gamma_ols))+(as.vector(psi(gamma_ols)%*%(solve(Omega)%*%colMeans(r_att))))
    avar_att <- (t(phi_att)%*%phi_att/n)
    att_se <- as.numeric(sqrt(avar_att/n))
    
    # ATU estimand and point estimates
    tSXm <- list() 
    for (k in 1:ncol(X)){tSXm[[k]] <- (1-tS)*X[,k]}
    tSXm <- do.call(cbind,tSXm)
    r_atu <- cbind(0, matrix(0,ncol=ncol(X),nrow=nrow(X)), 1-ps, X*(1-ps),1-tS,tSXm)/mean(A==0)
    R_atu <- cbind(0,matrix(0,ncol=ncol(X),nrow=nrow(X)),-1,-X,-tS1,-tSX1)/mean(A==0)
    atu_est <- as.numeric(t(colMeans(r_atu))%*%gamma_ols)
    atu_eff <- (as.numeric(t(colMeans(r_atu))%*%gamma_ols))+(as.numeric(t(colMeans(R_atu*(A-ps)))%*%gamma_ols))+(as.numeric(t(colMeans(r_atu))%*%solve(Omega)%*%colMeans(psi(gamma_ols))))
    phi_atu <- ((A-1)*atu_eff/mean(A==0))+(as.vector(r_atu%*%gamma_ols))+(as.vector((R_atu*(A-ps))%*%gamma_ols))+(as.vector(psi(gamma_ols)%*%(solve(Omega)%*%colMeans(r_atu))))
    avar_atu <- (t(phi_atu)%*%phi_atu/n)
    atu_se <- as.numeric(sqrt(avar_atu/n))
    
    # ASG estimand and point estimates
    asg_est <- att_est-atu_est
    asg_eff <- att_eff-atu_eff
    avar_asg <- (t(phi_att-phi_atu)%*%(phi_att-phi_atu)/n)
    asg_se <- as.numeric(sqrt(avar_asg/n))
    
    # Show table
    decimals <- function(number, places) format(round(number, places), nsmall = places)
    del_par <- c("ATE","ATT","ATU","ASG")
    del_est <- c(ate_eff,att_eff,atu_eff,asg_eff)
    del_se <- c(ate_se,att_se,atu_se,asg_se)
    del_lower <- c(ate_eff-qnorm(1-alpha/2)*ate_se,att_eff-qnorm(1-alpha/2)*att_se,atu_eff-qnorm(1-alpha/2)*atu_se,asg_eff-qnorm(1-alpha/2)*asg_se)
    del_upper <- c(ate_eff+qnorm(1-alpha/2)*ate_se,att_eff+qnorm(1-alpha/2)*att_se,atu_eff+qnorm(1-alpha/2)*atu_se,asg_eff+qnorm(1-alpha/2)*asg_se)
    res.table <- data.frame(del_par,decimals(del_est,3),decimals(del_se,3),decimals(del_lower,3),decimals(del_upper,3))
    colnames(res.table) <- c("param","estimate","std err","95% CI lower","95% CI upper")
    cat("\nMTE model: X and V have interaction effects\n")
    cat("\n---------------Efficient Target Parameter Estimates---------------\n")
    print(res.table,row.names=F)
    cat("------------------------------------------------------------------\n")
    
    if(ncol(X)>1){
      return(list("gamma_lr"=gamma_lr,"gamma_ols"=gamma_ols,"se_gamma_lr"=se_gamma_lr,"se_gamma_ols"=se_gamma_ols,"ate_est"=ate_est,"ate_eff"=ate_eff,"ate_se"=ate_se,"att_est"=att_est,"att_eff"=att_eff,"att_se"=att_se,"atu_est"=atu_est,"atu_eff"=atu_eff,"atu_se"=atu_se,"asg_est"=asg_est,"asg_eff"=asg_eff,"asg_se"=asg_se,"var_gamma"=var_gamma))
    }
    
    # MTE plot if X is a single binary variable
    if(ncol(X)==1){
      mte <- function(v0){
        mteS <- vector(,length=1)
        for (s in 1:S){mteS[s] <- (v0^s)*(s+1)}
        r_mte <- function(x0){c(0,0,1,x0,mteS,mteS*x0)}
        suppX <- unique(X[,1])
        mte_eff <- sum(sapply(suppX,function(x){mean(X[,1]==x)*as.numeric(t(r_mte(x))%*%gamma_lr)}))
        mte_var <- sum(sapply(suppX,function(x){mean(X[,1]==x)*as.numeric(t(r_mte(x))%*%var_gamma%*%r_mte(x))}))
        mte_lower <- mte_eff - qnorm(1-alpha/2)*sqrt(mte_var)
        mte_upper <- mte_eff + qnorm(1-alpha/2)*sqrt(mte_var)
        return(list("est"=mean(mte_eff),"lower"=mean(mte_lower),"upper"=mean(mte_upper)))
      }
      mte.sq <- seq(0,1,0.1)
      mte <- sapply(mte.sq,mte)
      mte.est <- function(v){mte[,v]$est}; mte.est <- sapply(1:length(mte.sq),mte.est)
      mte.lower <- function(v){mte[,v]$lower}; mte.lower <- sapply(1:length(mte.sq),mte.lower)
      mte.upper <- function(v){mte[,v]$upper}; mte.upper <- sapply(1:length(mte.sq),mte.upper)
      
      V.breaks <- seq(0,1,0.2)
      dat.mte <- data.frame(V=rep(mte.sq,3), method=factor(rep(c("point",paste0((1-alpha)*100,"% CI (lower)"),paste0((1-alpha)*100,"% CI (upper)")), each=(length(mte.sq)))), est=c(mte.est,mte.lower,mte.upper))
      plot.mte <- ggplot2::ggplot(data=dat.mte, ggplot2::aes(x=V, y=est, group=method)) + 
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0), colour="grey60", linetype = "dotted") + 
        ggplot2::theme_bw() + 
        ggplot2::scale_x_continuous(breaks=V.breaks) +  
        ggplot2::geom_line(ggplot2::aes(color = method, linetype = method)) + 
        ggplot2::scale_linetype_manual(values=c(2,2,1)) + 
        ggplot2::scale_color_manual(values=c('grey20','grey20','dodgerblue')) + 
        ggplot2::xlab("unobserved resistance to take treatment: V") + 
        ggplot2::ylab("MTE(V)") + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                       plot.subtitle = ggplot2::element_text(hjust = 0.5)) + 
        ggplot2::theme(panel.border = ggplot2::element_blank(), 
                       panel.grid.major = ggplot2::element_blank(), 
                       panel.grid.minor = ggplot2::element_blank(), 
                       axis.line = ggplot2::element_line(colour = "black")) + 
        ggplot2::theme(axis.text = ggplot2::element_text(size=12), 
                       legend.text = ggplot2::element_text(size=12), 
                       legend.title = ggplot2::element_text(size=14))
      return(list("gamma_lr"=gamma_lr,"gamma_ols"=gamma_ols,"se_gamma_lr"=se_gamma_lr,"se_gamma_ols"=se_gamma_ols,"ate_est"=ate_est,"ate_eff"=ate_eff,"ate_se"=ate_se,"att_est"=att_est,"att_eff"=att_eff,"att_se"=att_se,"atu_est"=atu_est,"atu_eff"=atu_eff,"atu_se"=atu_se,"asg_est"=asg_est,"asg_eff"=asg_eff,"asg_se"=asg_se,"var_gamma"=var_gamma,"plot"=plot.mte))
    }
  }
}
