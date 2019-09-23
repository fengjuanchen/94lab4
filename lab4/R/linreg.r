linreg <- function(formula, data){
  X <- model.matrix(formula,data)
  Z <- all.vars(formula)
  y <- Z[1]
  lnrg <- list()
  class(lnrg) <- "linreg"
  
  lnrg$regcoe <- solve(t(X) %*% X) %*% t(X) %*% data[,y]   
  lnrg$fitval <- X %*% lnrg$regcoe
  lnrg$residu <- data[,y]-lnrg$fitval
  lnrg$degfre <- nrow(data)-ncol(X) 
  lnrg$resvar <- (t(lnrg$residu) %*% lnrg$residu) /lnrg$degfre
  sca_resvar <- lnrg$resvar[1,1]
  
  regcoemat1 <- solve(t(X) %*% X)
  regcoemat2 <- vector(length = nrow(regcoemat1))
  for (i in 1:nrow(regcoemat1)) {
    regcoemat2[i] <- regcoemat1[i,i]
  }
  
  lnrg$varregcoe <- regcoemat2 * sca_resvar
  
  tvalue <- vector(length = nrow(lnrg$regcoe))
  for (i in 1:nrow(lnrg$regcoe)) {
    tvalue[i] <- lnrg$regcoe[i,1]/sqrt(lnrg$varregcoe[i])
  }
  lnrg$tval <- tvalue
  return(lnrg)
}

resid.linreg <- function(x){
  #print.default(x$residu)
  return(x$residu)
}


pred.linreg <- function(x){
  return(x$fitval)
}
coef.linreg <- function(x){
  named_coe <- vector(length = length(x$regcoe))
  vector_name <- vector(length = length(x$regcoe))
  for (i in 1:nrow(x$regcoe)) {
    named_coe[i] <- x$regcoe[i,1]
    vector_name[i] <- names(x$regcoe[i,1])
  }
  names(named_coe) <- vector_name
  return(named_coe)
}


summary.linreg <- function(x){
  if(length(coef(x))){
    cat("Coefficients:\n")
    stderror <- sqrt(x$varregcoe)
    
    summ <- cbind(x$regcoe,stderror,x$tval)
    colnames(summ) <- c("Estimate", "Std. Error", "t value")
    print.default(summ, print.gap = 2,quote = FALSE)
    cat("\n\nResidual standard error:", sqrt(x$resvar),"on",x$degfre, "degrees of freedom" )
  }
  else cat("No coefficients\n")
}


print.linreg <- function(x, digits=max(3,getOption("digits")-3)){
  cat("\nCall:\n", deparse(x), "\n\n", sep = "")
  if(length(coef(x))){
    cat("Coefficients:\n")
    print.default(format(coef(x),digits = digits), print.gap = 2,quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
  
}