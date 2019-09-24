
linreg <- function(formula, data){
  #if() stop("wrong dataset name")
  Z <- all.vars(formula)
  for(i in 1:length(Z)){
    if(!(Z[i] %in% names(data))) stop("wrong column name!")
  }
  
  X <- model.matrix(formula,data)
  
  y <- Z[1]
  
  cl <- match.call()
  
  lnrg <- list()
  class(lnrg) <- "linreg"
  
  new_qr_x <- qr(X)
  lnrg$regcoe <- solve(qr.R(new_qr_x)) %*% t(qr.Q(new_qr_x)) %*% data[,y]   
  lnrg$fitval <- qr.fitted(new_qr_x,data[,y])
  lnrg$residu <- qr.resid(new_qr_x,data[,y])
  lnrg$degfre <- nrow(data)-ncol(X) 
  lnrg$resvar <- (t(lnrg$residu) %*% lnrg$residu) /lnrg$degfre
  sca_resvar <- lnrg$resvar[1,1]
  
  regcoemat1 <- chol2inv(qr.R(new_qr_x))
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
  lnrg$pval <- pt(abs(lnrg$tval),lnrg$degfre, lower.tail = FALSE)
  lnrg$call <- cl
  
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
    
    summ <- cbind(x$regcoe,stderror,x$tval, x$pval)
    colnames(summ) <- c("Estimate", "Std. Error", "t value","Pr(>|t|)")
    print.default(summ, print.gap = 2,quote = FALSE)
    cat("\n\nResidual standard error:", sqrt(x$resvar),"on",x$degfre, "degrees of freedom" )
  }
  else cat("No coefficients\n")
}

plot.linreg <- function(x){
  
  plot(x$fitval,x$residu)
  sca_resvar <- x$resvar[1,1]
  st <- sqrt(sca_resvar)
  abs_residu <- abs(x$residu)/st
  plot(x$fitval,sqrt(abs_residu))
  
  
}

print.linreg <- function(x, digits=max(3,getOption("digits")-3)){
  cat("\nCall:\n", deparse(x$call),"\n\n", sep = "")
  if(length(coef(x))){
    cat("Coefficients:\n")
    print.default(format(coef(x),digits = digits), print.gap = 2,quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
  
}



