poly_degree <- function(xtab, ytab, prange = 0:20, type = "AIC", control = list("tm_limit" = 700)){
  # verification
  stopifnot(type %in% c("AIC","BIC"), length(xtab)==length(ytab), is.integer(prange))
  
  # initialisation     
  n <- length(xtab)
  n.deg <- length(prange)
  criteria <- numeric(n.deg)
 
  if(type=="AIC"){
    for(k in 1:n.deg){
      deg <- prange[k]
      est <- poly_est(xtab, ytab, xtab, deg, control = control)
      criteria[k] <- log(sum(abs(ytab - est))) + (deg+1)/n}
    }else
    {for(k in 1:n.deg){
      deg <- prange[k]
      est <- poly_est(xtab, ytab, xtab, deg, control = control)
      criteria[k] <- log(sum(abs(ytab - est))) + log(n)*(deg+1)/(2*n)}
    }
  
 return(prange[which.min(criteria)])
 }