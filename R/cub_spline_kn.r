cub_spline_kn <- function(xtab, ytab, method, krange = 1:20, type = "AIC", control = list("tm_limit" = 700)){
 # verification
  stopifnot(type %in% c("AIC","BIC"), length(xtab)==length(ytab), method%in%c("u","m","mc"), is.integer(krange))
 
  # initialisation       
  n.k <- length(krange)
  criteria <- numeric(n.k)
  n <- length(xtab)
  
  if(type=="AIC"){
    for(k in 1:n.k){
      est <- cub_spline_est(xtab, ytab, xtab, kn = krange[k], method = method, control = control)
      criteria[k] <- log(sum(abs(ytab - est))) + (krange[k]+3)/n}
    }else
    {for(k in 1:n.k){
      est = cub_spline_est(xtab, ytab, xtab, kn = krange[k], method = method, control = control)
      criteria[k] <- log(sum(abs(ytab - est))) + log(n)*(krange[k]+3)/(2*n)
      }
     }
  return(krange[which.min(criteria)])
}
