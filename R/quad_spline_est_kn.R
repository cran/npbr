 quad_spline_est_kn<-function(xtab,ytab,x,cv,krange=1:20,type="AIC")
 {
  stopifnot(type %in% c("AIC","BIC"),length(xtab)==length(ytab), cv%in%c(0,1), is.integer(krange))
       
  criteria<-NULL
 
  if(type=="AIC")
  {for (k in krange)
   {criteria<-c(criteria,log(sum(abs(ytab-quad_spline_est(xtab,ytab,xtab,k,cv))))+2*(k+3)/length(xtab))}
  } 
  else
  {for (k in krange)
   {criteria<-c(criteria,log(sum(abs(ytab-quad_spline_est(xtab,ytab,xtab,k,cv))))+log(length(xtab))*(k+3)/length(xtab))}
  }
  
 return(krange[which.min(criteria)])
 }