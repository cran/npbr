poly_est <- function(xtab, ytab, x, deg, control = list("tm_limit" = 700)){
  
 # verification
  stopifnot(length(xtab)==length(ytab))
  
 # initialisation  
  expo <- 1:(deg+1)
  rec <- 1/expo
  l_end <- min(xtab)
  r_end <- max(xtab)

  opt_coef <- rec*(r_end^expo - l_end^expo)
  desM <- c()

  for(i in 0:deg){
   desM <- cbind(desM, xtab^i)
  }

  # formulation of linear programming
  obj <- opt_coef
  mat <- desM
  rhs <- ytab
  dir <- c(rep(">=", length(xtab)))
  bounds <- list(lower = list(ind = 1:(deg+1), val = rep(-Inf, deg+1)),
                 upper = list(ind = 1:(deg+1), val = rep(Inf, deg+1)))

  Sol <- Rglpk_solve_LP(obj, mat, dir, rhs, bounds, types = NULL, max = FALSE, control = control)

  if(Sol$status!=0){
    # warning("It seems no optimal solution has been found by Glpk")
    OPT <- rep(NA, (deg+1)) 
  }else
    OPT <- Sol$solution
  
  gridM <- c()
 
  for(i in 0:deg){
   gridM <- cbind(gridM, x^i)
   }

  fitt <- gridM %*% OPT
  
  return(fitt)
 }