loc_est_bw <- function(xtab, ytab, x, hini, B = 5, method = "u", fix.seed = FALSE, control = list("tm_limit" = 700)){
   # verification
   stopifnot(length(xtab)==length(ytab))
  
   # initialiasation
   h_loc_min <- max(diff(sort(xtab)))
   h_loc_max <- (max(xtab)-min(xtab))/2
   h_loc_range <- seq(h_loc_min, h_loc_max, length = 30)
   BIterN <- B # number of bootstrap iteration

   ghat <- loc_est(xtab, ytab, x, hini, method = method, control = control) # a pilot estimator
   gpil <- length(ghat) # a second pilot estimator
   h1 <- 1.5*hini

   # internal function
   ckernel <- function(x, bandwidth){
    as.numeric((x/bandwidth<=0.5) & (x/bandwidth>=-0.5))*bandwidth^{-1}
   }

   for(i in 1:length(x)){
     gpil[i] <- sum(ghat*ckernel(x[i]-x, h1), na.rm = T)*unique(diff(x))[1]/(sum(ckernel(x[i]-x, h1), na.rm = T)*unique(diff(x))[1]) # smoothing
   }

   ####################################
   # STEA (a) in Hall and Park (2004) #
   ####################################

   rm_ind <- c()

   for (i in 1:length(xtab)){
     if(ytab[i]>gpil[which.min( abs(xtab[i]-x))]){
       rm_ind <- c(rm_ind, i)
     }
   }

   # L_1
   xtab_r <- xtab[-rm_ind]
   ytab_r <- ytab[-rm_ind]

   ####################################
   # STEA (b) in Hall and Park (2004) #
   ####################################

   # L
   xtab_e <- xtab_r
   ytab_e <- ytab_r

   for(i in 1:length(xtab_r)){
     ytab_e <- c(ytab_e, 2*gpil[which.min(abs(xtab_r[i]-x))]-ytab_r[i])
     xtab_e <- c(xtab_e, xtab_r[i])
   }

   ######################################
   # STEA (c,d) in Hall and Park (2004) #
   ######################################

   m <- as.matrix(dist(cbind(xtab_e, ytab_e)))
   sortm <- apply(m, 2, sort)
   k <- 10
   DZ <- sortm[k+1,]

   MSE <- rep(0, length(h_loc_range))
   count <- rep(0, length(h_loc_range))

   for(iterB in 1:BIterN){
    cat("Bootstrap Sample #", iterB, "\n")
     if(fix.seed)
       set.seed(iterB)
    NZ <- rpois(length(xtab_e), 1)
    bxtab <- c()
    bytab <- c()

    for (i in 1:length(xtab_e)){
     if(NZ[i]>0){
       if(fix.seed)
         set.seed(i)
       R <- runif(NZ[i], 0, DZ[i])
       if(fix.seed)
         set.seed(i+1)
       theta <- runif(NZ[i], 0, 2*pi)
       bxtab <- c(bxtab, xtab_e[i] + R*cos(theta))
       bytab <- c(bytab, ytab_e[i] + R*sin(theta))
      }
     }

    rm_ind_b <- c()

    for(i in 1:length(bxtab)){
      if(bytab[i]>gpil[which.min(abs(bxtab[i]-x))]){
        rm_ind_b <- c(rm_ind_b, i)
      }
    }

    bxtab_r <- bxtab[-rm_ind_b]
    bytab_r <- bytab[-rm_ind_b]

    h_chk_min <- max(diff(sort(bxtab_r)))

    if(h_chk_min < max(h_loc_range)){
      count <- count + as.numeric(h_loc_range>=h_chk_min)
      xb <- x[(x>=min(bxtab_r)) & (x<=max(bxtab_r))]
      gpilb <- gpil[(x>=min(bxtab_r)) & (x<=max(bxtab_r))]

     for(hi in min(which(h_loc_range>=h_chk_min)):length(h_loc_range)){
       est <- loc_est(bxtab_r, bytab_r, xb, h = h_loc_range[hi], method = method, control = control)
       MSE[hi] <- MSE[hi] + mean((est-gpilb)^2, na.rm = T)
     }
    }
   }

   h_loc_range_r <- h_loc_range[count>0]
   hbopt <- h_loc_range_r[which.min(MSE[count>0]*(1/count[count>0]))]
   return(hbopt)
 }