### R code from vignette source 'Z:/Thibault Pro/research/abdel/polynomial splines/npbr/inst/doc/ex-npbr.rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: ex-npbr.rnw:33-38
###################################################
owidth <- getOption("width")
options("width"=70)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0


###################################################
### code chunk number 2: bfig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 9, height = 7, pointsize = 14, bg = "white")


###################################################
### code chunk number 3: bfig2 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=name.file, width = 12, height = 7, pointsize = 14, bg = "white")


###################################################
### code chunk number 4: zfig (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 5: zfig2 (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 6: ex-npbr.rnw:100-103
###################################################
require("npbr")
data("green")
data("nuclear")


###################################################
### code chunk number 7: ex-npbr.rnw:106-109 (eval = FALSE)
###################################################
## plot(log(OUTPUT)~log(COST), data=green, pch=16,col='blue2')
## plot(ytab~xtab, data=nuclear, pch=16,col='blue2',
##  xlab="temperature of the reactor vessel", ylab="fracture toughness")


###################################################
### code chunk number 8: ex-npbr.rnw:114-121
###################################################
.PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=name.file, width = 12, height = 7, pointsize = 14, bg = "white")
op<-par(mfrow=c(1,2),oma=rep(0,4),cex.lab=1.2)
plot(log(OUTPUT)~log(COST), data=green, pch=16,col='blue2')
plot(ytab~xtab, data=nuclear, pch=16,col='blue2',
 xlab="temperature of the reactor vessel", ylab="fracture toughness")
par(op)
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 9: ex-npbr.rnw:153-155
###################################################
x.green <- seq(min(log(green$COST)), max(log(green$COST)), 
 length.out=1001)


###################################################
### code chunk number 10: ex-npbr.rnw:158-164
###################################################
y.dea<-dea_est(log(green$COST), log(green$OUTPUT), 
 x.green, type="dea")
y.fdh<-dea_est(log(green$COST), log(green$OUTPUT), 
 x.green, type="fdh")
y.lfdh=dea_est(log(green$COST), log(green$OUTPUT), 
 x.green, type="lfdh")


###################################################
### code chunk number 11: ex-npbr.rnw:168-174 (eval = FALSE)
###################################################
## plot(log(OUTPUT)~log(COST), data=green)
## lines(x.green, y.dea, lty=1, lwd=4, col="red")  
## lines(x.green, y.fdh, lty=2, lwd=4, col="blue")
## lines(x.green, y.lfdh, lty=3, lwd=4, col="green")   
## legend("topleft", legend=c("dea","fdh","lfdh"), 
##  col=c("red","blue","green"), lty=1:3, lwd=4)


###################################################
### code chunk number 12: ex-npbr.rnw:178-188
###################################################
.PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=name.file, width = 9, height = 7, pointsize = 14, bg = "white")
par<-c(oma=rep(0,4),cex.lab=1.2)
plot(log(OUTPUT)~log(COST), data=green)
lines(x.green, y.dea, lty=1, lwd=4, col="red")  
lines(x.green, y.fdh, lty=2, lwd=4, col="blue")
lines(x.green, y.lfdh, lty=3, lwd=4, col="green")   
legend("topleft", legend=c("dea","fdh","lfdh"), 
 col=c("red","blue","green"), lty=1:3, lwd=4)
par(op) 
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 13: ex-npbr.rnw:209-211
###################################################
x.nucl <- seq(min(nuclear$xtab), max(nuclear$xtab), 
 length.out=1001)


###################################################
### code chunk number 14: ex-npbr.rnw:216-218
###################################################
y.loc.opt<-loc_est(nuclear$xtab, nuclear$ytab, x.nucl, h=79.11877)
y.loc<-loc_est(nuclear$xtab, nuclear$ytab, x.nucl, h=40) 


###################################################
### code chunk number 15: ex-npbr.rnw:222-227 (eval = FALSE)
###################################################
## plot(ytab~xtab, data=nuclear)
## lines(x.nucl, y.loc.opt, lty=1, lwd=4, col="red")
## lines(x.nucl, y.loc, lty=2, lwd=4, col="blue") 
## legend("topleft",legend=c("h=79.11877", "h=40"), 
##  col=c("red","blue"), lwd=4, lty=c(1,2))


###################################################
### code chunk number 16: ex-npbr.rnw:231-238
###################################################
.PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=name.file, width = 9, height = 7, pointsize = 14, bg = "white")
plot(ytab~xtab, data=nuclear)
lines(x.nucl, y.loc.opt, lty=1, lwd=4, col="red")
lines(x.nucl, y.loc, lty=2, lwd=4, col="blue") 
legend("topleft",legend=c("h=79.11877", "h=40"), 
 col=c("red","blue"), lwd=4, lty=c(1,2))
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 17: ex-npbr.rnw:260-262
###################################################
y.poly.2<-poly_est(nuclear$xtab, nuclear$ytab, x.nucl, deg=2)
y.poly.4<-poly_est(nuclear$xtab, nuclear$ytab, x.nucl, deg=4) 


###################################################
### code chunk number 18: ex-npbr.rnw:265-270 (eval = FALSE)
###################################################
## plot(ytab~xtab, data=nuclear)
## lines(x.nucl, y.poly.2, lty=1, lwd=4, col="red") 
## lines(x.nucl, y.poly.4, lty=2, lwd=4, col="blue")
## legend("topleft",legend=c("degree=2", "degree=4"), 
##  col=c("red","blue"), lwd=4, lty=c(1,2))


###################################################
### code chunk number 19: ex-npbr.rnw:274-281
###################################################
.PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=name.file, width = 9, height = 7, pointsize = 14, bg = "white")
plot(ytab~xtab, data=nuclear)
lines(x.nucl, y.poly.2, lty=1, lwd=4, col="red") 
lines(x.nucl, y.poly.4, lty=2, lwd=4, col="blue")
legend("topleft",legend=c("degree=2", "degree=4"), 
 col=c("red","blue"), lwd=4, lty=c(1,2))
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 20: ex-npbr.rnw:375-377
###################################################
(kn.aic.mono<-quad_spline_est_kn(log(green$COST), log(green$OUTPUT), 
 x.green, cv=0, type="AIC"))


###################################################
### code chunk number 21: ex-npbr.rnw:381-383
###################################################
y.quad.1<-quad_spline_est(log(green$COST), log(green$OUTPUT), 
 x.green, kn=kn.aic.mono, cv=0)


###################################################
### code chunk number 22: ex-npbr.rnw:388-392
###################################################
(kn.bic.conca<-quad_spline_est_kn(log(green$COST), log(green$OUTPUT), 
 x.green, cv=1, type="BIC"))
y.quad.2<-quad_spline_est(log(green$COST), log(green$OUTPUT), 
 x.green, kn=kn.bic.conca, cv=1)


###################################################
### code chunk number 23: ex-npbr.rnw:395-398
###################################################
y.quad.3<-quad_spline_est(log(green$COST), log(green$OUTPUT), 
 x.green, cv=1, 
 all.dea=TRUE)


###################################################
### code chunk number 24: ex-npbr.rnw:402-409 (eval = FALSE)
###################################################
## plot(log(OUTPUT)~log(COST), data=green)
## lines(x.green, y.quad.1, lty=2, lwd=4, col="red")   
## lines(x.green, y.quad.2, lty=2, lwd=4, col="blue")
## lines(x.green, y.quad.3, lwd=4, lty=1)   
## legend("topleft", col=c("red","blue","black"), lty=c(2,2,1), 
##  legend=c("mono(kn=6)", "mono + concav (kn=1)", 
##  "mono + concav (kn=all DEA points)"), lwd=4, cex=0.8)


###################################################
### code chunk number 25: ex-npbr.rnw:413-422
###################################################
.PngNo <- .PngNo + 1; name.file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=name.file, width = 9, height = 7, pointsize = 14, bg = "white")
plot(log(OUTPUT)~log(COST), data=green)
lines(x.green, y.quad.1, lty=2, lwd=4, col="red")   
lines(x.green, y.quad.2, lty=2, lwd=4, col="blue")
lines(x.green, y.quad.3, lwd=4, lty=1)   
legend("topleft", col=c("red","blue","black"), lty=c(2,2,1), 
 legend=c("mono(kn=6)", "mono + concav (kn=1)", 
 "mono + concav (kn=all DEA points)"), lwd=4, cex=0.8)
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", name.file, "}\n\n", sep="")


###################################################
### code chunk number 26: ex-npbr.rnw:446-462
###################################################
simulate.data<- function(n, funs=1, betav=0.5)
{
 # internal function
 Fron<-function(x,funs)
 {
  if (funs==1) { return(exp (-5 + 10*x)/(1 + exp(-5 + 10*x)))} 
  if (funs==2) { return(sqrt(x))}                              
  if (funs==3) { return(x)}                                    
 }

  xtab <- runif(n, 0, 1)
  V <-rbeta(n, betav, betav)
  ytab <- Fron(xtab, funs)*V

  return(data.frame(xtab=xtab, ytab=ytab))
}


###################################################
### code chunk number 27: ex-npbr.rnw:468-476 (eval = FALSE)
###################################################
## N<-200 
## x.sim <- seq(0, 1, length.out=1000)
## y.dea<-matrix(0, N, 1000)
## for(k in 1:N)
## {
##  don<-simulate.data(25)
##  y.dea[k,]<-dea_est(don$xtab, don$ytab, x.sim, type="dea")
## }


###################################################
### code chunk number 28: ex-npbr.rnw:479-480 (eval = FALSE)
###################################################
## y.dea[is.infinite(y.dea)]<-0


###################################################
### code chunk number 29: ex-npbr.rnw:483-484 (eval = FALSE)
###################################################
## error.dea<-matrix(x.sim, N, 1000, byrow=TRUE)-y.dea


###################################################
### code chunk number 30: ex-npbr.rnw:489-492 (eval = FALSE)
###################################################
## (IBIAS2<-mean((x.sim-apply(y.dea,2,mean))^2))
## (IVAR<-mean((y.dea-matrix(apply(y.dea,2,mean),N,1000,byrow=TRUE))^2))
## (MISE<-IBIAS2+IVAR)


