library(survival)
library(geepack)

###function to estimate CDF of L
###function to estimate CDF of L
fn_est <- function(a,l,x,r,w){
  ## to estimate P(L<u)
  ## L<X<R,
  ## observable: (L,X,R)|X<R
  n<-length(x)
  fit <-  survival::survfit(survival::Surv(x, r, event=rep(1,n)) ~ 1,weights = w,
                            timefix = FALSE)
  s.r.x <- summary(fit,time=x)$surv
  s.r.x <- s.r.x[match(x,sort(x))]
  F.l <- 1/sum(w/s.r.x)*sum((l<=a)*w/s.r.x)
  #  bet <- sum(w/s.r.x)*n
  bet <- 1-n/sum(w/s.r.x) #truncation probability
  return(list(F.l,bet,s.r.x))
}

###functions to apply
fun_test <- function(a,l,x,r,w){
  fit <- fn_est(a,l,x,r,w)
  return(fit[[1]])
}

###function to estimate CDF of L at all the obs points
fn_est_all <- function(l,x,r,w){
  ## to estimate P(L<u)
  ## L<X<R,
  ## observable: (L,X,R)|X<R
  n<-length(x)
  l.sort<-sort(l)
  fit <-  survival::survfit(survival::Surv(x, r, event=rep(1,n)) ~ 1,weights = w,
                            timefix = FALSE)
  s.r.x <- summary(fit,time=x)$surv
  s.r.x <- s.r.x[match(x,sort(x))]
  # l.rank<- rank(l, ties.method="max") #obtain sum(l<=a)
  F.l<-rep(0,n)
  for(i in 1:n){
    F.l[i]<- 1/sum(w/s.r.x)*sum((l<=l.sort[i])*w/s.r.x)
  }
  #  bet <- sum(w/s.r.x)*n
  bet <- 1-n/sum(w/s.r.x) #truncation probability
  return(list(F.l,bet,s.r.x))
}

###
#function to estimate E(log T)
fun.E.logT<-function(tt,surv){
  ## tt is already sorted from small to large
  m<-length(tt)
  h <- log(tt[-1])-log(tt[-m])
  surv.L <- surv[-m]
  surv.R<- surv[-1]
  est<-log(tt[1]) +  sum((surv.L+surv.R)*h/2)   
  return(est)                 
}

# Simple Pseudo-Observation Regression with Sequentially Truncated Data
seqTrun.POreg.aft <- function(data){
  ##vector to be used for estimation
  temp.dat <- data
  n <- dim(temp.dat)[1]
  l <- temp.dat[,1]
  x <- temp.dat[,2]
  r <- temp.dat[,3]
  t0 <- quantile(l, seq(0.1, 0.8, length.out=5))
  po.logT<-rep(NA, n) #vector to store pseudo obs
  
  #calculate pseudo obs
  tt<-sort(l)
  w <- rep(1,length(x))
  fit.surv<- 1-fn_est_all(l,x,r,w)[[1]]
  theta<-fun.E.logT(tt, fit.surv)
  
  for (i in 1:n){
    dat2<-dat[-i,]
    w<-rep(1,(n-1))
    l=dat2$L
    x=dat2$X
    r=dat2$R
    tt2 <- sort(dat2$L)
    fit2.surv<-  1-fn_est_all(l,x,r,w)[[1]]
    theta.i<-fun.E.logT(tt2, fit2.surv)  
    po.logT[i]<-n*theta-(n-1)*theta.i
  }

#rearrange the data into a long data set
dat.long <- NULL
#for(j in 1:m){
  dat.long <- rbind(dat.long, cbind(temp.dat, id=1:nrow(temp.dat),row.names = NULL))
#}
dat.long <- dat.long[order(dat.long$id),]
#fit an AFT model using geese() function
fit <- geese(po.logT~Z,data=dat.long,scale.fix=TRUE,family=gaussian,jack=TRUE, mean.link="identity",corstr="independence")
output<-summary(fit) 
b0.est<-output$mean[nrow(output$mean),1]#point estimate of b0
b0.ajs.se<-output$mean[nrow(output$mean),3] 
out <- list(b0.est, b0.ajs.se)
names(out) <- c("Coefficient estimate","Standard Error")
out
}

##
#Coefficient estimate: The estimate of the regression coefficients.
#Standard Error: The standard error estimates.


# Modified Pseudo-Observation Regression with Sequentially Truncated Data
seqTrun.modPOreg.aft <- function(data){
  ##vector to be used for estimation
  temp.dat <- data
  n <- dim(temp.dat)[1]
  l <- temp.dat[,1]
  x <- temp.dat[,2]
  r <- temp.dat[,3]
  logL=log(l)
  t0 <- quantile(l, seq(0.1, 0.8, length.out=5))
  
  w <- rep(1,length(x))
  S.r.x <- fn_est(t0[1],l,x,r,w)[[3]] 
  po.wt <- 1/pmax(S.r.x, 1e-4)#weight in ipw reg
  #rearrange the data into a long data set
  dat.long <- NULL
  dat.long <- rbind(dat.long, cbind(temp.dat, id=1:nrow(temp.dat), po.wt=po.wt, row.names = NULL))
  dat.long <- dat.long[order(dat.long$id),]
  #fit an AFT model using geese() function
  fit <- geese(logL~Z,data=dat.long,scale.fix=TRUE,family=gaussian,jack=TRUE, weights=po.wt, mean.link="identity",corstr="independence")
  output<-summary(fit) 
  # b0.est<-output$mean[nrow(output$mean),1]#point estimate of b0
  # b0.ajs.se<-output$mean[nrow(output$mean),3] 
  # out <- list(b0.est, b0.ajs.se)
  # names(out) <- c("Coefficient estimate","Standard Error")
  # out
}

##
#Coefficient estimate: The estimate of the regression coefficients.
#Standard Error: The standard error estimates based on sandwich variance estimator.

