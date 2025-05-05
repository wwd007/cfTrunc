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


# Simple Pseudo-Observation Regression with Sequentially Truncated Data
seqTrun.POreg.cox <- function(data){
  ##vector to be used for estimation
  temp.dat <- data
  n <- dim(temp.dat)[1]
  l <- temp.dat[,1]
  x <- temp.dat[,2]
  r <- temp.dat[,3]
  t0 <- quantile(l, seq(0.1, 0.8, length.out=5))
  m<-length(t0)
  rho<-matrix(NA, n, m) #matrix to store pseudo obs
  
  #calculate pseudo obs
  w <- rep(1,length(x))
  F.npmle <- sapply(sort(t0),fun_test,l=l,x=x,r=r,w=w)
  
  fit.1<-F.npmle
  km.est<-1-fit.1
  for (i in 1:n){
    temp.dat.2<-temp.dat[-i,]
    l<-temp.dat.2[,1]
    x<-temp.dat.2[,2]
    r<-temp.dat.2[,3]
    w <- rep(1,length(x))
    fit.2<-sapply(sort(t0),fun_test,l=l,x=x,r=r,w=w)
    km.est.2<-1-fit.2
    for (j in 1:m){ 
      rho[i,j]<-n*km.est[j]-(n-1)*km.est.2[j]
    }
  }
###run a regression analysis using GEE with complementary log-log link function
#rho.fwd <- 1 - rho #1-rho = Pr(X>tau0-t0)
#rearrange the data into a long data set
dat.long <- NULL
for(j in 1:m){
  dat.long <- rbind(dat.long, cbind(temp.dat, pseudo = rho[,j], ipseudo=1-rho[,j], tpseudo = t0[j], id=1:nrow(temp.dat),row.names = NULL))
}
dat.long <- dat.long[order(dat.long$id),]
#fit a Cox model using GEE
fit <- geese(ipseudo~as.factor(tpseudo)+Z,data=dat.long,scale.fix=TRUE,family=gaussian,jack=TRUE, mean.link="cloglog",corstr="independence")
output<-summary(fit) 
b0.est<-output$mean[nrow(output$mean),1]#point estimate of b0
b0.ajs.se<-output$mean[nrow(output$mean),3] 
out <- list(b0.est, b0.ajs.se)
names(out) <- c("Coefficient estimate","Standard Error")
out
}

##
#Coefficient estimate: The estimate of the regression coefficients.
#Standard Error: The standard error estimates using approximate jackknife variance estimate in geepack.


# Modified Pseudo-Observation Regression with Sequentially Truncated Data
seqTrun.modPOreg.cox <- function(data){
  ##vector to be used for estimation
  temp.dat <- data
  n <- dim(temp.dat)[1]
  l <- temp.dat[,1]
  x <- temp.dat[,2]
  r <- temp.dat[,3]
  t0 <- quantile(l, seq(0.1, 0.8, length.out=5))
  m<-length(t0)
  rho<-matrix(NA, n, m) #matrix to store pseudo obs
  
  #calculate pseudo obs
  w <- rep(1,length(x))
  S.r.x <- fn_est(t0[1],l,x,r,w)[[3]] 
  po.wt <- 1/pmax(S.r.x, 1e-4)#weight in binomial reg
  
  for (j in 1:m){ 
    rho[,j]<-(l>t0[j])*1
  }
  ###run a regression analysis using GEE with complementary log-log link function
  #rho.fwd <- 1 - rho #1-rho = Pr(X>tau0-t0)
  #rearrange the data into a long data set
  dat.long <- NULL
  for(j in 1:m){
    dat.long <- rbind(dat.long, cbind(temp.dat, irho=1-rho[,j], trho = t0[j], po.wt=po.wt, id=1:nrow(temp.dat),row.names = NULL))
  }
  dat.long <- dat.long[order(dat.long$id),]
  #fit a Cox model using GEE
  fit <- geese(irho~as.factor(trho)+Z,data=dat.long,scale.fix=TRUE,family=binomial,jack=TRUE, weights=po.wt, mean.link="cloglog",corstr="independence")
  output<-summary(fit) 
  b0.est<-output$mean[nrow(output$mean),1]#point estimate of b0
  b0.ajs.se<-output$mean[nrow(output$mean),3] 
  out <- list(b0.est, b0.ajs.se)
  names(out) <- c("Coefficient estimate","Standard Error")
  out
}

##
#Coefficient estimate: The estimate of the regression coefficients.
#Standard Error: The standard error estimates using approximate jackknife variance estimate in geepack.

