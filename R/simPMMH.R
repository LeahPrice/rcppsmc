len <- 500
sigv <- sqrt(10)
sigw <- sqrt(1)

x <- rep(0,len)
y <- rep(0,len)
x[1] <- rnorm(1,mean=0,sd=sqrt(5))
y[1] <- rnorm(1,mean=x[1]^2/20,sd=sigw)

for (i in 2:len){
  x[i] <- rnorm(1,mean=x[i-1]/2.0 + 25.0*x[i-1]/(1+x[i-1]^2) + 8*cos(1.2*i),sd=sigv)
  y[i] <- rnorm(1,mean=x[i]^2/20,sd=sigw)
  
}
y = as.data.frame(x=y)
names(y) <- toString(len)

write.csv(y,file="PMMH-data.csv",row.names = FALSE, quote=FALSE) 
