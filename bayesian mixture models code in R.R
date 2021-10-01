#R Codes (For 4 component Normal Mixture Model)
#Creating Gibbs Sampler
gibbsnorm <- function (dat, k, niter, alpha = 1.28, beta = 0.36 * var(dat),
                       lam = mean(dat), tau = 2.6/(diff(range(dat)))^2, g = 1)
{
  rigamma <- function(n, a, b) { return(1/rgamma(n, shape = a, rate = b))
  }
  rdirichlet <- function(n, par) { 
    k = length(par)
    z =   array (0,dim=c (n, k)) 
    s = array (0, dim =   c (n, 1) ) 
    for (i in 1 : k) {
      z [,i]= rgamma (n,shape= par[i] ) 
      s = s + z [ , i]
    }
    for (i in 1:k) {
      z[,i] = z[, i]/s 
    }
    return(z)
  }
  n <- length(dat)
  mu <- rnorm(k, mean = mean(dat), sd = sd(dat)) 
  sig <- sd(dat)/k
  p <- rep(1/k, k)
  mixparam <- list(p = p, mu = mu, sig = sig) 
  z <- rep(0, k)
  nj <- z 
  sj <- z 
  sj2 <- z
  gibbsmu <- matrix(0, nrow = niter, ncol = k) 
  gibbssig <- gibbsmu
  gibbsp <- gibbsmu
  for (i in 1:niter) 
  { for (t in 1:n)
  { prob <- mixparam$p * dnorm(dat[t], mean = mixparam$mu,
                               sd = mixparam$sig)
  z[t] <- sample(x = 1:k, size = 1, prob = prob)
  }
    for	(j	in	1 : k) {
      nj [j] <-sum (z == j )
      sj [j] <-sum(as.numeric(z == j) * dat)
    }
    repeat {
      gibbsmu[i, ]<- rnorm(k, mean = (lam * tau + sj)/(nj + tau),
                           sd = sqrt(mixparam$sig^2/(tau + nj)))
      if (max(gibbsmu[i, ]) < max(dat) & min(gibbsmu[i, ]) > min(dat))
        break
    }
    mixparam$mu<-gibbsmu[ i , ]
    for(j	in	1 : k) {
      sj2[j] =   sum(as.numeric(z==j ) *  (dat - mixparam$mu[j] )^2)
    }
    gibbssig[i,]<- sqrt(rigamma(k ,alpha  +0.5 *(nj + 1), 
                                beta + 0.5* tau *(mixparam$mu-lam)^2+ 0.5*sj2))
    mixparam$sig <- gibbssig[i, ]
    gibbsp[i, ] <- rdirichlet(1, par = nj + g)
    
    mixparam$p <- gibbsp[i, ]
  }    
  data.frame(p = gibbsp, mu = gibbsmu, sigma = gibbssig)
} 

#Running the Code for Required number of Components and Iterations
m4<-gibbsnorm(Boston$medv,4,50000)

#Code for Mixture Density Using Created Dataframe
mixdensity <- function(x,mu1,mu2,mu3,mu4,s1,s2,s3,s4,p1,p2,p3,p4) {
  den<-((p1)*dnorm(x,mean=mu1,sd=s1))+((p2)*dnorm(x,mean=mu2,sd=s2))+((p3)*dnorm(x,mean=mu3,sd=s3))
  +((p4)*dnorm(x,mean=mu4,sd=s4))
  return(den)
}

#Plotting the Predictive Density over the Histogram of the Data
p1<-mean(m4$p.1)
p2<-mean(m4$p.2)
p3<-mean(m4$p.3)
p4<-mean(m4$p.4)

mu1<-mean(m4$mu.1)
mu2<-mean(m4$mu.2)
mu3<-mean(m4$mu.3)
mu4<-mean(m4$mu.4)

s1<-mean(m4$sigma.1)
s2<-mean(m4$sigma.2)
s3<-mean(m4$sigma.3)
s4<-mean(m4$sigma.4)
par(mfrow=c(1,1))
hist(Boston$medv,xlab="housing price",main="Predictive Density based on 4 component normal mixture",
     probability = TRUE,breaks=seq(0,50,by=2))
curve(mixdensity(x,mu1,mu2,mu3,mu4,s1,s2,s3,s4,p1,p2,p3,p4),range(0,50),add=TRUE)

#Plotting Traceplots
library(ggplot2)
p <- ggplot(m4, aes(x=seq(from=1,to=50000,by=1), y=mu.2)) +
  geom_line()
p

#Plotting Mean vs Sigma Graphs
m<-c(m42$mu.1,m42$mu.2,m42$mu.3,m42$mu.4)
s<-c(m42$sigma.1,m42$sigma.2,m42$sigma.3,m42$sigma.4)
par(mfrow=c(1,1))
plot(m,s)

#Finding DIC value for the fitted model
zm4=matrix(NA,nrow=50000,ncol=506)
for(i in 1:50000){
  for(j in 1:506)
  { zm4[i,j]<-mixdensity4(x[j],m4$mu.1[i],m4$mu.2[i],m4$mu.3[i],m42$mu.4[i],
                          m4$sigma.1[i],m4$sigma.2[i],m4$sigma.3[i],m4$sigma.4[i],
                          m4$p.1[i],m4$p.2[i],m4$p.3[i],m4$p.4[i])
  }}
logzm4<-log(zm4)
lik4<-c()
for(i in 1:50000){
  for(j in 1:506){
    lik4[i]=sum(logzm4[i,j])
  }}
mean(lik4)
l4<-c()
for(i in 1:50000){
  for(j in 1:506){
    l4[j]=sum(zm4[i,j])/50000
  }}
ll4<-log(l4)
mean(ll4)


DIC<-(-4*mean(lik4))+(2*sum(ll4))
