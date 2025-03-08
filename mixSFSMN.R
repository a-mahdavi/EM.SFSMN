 # A function to fit a g component mixtures of SFSMN distributions
 # family argument detrmine the fitted distribution including: "SFN","SFT","SFL","SFSL","SFEP","SFCN" 

mix.SFSMN <- function(y, g=1, w=1, mu, s, del, nu=1, la, family="SFT", iter.max=100, tol=10^-6, get.init = TRUE, group=T){  
   begin <- proc.time()[3]
	skewness <- function (x, na.rm = FALSE) 
		{
    	if (is.matrix(x)) 
        apply(x, 2, skewness, na.rm = na.rm)
    	else if (is.vector(x)) {
        if (na.rm) 
            x <- x[!is.na(x)]
        n <- length(x)
        (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  	  	}
    	else if (is.data.frame(x)) 
        sapply(x, skewness, na.rm = na.rm)
    	else skewness(as.vector(x), na.rm = na.rm)
		}
	if(family=="SFN"){
	dSFN <- function(y, mu, s, del, la)
		  dnorm(abs(y-mu)/s-del)/(s*pnorm(del))*pnorm(la*(y-mu)/s)
	dmixSFN <- function(y, w, mu, s, del, la){
			d <- 0 ; g <- length(w)
			for ( j in 1:g)
			d <- d + w[j]*dSFN(y, mu[j], s[j], del[j], la[j])
			return(d) }
  n <- length(y)   ;     dif <- 1   ;     count <- 0 
  if (get.init == TRUE) {
	if(g==1){
		mu <- runif(1,summary(y)[2],summary(y)[5])
		s <- runif(1,0,2)*sd(y)
		la <- skewness(y)+runif(1,-1,1)
		del <- runif(1,-1,1)
		}
		else {
                 init <- kmeans(y, g,  algorithm="Hartigan-Wong")
                 w <- init$size/n ; del <-la <- rep(0,g)
                 mu <- as.vector(init$centers)
                 s <- sqrt(init$withinss/init$size) 
			nu <- NULL
			for( j in 1:g)
			la[j] <- skewness(y[init$cluster==j])				
			  }}
   LL <- 1  ; al=la/s
  while ((dif > tol) && (count <= iter.max)) {
  z.hat  <- matrix(0,n,g)
  # E step
  for (j in 1:g){
  z.hat[,j] <- w[j]*dSFN(y, mu[j], s[j], del[j], la[j])/dmixSFN(y, w, mu, s , del, la)
	# eta <- abs(y-mu[j])/s[j]-del[j]
	gamh1 <- rep(1,n) ; gamh2 <- al[j]*(y-mu[j])+ dnorm(al[j]*(y-mu[j]))/pnorm(al[j]*(y-mu[j]))
  # MCE steps
  w[j] <- sum(z.hat[,j])/n
	mu[j] <- (sum(z.hat[,j]*gamh1*(y-s[j]*del[j]*sign(y-mu[j])))-al[j]*s[j]^2*sum(z.hat[,j]*gamh2)+al[j]^2*s[j]^2*sum(z.hat[,j]*y))/sum(z.hat[,j]*(gamh1+s[j]^2*al[j]^2))
	a <- del[j]*sum(z.hat[,j]*gamh1*abs(y-mu[j]))+la[j]*sum(z.hat[,j]*gamh2*(y-mu[j]))
	b <- sum(z.hat[,j]*gamh1*(y-mu[j])^2)+la[j]^2*sum(z.hat[,j]*(y-mu[j])^2)
	s1 <- (-a+sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j])) ; s2 <- (-a-sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j]))
		if ( s1 >0) 
		s[j] <- s1 else 
		s[j] <- s2
	 del[j] <- (sum(z.hat[,j]*gamh1*abs(y-mu[j])/s[j])-sum(z.hat[,j])*dnorm(del[j])/pnorm(del[j]))/sum(z.hat[,j]*gamh1)
	la [j] <- sum(z.hat[,j]*gamh2*(y-mu[j])/s[j])/sum(z.hat[,j]*((y-mu[j])/s[j])^2)
	}
  LL.new <- sum(log(dmixSFN(y,w,mu,s,del,la))) # log-likelihood function
  count <- count +1 
  dif <- abs(LL.new/LL-1)
	LL <- LL.new; al=la/s
	cat('iter =', count, '\tloglike =', LL.new, '\n')
  } }
	if(family=="SFT"){
	dSFt <- function(y, mu, s, del, nu, la)
		  dt(abs(y-mu)/s-del,nu)/(s*pt(del,nu))*pnorm(la*(y-mu)/s)
	dmixSFt <- function(y, w, mu, s, del, nu, la){
			d <- 0 ; g <- length(w)
			for ( j in 1:g)
			d <- d + w[j]*dSFt(y, mu[j], s[j], del[j], nu[j], la[j])
			return(d) }
  n <- length(y)   ;     dif <- 1   ;     count <- 0 
  if (get.init == TRUE) {	
		if(g==1){
		mu <- runif(1,summary(y)[2],summary(y)[5])
		s <- runif(1,0,2)*sd(y)
		la <- skewness(y)+runif(1,-1,1)
		del <- runif(1,-1,1)
		nu <- runif(1,1,10)
		} else{
                 init <- kmeans(y, g,  algorithm="Hartigan-Wong")
                 w <- init$size/n ; del <-la <- rep(0,g)
                 mu <- as.vector(init$centers)
                 s <- sqrt(init$withinss/init$size) 
			nu <- runif(g,1,10)
			for( j in 1:g)
			la[j] <- skewness(y[init$cluster==j])				
			  }}
   LL <- 1  ; al=la/s
  while ((dif > tol) && (count <= iter.max)) {
  z.hat  <- matrix(0,n,g)
  # E step
  for (j in 1:g){
  z.hat[,j] <- w[j]*dSFt(y, mu[j], s[j], del[j], nu[j], la[j])/dmixSFt(y, w, mu, s , del, nu, la)
	eta <- abs(y-mu[j])/s[j]-del[j]
	gamh1 <- (nu[j]+1)/(nu[j]+eta^2); gamh2 <- al[j]*(y-mu[j])+ dnorm(al[j]*(y-mu[j]))/pnorm(al[j]*(y-mu[j]))
  # MCE steps
  w[j] <- sum(z.hat[,j])/n
	mu[j] <- (sum(z.hat[,j]*gamh1*(y-s[j]*del[j]*sign(y-mu[j])))-al[j]*s[j]^2*sum(z.hat[,j]*gamh2)+al[j]^2*s[j]^2*sum(z.hat[,j]*y))/sum(z.hat[,j]*(gamh1+s[j]^2*al[j]^2))
	a <- del[j]*sum(z.hat[,j]*gamh1*abs(y-mu[j]))+la[j]*sum(z.hat[,j]*gamh2*(y-mu[j]))
	b <- sum(z.hat[,j]*gamh1*(y-mu[j])^2)+la[j]^2*sum(z.hat[,j]*(y-mu[j])^2)
	s1 <- (-a+sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j])) ; s2 <- (-a-sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j]))
		if ( s1 >0) 
		s[j] <- s1 else 
		s[j] <- s2
	del[j] <- (sum(z.hat[,j]*gamh1*abs(y-mu[j])/s[j])-sum(z.hat[,j])*dt(del[j],nu[j])/pt(del[j],nu[j]))/sum(z.hat[,j]*gamh1)
	la [j] <- sum(z.hat[,j]*gamh2*(y-mu[j])/s[j])/sum(z.hat[,j]*((y-mu[j])/s[j])^2)
	nu[j] <- optim(nu[j],function(x){
			nu[j] <- x
		-sum(log(dmixSFt(y,w,mu,s,del,nu,la)))
		},method="L-BFGS-B",lower=0.01,upper=30)$par
	}
  LL.new <- sum(log(dmixSFt(y,w,mu,s,del,nu,la))) # log-likelihood function
  count <- count +1 
  dif <- abs(LL.new/LL-1)
	LL <- LL.new ;	al=la/s
	cat('iter =', count, '\tloglike =', LL.new, '\n')
  } }
	if(family=="SFL"){
	dL <- function(y)
		exp(-abs(y))/2
	dSFL <- function(y, mu, s, del, la)
		  dL(abs(y-mu)/s-del)/(s*integrate(dL,-Inf,del,stop.on.error=F)$value)*pnorm(la*(y-mu)/s)
	dmixSFL <- function(y, w, mu, s, del, la){
			d <- 0 ; g <- length(w)
			for ( j in 1:g)
			d <- d + w[j]*dSFL(y, mu[j], s[j], del[j], la[j])
			return(d) }
  n <- length(y)   ;     dif <- 1   ;     count <- 0 
  if (get.init == TRUE) {
		if(g==1){
		mu <- runif(1,summary(y)[2],summary(y)[5])
		s <- runif(1,0,2)*sd(y)
		la <- skewness(y)+runif(1,-1,1)
		del <- runif(1,-1,1)
		nu <- NULL
		} else{
                 init <- kmeans(y, g,  algorithm="Hartigan-Wong")
                 w <- init$size/n ; del <-la <- rep(0,g)
                 mu <- as.vector(init$centers)
                 s <- sqrt(init$withinss/init$size) 
			for( j in 1:g)
			la[j] <- skewness(y[init$cluster==j])
			nu <- NULL
			  }}
   LL <- 1 ;al<-la/s
  while ((dif > tol) && (count <= iter.max)) {
  z.hat  <- matrix(0,n,g)
  # E step
  for (j in 1:g){
  z.hat[,j] <- w[j]*dSFL(y, mu[j], s[j], del[j], la[j])/dmixSFL(y, w, mu, s , del, la)
	eta <- abs(y-mu[j])/s[j]-del[j]
	gamh1 <- 1/abs(eta)*besselK(abs(eta),1/2)/besselK(abs(eta),-1/2)
	gamh1[is.na(gamh1)] <- 0 ; gamh2 <- al[j]*(y-mu[j])+ dnorm(al[j]*(y-mu[j]))/pnorm(al[j]*(y-mu[j]))
  # MCE steps
  w[j] <- sum(z.hat[,j])/n
mu[j] <- (sum(z.hat[,j]*gamh1*(y-s[j]*del[j]*sign(y-mu[j])))-al[j]*s[j]^2*sum(z.hat[,j]*gamh2)+al[j]^2*s[j]^2*sum(z.hat[,j]*y))/sum(z.hat[,j]*(gamh1+s[j]^2*al[j]^2))
	a <- del[j]*sum(z.hat[,j]*gamh1*abs(y-mu[j]))+la[j]*sum(z.hat[,j]*gamh2*(y-mu[j]))
	b <- sum(z.hat[,j]*gamh1*(y-mu[j])^2)+la[j]^2*sum(z.hat[,j]*(y-mu[j])^2)	
     s1 <- (-a+sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j])) ; s2 <- (-a-sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j]))
		if ( s1 >0) 
		s[j] <- s1 else 
		s[j] <- s2
	 del[j] <- (sum(z.hat[,j]*gamh1*abs(y-mu[j])/s[j])-sum(z.hat[,j])*dL(del[j])/integrate(dL,-Inf,del[j],stop.on.error=F)$value)/sum(z.hat[,j]*gamh1)
	la [j] <- sum(z.hat[,j]*gamh2*(y-mu[j])/s[j])/sum(z.hat[,j]*((y-mu[j])/s[j])^2)
		}
  LL.new <- sum(log(dmixSFL(y,w,mu,s,del,la))) # log-likelihood function
  count <- count +1 
  dif <- abs(LL.new/LL-1)
	LL <- LL.new; al=la/s
	cat('iter =', count, '\tloglike =', LL.new, '\n')
  } }
	if(family=="SFSL"){
	dSL <- function(y,nu){
		y[which(y==0)] <- 10^-10
		nu*gamma(nu+1/2)*2^nu*pgamma(y^2/2,nu+1/2)/(sqrt(pi)*abs(y)^(2*nu+1))
		}
	dSFSL <- function(y, mu, s, del, nu, la)
		  dSL(abs(y-mu)/s-del,nu)/(s*integrate(dSL,-Inf,del,nu=nu,stop.on.error=F)$value)*pnorm(la*(y-mu)/s)
	dmixSFSL <- function(y, w, mu, s, del, nu, la){
			d <- 0 ; g <- length(w)
			for ( j in 1:g)
			d <- d + w[j]*dSFSL(y, mu[j], s[j], del[j], nu[j], la[j])
			return(d) }
  n <- length(y)   ;     dif <- 1   ;     count <- 0 
  if (get.init == TRUE) {
		if(g==1){
		mu <- runif(1,summary(y)[2],summary(y)[5])
		s <- runif(1,0,2)*sd(y)
		la <- skewness(y)+runif(1,-1,1)
		del <- runif(1,-1,1)
		nu <- runif(1,1,10)
		} else{
                 init <- kmeans(y, g,  algorithm="Hartigan-Wong")
                 w <- init$size/n ; del <-la <- rep(0,g)
                 mu <- as.vector(init$centers)
                 s <- sqrt(init$withinss/init$size) 
			nu <- runif(g,1,10)
			for( j in 1:g)
			la[j] <- skewness(y[init$cluster==j])	
			  }}
   LL <- 1; al <- la/s
  while ((dif > tol) && (count <= iter.max)) {
  z.hat  <- matrix(0,n,g)
  # E step
  for (j in 1:g){
  z.hat[,j] <- w[j]*dSFSL(y, mu[j], s[j], del[j], nu[j], la[j])/dmixSFSL(y, w, mu, s , del, nu, la)
	eta <- abs(y-mu[j])/s[j]-del[j]
	gamh1 <-  ((2*nu[j]+1)/eta^2)*(pgamma(eta^2/2,nu[j]+3/2)/pgamma(eta^2/2,nu[j]+1/2))
	gamh1[which(gamh1==0)] <- (2*nu[j]+1)/(2*nu[j]+3)
	gamh2 <- al[j]*(y-mu[j])+ dnorm(al[j]*(y-mu[j]))/pnorm(al[j]*(y-mu[j]))
  # MCE steps
  	w[j] <- sum(z.hat[,j])/n
	mu[j] <- (sum(z.hat[,j]*gamh1*(y-s[j]*del[j]*sign(y-mu[j])))-al[j]*s[j]^2*sum(z.hat[,j]*gamh2)+al[j]^2*s[j]^2*sum(z.hat[,j]*y))/sum(z.hat[,j]*(gamh1+s[j]^2*al[j]^2))
	a <- del[j]*sum(z.hat[,j]*gamh1*abs(y-mu[j]))+la[j]*sum(z.hat[,j]*gamh2*(y-mu[j]))
	b <- sum(z.hat[,j]*gamh1*(y-mu[j])^2)+la[j]^2*sum(z.hat[,j]*(y-mu[j])^2)
	s1 <- (-a+sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j])) ; s2 <- (-a-sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j]))
		if ( s1 >0) 
		s[j] <- s1 else 
		s[j] <- s2
	del[j] <- (sum(z.hat[,j]*gamh1*abs(y-mu[j])/s[j])-sum(z.hat[,j])*dSL(del[j],nu[j])/integrate(dSL,-Inf,del[j],nu=nu[j],stop.on.error=F)$value)/sum(z.hat[,j]*gamh1)
	la [j] <- sum(z.hat[,j]*gamh2*(y-mu[j])/s[j])/sum(z.hat[,j]*((y-mu[j])/s[j])^2)
	nu[j] <- optim(nu[j],function(x){
			nu[j] <- x
		-sum(log(dmixSFSL(y,w,mu,s,del,nu,la)))
		},method="L-BFGS-B",lower=.01,upper=30)$par
	}
  LL.new <- sum(log(dmixSFSL(y,w,mu,s,del, nu, la))) # log-likelihood function
  count <- count +1 
  dif <- abs(LL.new/LL-1)
	LL <- LL.new; al=la/s
	cat('iter =', count, '\tloglike =', LL.new, '\n')
  } }
if(family=="SFEP"){
	dEP <- function(y,nu)
		nu*exp(-1/2*abs(y)^(2*nu))/(sqrt(2^(1/nu))*gamma(1/(2*nu)))
	dSFEP <- function(y, mu, s, del, nu, la)
		  dEP(abs(y-mu)/s-del,nu)/(s*integrate(dEP,-Inf,del,nu=nu,stop.on.error=F)$value)*pnorm(la*(y-mu)/s)
	dmixSFEP <- function(y, w, mu, s, del, nu, la){
			d <- 0 ; g <- length(w)
			for ( j in 1:g)
			d <- d + w[j]*dSFEP(y, mu[j], s[j], del[j], nu[j], la[j])
			return(d) }
  n <- length(y)   ;     dif <- 1   ;     count <- 0 
  if (get.init == TRUE) {
		if(g==1){
		mu <- runif(1,summary(y)[2],summary(y)[5])
		s <- runif(1,0,2)*sd(y)
		la <- skewness(y)+runif(1,-1,1)
		del <- runif(1,-1,1)
		nu <- runif(1,0,1)
		} else{
                 init <- kmeans(y, g,  algorithm="Hartigan-Wong")
                 w <- init$size/n ; del <-la <- rep(0,g)
                 mu <- as.vector(init$centers)
                 s <- sqrt(init$withinss/init$size) 
			nu <- runif(g,0,1)
			for( j in 1:g)
			la[j] <- skewness(y[init$cluster==j])
			  }}
   LL <- 1 ; al <- la/s
  while ((dif > tol) && (count <= iter.max)) {
  z.hat  <- matrix(0,n,g)
  # E step
  for (j in 1:g){
  z.hat[,j] <- w[j]*dSFEP(y, mu[j], s[j], del[j], nu[j], la[j])/dmixSFEP(y, w, mu, s , del, nu, la)
	eta <- abs(y-mu[j])/s[j]-del[j]
	gamh1 <- nu[j]*(abs(eta))^(2*nu[j]-2); gamh2 <- al[j]*(y-mu[j])+ dnorm(al[j]*(y-mu[j]))/pnorm(al[j]*(y-mu[j]))
  # MCE steps
 	w[j] <- sum(z.hat[,j])/n
	mu[j] <- (sum(z.hat[,j]*gamh1*(y-s[j]*del[j]*sign(y-mu[j])))-al[j]*s[j]^2*sum(z.hat[,j]*gamh2)+al[j]^2*s[j]^2*sum(z.hat[,j]*y))/sum(z.hat[,j]*(gamh1+s[j]^2*al[j]^2))
	a <- del[j]*sum(z.hat[,j]*gamh1*abs(y-mu[j]))+la[j]*sum(z.hat[,j]*gamh2*(y-mu[j]))
	b <- sum(z.hat[,j]*gamh1*(y-mu[j])^2)+la[j]^2*sum(z.hat[,j]*(y-mu[j])^2)
	s1 <- (-a+sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j])) ; s2 <- (-a-sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j]))
		if ( s1 >0) 
		s[j] <- s1 else 
		s[j] <- s2
	del[j] <- (sum(z.hat[,j]*gamh1*abs(y-mu[j])/s[j])-sum(z.hat[,j])*dEP(del[j],nu[j])/integrate(dEP,-Inf,del[j],nu=nu[j],stop.on.error=F)$value)/sum(z.hat[,j]*gamh1)
	la [j] <- sum(z.hat[,j]*gamh2*(y-mu[j])/s[j])/sum(z.hat[,j]*((y-mu[j])/s[j])^2)
	nu[j] <- optim(nu[j],function(x){
			 nu[j] <- x
		-sum(log(dmixSFEP(y,w,mu,s,del,nu, la)))
		},method="L-BFGS-B",lower=0.01,upper=.99)$par
	}
  LL.new <- sum(log(dmixSFEP(y,w,mu,s,del,nu, la))) # log-likelihood function
  count <- count +1 
  dif <- abs(LL.new/LL-1)
	LL <- LL.new; al=la/s
	cat('iter =', count, '\tloglike =', LL.new, '\n')
  } }
	if(family=="SFCN"){
	dCN <- function(y,nu1,nu2)
		nu1*nu2^(1/2)*dnorm(nu2^(1/2)*y)+(1-nu1)*dnorm(y)
	dSFCN <- function(y, mu, s, del, nu1, nu2, la)
		  dCN(abs(y-mu)/s-del,nu1,nu2)/(s*integrate(dCN,-Inf,del,nu1=nu1,nu2=nu2,stop.on.error=F)$value)*pnorm(la*(y-mu)/s)
	dmixSFCN <- function(y, w, mu, s, del, nu1, nu2, la){
			d <- 0 ; g <- length(w)
			for ( j in 1:g)
			d <- d + w[j]*dSFCN(y, mu[j], s[j], del[j], nu1[j], nu2[j], la[j])
			return(d) }
  n <- length(y)   ;     dif <- 1   ;     count <- 0 
  if (get.init == TRUE) {
		if(g==1){
		mu <- runif(1,summary(y)[2],summary(y)[5])
		s <- runif(1,0,2)*sd(y)
		la <- skewness(y)+runif(1,-1,1)
		del <- runif(1,-1,1)
		nu1 <- runif(1,0,1);nu2 <- runif(1,0,1)
		} else{
                 init <- kmeans(y, g,  algorithm="Hartigan-Wong")
                 w <- init$size/n ; del <-la <- rep(0,g)
                 mu <- as.vector(init$centers)
                 s <- sqrt(init$withinss/init$size) 
			nu1 <- runif(g,0,1);nu2 <- runif(g,0,1)
			for( j in 1:g)
			la[j] <- skewness(y[init$cluster==j])
			  }}
	else{
		nu1 <- nu[[1]] ; nu2 <- nu[[2]] }
   LL <- 1 ; al <- la/s
  while ((dif > tol) && (count <= iter.max)) {
  z.hat  <- matrix(0,n,g)
  # E step
  for (j in 1:g){
  z.hat[,j] <- w[j]*dSFCN(y, mu[j], s[j], del[j], nu1[j], nu2[j], la[j])/dmixSFCN(y, w, mu, s , del, nu1, nu2, la)
	eta <- abs(y-mu[j])/s[j]-del[j]
	gamh1 <- (1-nu1[j]+nu1[j]*nu2[j]^(3/2)*exp((1-nu2[j])*eta^2/2))/(1-nu1[j]+nu1[j]*nu2[j]^(1/2)*exp((1-nu2[j])*eta^2/2))
	gamh2 <- al[j]*(y-mu[j])+ dnorm(al[j]*(y-mu[j]))/pnorm(al[j]*(y-mu[j]))
  # MCE steps
  w[j] <- sum(z.hat[,j])/n
	mu[j] <- (sum(z.hat[,j]*gamh1*(y-s[j]*del[j]*sign(y-mu[j])))-al[j]*s[j]^2*sum(z.hat[,j]*gamh2)+al[j]^2*s[j]^2*sum(z.hat[,j]*y))/sum(z.hat[,j]*(gamh1+s[j]^2*al[j]^2))
	a <- del[j]*sum(z.hat[,j]*gamh1*abs(y-mu[j]))+la[j]*sum(z.hat[,j]*gamh2*(y-mu[j]))
	b <- sum(z.hat[,j]*gamh1*(y-mu[j])^2)+la[j]^2*sum(z.hat[,j]*(y-mu[j])^2)
	s1 <- (-a+sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j])) ; s2 <- (-a-sqrt(a^2+4*sum(z.hat[,j])*b))/(2*sum(z.hat[,j]))
		if ( s1 >0) 
		s[j] <- s1 else 
		s[j] <- s2
	del[j] <- (sum(z.hat[,j]*gamh1*abs(y-mu[j])/s[j])-sum(z.hat[,j])*dCN(del[j],nu1[j],nu2[j])/integrate(dCN,-Inf,del[j],nu1=nu1[j],nu2=nu2[j],stop.on.error=F)$value)/sum(z.hat[,j]*gamh1)
	la [j] <- sum(z.hat[,j]*gamh2*(y-mu[j])/s[j])/sum(z.hat[,j]*((y-mu[j])/s[j])^2)
	nu1 <- optim(nu1,function(x){
		-sum(log(dmixSFCN(y,w,mu,s,del,x,nu2, la)))
		},method="L-BFGS-B",lower=0.01,upper=.99)$par
	nu2 <- optim(nu2,function(x){
		-sum(log(dmixSFCN(y,w,mu,s,del,nu1,x, la)))
		},method="L-BFGS-B",lower=0.01,upper=.99)$par 
	}
  LL.new <- sum(log(dmixSFCN(y,w,mu,s,del, nu1 ,nu2, la))) # log-likelihood function
  count <- count +1 
  dif <- abs(LL.new/LL-1)
	LL <- LL.new; al=la/s
	cat('iter =', count, '\tloglike =', LL.new, '\n')
  }
	nu <- list(nu1,nu2)
 }
  if(family=="SFL"|family=="SFN"){
  aic <- -2 * LL.new + 2 * (4*g+g-1)
  bic <- -2 * LL.new + log(n) * (4*g+g-1)
  edc <- -2 * LL.new + 0.2*sqrt(n) *(4*g+g-1)
	 }
  if(family=="SFCN"){
  aic <- -2 * LL.new + 2 * (6*g+g-1)
  bic <- -2 * LL.new + log(n) * (6*g+g-1) 
  edc <- -2 * LL.new + 0.2*sqrt(n) *(6*g+g-1)
 }
 else {
  aic <- -2 * LL.new + 2 * (5*g+g-1)
  bic <- -2 * LL.new + log(n) * (5*g+g-1) 
  edc <- -2 * LL.new + 0.2*sqrt(n) *(5*g+g-1)
	}
  end <- proc.time()[3]
  time <- end-begin
  obj.out <- list(w=w, mu=mu, s=s, del=del, nu=nu, la=la , loglik=LL.new, aic=aic, bic=bic, edc=edc, iter=count,elapsed=as.numeric(time),group = apply(z.hat, 
                  1, which.max))
  if (group==FALSE)
  obj.out <- obj.out[names(obj.out)!="group"]
  obj.out
  }




