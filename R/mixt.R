#' @import cubature

#' @import mvtnorm

#' @import pryr

#' @import Rmixmod

var.grid.pond=function(bin)
{
  l=nrow(bin$grid)
  g=apply(bin$grid[2:(l-1),],1,mean)
  m=sum(g*bin$bin[2:(l-1)])/sum(bin$bin[2:(l-1)])
  return(sum(((g-m)^2)*bin$bin[2:(l-1)])/sum(bin$bin[2:(l-1)]))
}

#' @export
log_sum=function(x) #x=c(log(a),log(b)), dà log(a+b)
{
  l=which.max(x)
  return((x[l])+log(1+sum(exp(x[-l]-x[l]))))
}

#' @export
log_diff=function(x) #x=c(log(a),log(b)), dà log(a-b)
{
  l=which.max(x)
  return((x[l])+log(1-sum(exp(x[-l]-x[l]))))
}

#' @export
logdiff_mixt=function(x,pi,mu,v,lower.tail=F)
{
  if(lower.tail==T)
  {
    d=pnorm(x[2],mu,sqrt(v),log.p = T,lower.tail = lower.tail)+log(pi)
    e=pnorm(x[1],mu,sqrt(v),log.p = T,lower.tail = lower.tail)+log(pi)
    f=apply(cbind(d,e),1,log_diff)
    log_sum(f)} else {
      d=pnorm(x[1],mu,sqrt(v),log.p = T,lower.tail = lower.tail)+log(pi)
      e=pnorm(x[2],mu,sqrt(v),log.p = T,lower.tail = lower.tail)+log(pi)
      f=apply(cbind(d,e),1,log_diff)
      log_sum(f)}
}
ellipse=function(x,mu,v,a)
{
  ((x-mu))%*%solve(v)%*%((x-mu))-qchisq(a,length(mu))
}

#' Generates a univariate gaussian mixture
#'
#' Generates a univariate gaussian mixture
#' @param p Needs to a vector of probabilities
#' @param mu Needs to be a vector of means
#' @param v Needs to be a vector of variances
#' @param seed Fix a random seed (optional)
#' @export
gen_mixt=function(n,p,mu,v,seed)
{
  if(missing(seed)==F) set.seed(seed)
  ind=rmultinom(n,1,prob=p)
  ind=apply(ind,2,function(x) which(x==1))
  x=rnorm(n,mu[ind],sqrt(v[ind]))
  return(list(x=x,ind=ind))
}

#' Generates a gaussian multivariate mixture
#'
#' Generates a gaussian multivariate mixture
#' @param p Needs to a vector of probabilities
#' @param mu Needs to be a list of vector
#' @param sigma Needs to be a list of matrix
#' @param seed Fix a random seed (optional)
#' @export
gen_mixt_m=function(n,p,mu,sigma,seed)
{
  if(missing(seed)==F) set.seed(seed)
  dim=length(p)
  ind=rmultinom(n,1,prob=p)
  ind=apply(ind,2,function(x) which(x==1))
  t=factor(c(1:dim))
  tavola=table(ind)
  na=matrix(as.numeric(t[t %in% names(tavola)]),ncol=1)
  if(nrow(na)==1)
  {
    tab=rmvnorm(tavola[1],mu[[na]],sigma[[na]])
    ind1=(unlist(sapply(as.numeric(names(tavola)),function(x) rep(x,times=tavola[x]),simplify = T)))
    return(list(x=tab,ind=ind1))
  } else {
    tab=apply(na,1,function(x) rmvnorm(tavola[x],mu[[x]],sigma[[x]]))
    if(class(tab)=="list")
    {
      ind1=(unlist(sapply(as.numeric(names(tavola)),function(x) rep(x,times=tavola[x]),simplify = T)))
      return(list(x=do.call(rbind,tab),ind=ind1))} else {
        ind1=(unlist(sapply(as.numeric(names(tavola)),function(x) rep(x,times=tavola[x]),simplify = T)))
        return(list(x=tab,ind=ind1))}
    }
}

#' Build multivariate/univariate binned data using a square grid
#'
#' Build multivariate/univariate binned data using a square grid
#' @param x Needs to a matrix or a vector of data
#' @param R Needs to be a number. It indicates the density of grid along each dimension
#' @param sigma Needs to be a list of matrix
#' @return The binned data and the square grid
#' @export
build_bin_mext=function (x, R)
{
  d = NCOL(x)
  if (d == 1) {
    if (R == 1) {
      b = c(-Inf, mean(x), Inf)
      n = numeric(R + 1)
      n[1] = sum(x <= b[2])
      n[2] = length(x) - sum(n)
    }
    else {
      b = c(-Inf, seq(min(x), max(x),
                      length = R), Inf)
      n = numeric(R + 1)
      n[1] = sum(x <= b[2])
      for (i in 2:R) n[i] = sum(x <= b[i + 1] & x > b[i])
      n[R + 1] = length(x) - sum(n)
    }
    b = cbind(b[1:(R + 1)], b[2:(R + 2)])
  }
  else {
    if (R == 1) {
      b = apply(x, 2, function(y) mean(y))
      n = numeric((R + 1)^d)
      index = matrix(0, nrow(x), d)
      for (i in 1:d) index[, i] = sapply(x[, i], function(x) sum(x >=
                                                                   b[i]))
    }
    else {
      b = apply(x, 2, function(y) seq(min(y),
                                      max(y), length = R))
      n = numeric((R + 1)^d)
      index = matrix(0, nrow(x), d)
      for (i in 1:d) index[, i] = sapply(x[, i], function(x) sum(x >=
                                                                   b[, i]))
    }
    index = as.numeric((index[, sort(1:d, decreasing = T)]) %*%
                         (R + 1)^c(0:(d - 1))) + 1
    names(n) = c(1:((R + 1)^d))
    n[as.numeric(names(table(index)))] = table(index)
    b = rbind(rep(-Inf, d), b, rep(Inf, d))
    ext = list(as.data.frame(matrix(b[1:(R + 1), ], ncol = d)),
               as.data.frame(matrix(b[2:(R + 2), ], ncol = d)))
    b = cbind(expand.grid(rev(ext[[1]]))[sort(c(1:d), decreasing = T)],
              expand.grid(rev(ext[[2]]))[sort(c(1:d), decreasing = T)])
  }
  return(list(bin = n, grid = as.matrix(b)))
}

#' Binned-EM algorithm for univariate Gaussian mixtures
#'
#' Binned-EM algorithm for univariate Gaussian mixtures
#' @param bin Univariate binned data
#' @param pi0 Starting point for proportions. If missing, a binned_EM algorithm for a single Gaussian is performed
#' @param mu0 Starting point for means (numerical vector)
#' @param sigma0 Starting point for variances (numerical vector)
#' @param tim Maximum number of iterations
#' @param tol Tolerance of the algorithm
#' @return List with final estimates, log-likelihood sequence, parameters sequence, BIC and medium time per iteration
#' @export
em_univ=function (bin, pi0, m0, s0,tim,tol)
{
  if(missing(pi0))
  {
    t0=Sys.time()
    bin1 = nozero(bin)
    n=sum(bin1$bin)
    ext1 = bin1$grid
    logli=numeric(tim)
    parm=matrix(0,tim+1,2)
    parm[1,]=c(m0,s0)
    flag=1
    j=1
    logli[j]=logli_mult(bin1,1,m0,s0)
    while(flag==1)
    {
      rip1=apply(ext1,1,function(x) (pnorm(x[2],mean=m0,sqrt(s0))-pnorm(x[1],mean=m0,sqrt(s0))))
      nor1=apply(ext1,1,function(x) (dnorm(x[2],mean=m0,sqrt(s0))-dnorm(x[1],mean=m0,sqrt(s0))))
      c1=apply(ext1,1,function(x)
      {
        if(x[1]==-Inf)
          return(x[2]*dnorm(x[2],mean=m0,sqrt(s0)))
        if(x[2]==Inf)
          return(-x[1]*dnorm(x[1],mean=m0,sqrt(s0)))
        if(x[1]!=-Inf & x[2]!=Inf)
          return(x[2]*dnorm(x[2],mean=m0,sqrt(s0))-x[1]*dnorm(x[1],mean=m0,sqrt(s0)))
      }
      )
      mu1=sum(bin1$bin*(m0-s0*nor1/rip1))/n
      sigma1=sum(bin1$bin*((s0*(rip1+(2*mu1-m0)*nor1-c1)+rip1*(mu1-m0)^2)/rip1))/n
      j=j+1
      parm[j,]=c(mu1,sigma1)
      logli[j]=logli_mult(bin1,1,mu1,sigma1)
      if(abs((logli[j]-logli[j-1])/logli[j-1])>tol &
         j <= tim)
      {
        s0=v1
        m0=mu1
        flag=1
      }
    }
    conv = TRUE
    if (j > tim)
      conv = FALSE
    logli = logli[1:min(j, tim + 1)]
    parm=parm[1:min(j, tim + 1),]
    t1=Sys.time()
    return(list(mu = mu1, v = sigma1,logli=lo,parm=parm,bic=bic,timem=as.numeric(difftime(t1,t0,units="sec"))/(length(logli)-1)))} else {
      t0=Sys.time()
      class = length(pi0)
      bin1 = nozero(bin)
      ext1 = bin1$grid
      n=sum(bin1$bin)
      logli=numeric(tim)
      parm=matrix(0,tim+1,3*class)
      parm[1,]=c(pi0,m0,s0)
      flag=1
      j=1
      logli[j]=logli_mult(bin1,pi0,m0,s0)

      while(flag==1)
      {
        mat=cbind(m0,s0)
        s = apply(ext1, 1, function(x) apply(mat, 1, function(y) {
          max(pnorm(x[2],y[1], sqrt(y[2])) - pnorm(x[1], y[1], sqrt(y[2])),-pnorm(x[2],y[1], sqrt(y[2]),lower.tail=F)+ pnorm(x[1], y[1], sqrt(y[2]),lower.tail=F))
        }))
        s1 = as.numeric(pi0 %*% s)
        sm = apply(ext1, 1, function(x) apply(mat, 1, function(y) (dnorm(x[2], y[1], sqrt(y[2])) - dnorm(x[1], y[1], sqrt(y[2])))))
        sv = apply(ext1, 1, function(x) apply(mat, 1, function(y) {
          if (x[1] == -Inf)
            return(x[2] * dnorm(x[2], mean = y[1], sqrt(y[2])))
          if (x[2] == Inf)
            return(-x[1] * dnorm(x[1], mean = y[1], sqrt(y[2])))
          if (x[1] != -Inf & x[2] != Inf)
            return(x[2] * dnorm(x[2], mean = y[1], sqrt(y[2])) -
                     x[1] * dnorm(x[1], mean = y[1], sqrt(y[2])))
        }))
        int = t(t(s * pi0)/s1)
        pp = NULL
        for (i in 1:class) pp = rbind(pp, s1)
        int3 = (1/(pp/pi0)) * (s * mat[, 1] - sm * mat[, 2])
        sum1 = apply(int, 1, function(x) sum(bin1$bin * x))
        pi1=sum1/sum(bin1$bin)
        sum2 = apply(int3, 1, function(x) sum(bin1$bin * x))
        mu1 = sum2/sum1
        int4 = ((1/(pp/pi0)) * ((s + sm * (2 * mu1 - mat[, 1]) -
                                   sv) * mat[, 2] + s * (mu1 - mat[, 1])^2))
        sum3 = apply(int4, 1, function(x) sum(bin1$bin * x))
        v1 = sum3/sum1
        j=j+1
        parm[j,]=c(pi1,mu1,v1)
        logli[j]=logli_mult(bin1,pi1,mu1,v1)
        flag=0
        if(abs((logli[j]-logli[j-1])/logli[j-1])>tol &
           j <= tim)
        {
          s0=v1
          m0=mu1
          pi0=pi1
          flag=1
        }
      }
      conv = TRUE
      if (j > tim)
        conv = FALSE
      logli = logli[1:min(j, tim + 1)]
      parm=parm[1:min(j, tim + 1),]
      t1=Sys.time()
      return(list(bin=bin1,mu = mu1, v = v1, pi1 = pi1,logli=logli,parm=parm,timem=as.numeric(difftime(t1,t0,units="sec"))/(length(logli)-1)))}
}

d_mixt=function(x,pi,mu,v)
{
  sum(pi*dnorm(x,mu,sqrt(v)))
}

#' Density of univariate Gaussian mixture
#'
#' Density of univariate Gaussian mixture
#' @param x real vector
#' @param pi proportions
#' @param mu means
#' @param v variances
#' @return Density of univariate Gaussian mixture
#' @export
d_mixt=Vectorize(d_mixt,"x")


#' Density of multivariate Gaussian mixture
#'
#' Density of multivariate Gaussian mixture
#' @param x real matrix
#' @param pi proportions
#' @param mu means
#' @param v variances
#' @return Density of multivariate Gaussian mixture
#' @export
d_mixt_m=function(x,pi,mu,v)
{
  as.numeric(sapply(1:length(pi),function(y)dmvnorm(x,mean=mu[[y]],sigma=v[[y]]))%*%pi)
}

cumul_mixt=function(x,pi,mu,v,log=T,lower.tail=TRUE)
{
  if(log==F)
  sum(pi*(pnorm(x,mu,sqrt(v),lower.tail = lower.tail))) else
  {
    a=log(pi)+pnorm(x,mu,sqrt(v),log.p=T,lower.tail = lower.tail)
    l=which.max(a)
    if(min(a)==-Inf) return(-Inf)
    a[l]+log(1+sum(exp(a[-l]-a[l])))
  }
}

#' Univariate Gaussian mixture cdf
#'
#' Univariate Gaussian mixture cdf
#' @param x real vector
#' @param pi proportions
#' @param mu means
#' @param v variances
#' @param log Logical. If TRUE (default), it returns the logarithm
#' @param lower.tails Logical. If TRUE (default), it returns P[X<=x], if FALSE P[X>=x]
#' @return Univariate Gaussian mixture cdf
#' @export
cumul_mixt=Vectorize(cumul_mixt,"x")

#' Univariate Gaussian mixture log-likelihood
#'
#' Univariate Gaussian mixture log-likelihood
#' @param x real vector
#' @param pi proportions
#' @param mu means
#' @param v variances
#' @return Univariate Gaussian mixture log-likelihood
#' @export
logli_mixt=function(x,pi,mu,v)
{
  if(NCOL(x)==1)
  sum(log(d_mixt(x,pi,mu,v))) else
    sum(log(d_mixt_m(x,pi,mu,v)))
}

#' Multinomial log-likelihood for binned univariate Gaussian mixture
#'
#' Multinomial log-likelihood for binned univariate Gaussian mixture
#' @param bin Binned data
#' @param pi proportions
#' @param mu means
#' @param v variances
#' @return Multinomial log-likelihood for binned univariate Gaussian mixture
#' @export
logli_mult=function (bin, pi, mu, v)
{
  R = length(bin$bin)
  if (R == 2) {
    return(bin$bin[1] * cumul_mixt(bin$grid[1,2],pi,mu,v,log=T,lower.tail = T) + bin$bin[R] * cumul_mixt(bin$grid[1,2],pi,mu,v,log=T,lower.tail = F))
  } else {
    a=apply(bin$grid, 1, function(x) logdiff_mixt(x,pi, mu, v,lower.tail = T))
    a[a==-Inf | is.na(a) | is.nan(a)]=apply(matrix(bin$grid[a==-Inf | is.na(a) | is.nan(a),],ncol=2), 1, function(x) logdiff_mixt(x,pi, mu, v,lower.tail = F))
    return(sum(bin$bin *a))
  }}
#logli_mult=function (bin, pi, mu, v)
#{
 # R = length(bin$bin)
  #if (R == 2) {
   # return(bin$bin[1] * cumul_mixt(bin$grid[1,2],pi,mu,v,log=T,lower.tail = T) + bin$bin[R] * cumul_mixt(bin$grid[1,2],pi,mu,v,log=T,lower.tail = F))
  #} else {
    #a=apply(bin$grid, 1, function(x) {
    #  cumul_mixt(x[2],pi, mu, v,lower.tail = T)+log(1-exp((cumul_mixt(x[1],pi, mu, v,lower.tail = T) - cumul_mixt(x[2], pi, mu, v,lower.tail = T))))
    #})
    #a[a==-Inf]=apply(matrix(bin$grid[a==-Inf,],ncol=2), 1, function(x) cumul_mixt(x[1],pi, mu, v,lower.tail = F)+log(1-exp((cumul_mixt(x[2],pi, mu, v,lower.tail = F) - cumul_mixt(x[1], pi, mu, v,lower.tail = F)))))
    #return(sum(bin$bin *a))
  #}}

#' EM algorithm for univariate Gaussian mixture
#'
#' EM algorithm for univariate Gaussian mixture
#' @param x  Data
#' @param pi0 Starting point for proportions
#' @param mu0 Starting point for means
#' @param v0 Starting point for variances
#' @param tim Maximum number of iterations
#' @param tol Tolerance of the algorithm
#' @return List with final estimates, log-likelihood sequence, BIC and medium time per iteration
#' @export
em_subs=function(x,pi0,mu0,v0,tol,tim)
{
  time0=Sys.time()
  flag=1
  j=1
  m=length(x)
  lo=numeric(tim+1)
  lo[1]=logli_mixt(x,pi0,mu0,v0)
  while(flag==1)
  {
    par=cbind(pi0,mu0,v0)
    t=apply(par,1,function(y) y[1]*dnorm(x,y[2],sqrt(y[3]))/d_mixt(x,pi0,mu0,v0))
    sum1=apply(t,2,sum)
    pi1=sum1/m
    mu1=(apply(t,2,function(y) sum(y*x)))/sum1
    z=sapply(mu1,FUN=function(y) (x-y)^2)
    v1=diag(t(t)%*%z)/sum1
    flag=0
    print(paste(paste("Iter ",j),paste(c(pi1, mu1,v1),collapse=" "),sep=": "))
    lo[j+1]=logli_mixt(x,pi1,mu1,v1)
    if(abs(lo[j+1]-lo[j])>tol)
    {
      flag=1
      pi0=pi1
      mu0=mu1
      v0=v1
      j=j+1
    }
    if(j>tim) flag=0
  }

  conv=TRUE
  if(j>tim) conv=FALSE
  lo=lo[2:j]
  bic=max(lo)-(2*length(pi0)+length(pi0)-1)*log(m)/2
  return(list(sublength=m,pi=pi1,mu=mu1,v=v1,convergence=conv,logseq=lo,bic=bic,timem=as.numeric(difftime(Sys.time(),time0,units="sec"))/length(lo)))
}

tab1=function(bin1,bin2)
{
  a=sample(c(1,2),bin1$bin[1],replace=T)
  i1=c(length(a[a==1]),length(a[a==2]))
  b=sample(c(1,2),bin1$bin[2],replace=T)
  i2=c(length(b[b==1]),length(b[b==2]))
  if(any(i1+i2-bin2$bin!=c(0,0))==T)
    return((c(i1,i2))) else
      return(NULL)
}
tab2=function(bin1,bin2)
{
  a=sample(c(1,2),bin2$bin[1],replace=T)
  i1=c(length(a[a==1]),length(a[a==2]))
  b=sample(c(1,2),bin2$bin[2],replace=T)
  i2=c(length(b[b==1]),length(b[b==2]))
  if(any(i1+i2-bin1$bin!=c(0,0))==T)
    return((c(i1[1],i2[1],i1[2],i2[2]))) else
      return(NULL)
}
fatt=function(n)
{
  if(n==0) {return(1)} else {
    return( n*fatt(n-1))
  }
}


#' Marginal composite log-likelihood of binned diagonal Gaussian mixture
#'
#' Marginal composite log-likelihood of binned diagonal Gaussian mixture
#' @param bin Binned data
#' @param pi0 Proportions
#' @param mu0 List of means
#' @param sigma0 List of diagonal matrices
#' @return Marginal composite log-likelihood for binned diagonal Gaussian mixture
#' @export
logli_marg=function(bin,pi0,mu0,sigma0)
{
  d=length(mu0[[1]])
  mu0=sapply(1:d,function(y) as.numeric(lapply(mu0,function(x) x[y])),simplify=F)
  sigma0=sapply(1:d,function(y) as.numeric(lapply(sigma0,function(x) x[y,y])),simplify=F)
  return(
    sum(sapply(1:d,function(x) logli_mult(bin[[x]],pi0,mu0[[x]],sigma0[[x]]))))
}
nozero=function(bin1)
{
  ind1=bin1$bin>0
  bin1$bin=bin1$bin[ind1]
  bin1$grid=bin1$grid[ind1,]
  return(bin1)
}
em_dim=function(bin,pi0,mu0,sigma0,d)
{
  class=length(pi0)
  bin1=bin[[d]]
  ext1=bin1$grid
  mat=cbind(as.numeric(lapply(mu0,function(x) x[d])),as.numeric(lapply(sigma0,function(x) x[d,d])))
  s=apply(ext1, 1, function(x) apply(mat, 1, function(y) {
    max(pnorm(x[2],y[1], sqrt(y[2])) - pnorm(x[1], y[1], sqrt(y[2])),-pnorm(x[2],y[1], sqrt(y[2]),lower.tail=F)+ pnorm(x[1], y[1], sqrt(y[2]),lower.tail=F))
  }))
  s1=as.numeric(pi0%*%s)
  sm=apply(ext1,1,function(x) apply(mat,1,function(y) (dnorm(x[2],y[1],sqrt(y[2]))-dnorm(x[1],y[1],sqrt(y[2])))))
  sv=apply(ext1,1,function(x) apply(mat,1,function(y)
  {
    if(x[1]==-Inf)
      return(x[2]*dnorm(x[2],mean=y[1],sqrt(y[2])))
    if(x[2]==Inf)
      return(-x[1]*dnorm(x[1],mean=y[1],sqrt(y[2])))
    if(x[1]!=-Inf & x[2]!=Inf)
      return(x[2]*dnorm(x[2],mean=y[1],sqrt(y[2]))-x[1]*dnorm(x[1],mean=y[1],sqrt(y[2])))
  }
  ))
  int=t(t(s*pi0)/s1)
  pp=NULL
  for(i in 1:class)
    pp=rbind(pp,s1)
  int3=(1/(pp/pi0))*(s*mat[,1]-sm*mat[,2])
  sum1=apply(int,1,function(x) sum(bin1$bin*x))
  sum2=apply(int3,1,function(x) sum(bin1$bin*x))
  mu1=sum2/sum1
  int4=((1/(pp/pi0))*((s+sm*(2*mu1-mat[,1])-sv)*mat[,2]+s*(mu1-mat[,1])^2))
  sum3=apply(int4,1,function(x) sum(bin1$bin*x))
  v1=sum3/sum1
  return(list(mu=mu1,v=v1,sum1=sum1))
}

em_dim_nomixt=function(bin,mu0,sigma0,d)
{
  bin1=bin[[d]]
  ext1=bin1$grid
  mu0=mu0[d]
  sigma0=sigma0[d,d]

  rip1=apply(ext1,1,function(x) (pnorm(x[2],mean=mu0,sqrt(sigma0))-pnorm(x[1],mean=mu0,sqrt(sigma0))))
  nor1=apply(ext1,1,function(x) (dnorm(x[2],mean=mu0,sqrt(sigma0))-dnorm(x[1],mean=mu0,sqrt(sigma0))))
  c1=apply(ext1,1,function(x)
  {
    if(x[1]==-Inf)
      return(x[2]*dnorm(x[2],mean=mu0,sqrt(sigma0)))
    if(x[2]==Inf)
      return(-x[1]*dnorm(x[1],mean=mu0,sqrt(sigma0)))
    if(x[1]!=-Inf & x[2]!=Inf)
      return(x[2]*dnorm(x[2],mean=mu0,sqrt(sigma0))-x[1]*dnorm(x[1],mean=mu0,sqrt(sigma0)))
  }
  )

  mu1=sum(bin1$bin*(mu0-sigma0*nor1/rip1))/n
  sigma1=sum(bin1$bin*((sigma0*(rip1+(2*mu1-mu0)*nor1-c1)+rip1*(mu1-mu0)^2)/rip1))/n
  return(list(mu=mu1,sigma=sigma1))
}

#' Binned CL-EM algorithm for multivariate binned diagonal Gaussian mixture
#'
#' Binned CL-EM algorithm for multivariate binned diagonal Gaussian mixture
#' @param bin  Binned data
#' @param pi0 Starting point for proportions
#' @param mu0 Starting point for means
#' @param sigma0 Starting point for variances
#' @param tol Tolerance of the algorithm
#' @param tim Maximum number of iterations
#' @return List with original binned data, final estimates, log-likelihood and parameters sequences and medium time per iteration
#' @export
em_bin_gaus_comp=function(bin,pi0,mu0,sigma0,tol,tim,print=F)
{
  t0=Sys.time()
  flag=1
  j=1
  R=as.numeric(lapply(bin,function(x) length(x$bin)-1))
  n=sum(bin[[1]]$bin)
  bin=lapply(bin,function(x) nozero(x))
  class=length(pi0)
  d=length(mu0[[1]])
  kl=numeric(tim)
  parm=matrix(0,tim+1,2*d*class+class)
  logli=numeric(tim+1)
  logli[j]=logli_marg(bin,pi0,mu0,sigma0)
  parm[j,]=c(pi0,unlist(mu0),unlist(sigma0)[unlist(sigma0)!=0])
  while(flag==1)
  {
    if(print==T) cat(j," ")
    dim1=sapply(1:d,function(x) em_dim(bin,pi0,mu0,sigma0,x),simplify = F)
    pi1=colSums(matrix(unlist(lapply(dim1,function(x) x$sum)),byrow=T,ncol=class))/(d*n)

    mu1=list()
    sigma1=list()
    for(i in 1:class)
    {
      mu1[[i]]=as.numeric(lapply(dim1,function(x) x$mu[[i]]))
      sigma1[[i]]=diag(as.numeric(lapply(dim1,function(x) x$v[[i]])))
    }
    j=j+1
    logli[j]=logli_marg(bin,pi1,mu1,sigma1)
    mu1_ul=unlist(mu1)
    sigma1_ul=unlist(sigma1)
    sigma1_ul=sigma1_ul[sigma1_ul!=0]

    parm[j,]=c(pi1,mu1_ul,sigma1_ul)
    flag=0
    #print(paste(paste("Iter ",j),paste(c(pi1,mu1, sigma1_ul),collapse=" "),sep=": "))
    if(abs((logli[j]-logli[j-1])/logli[j-1])>tol & j<=tim)
    {
      flag=1
      pi0=pi1
      mu0=mu1
      sigma0=sigma1
    }
  }
  conv=TRUE
  if(j>tim) conv=FALSE
  logli=logli[1:min(j,tim+1)]
  return(list(bin=bin,density=R,pi=pi1,mu=mu1,v=sigma1,logli=logli,parm=parm[1:min(j,tim+1),],convergence=conv,timem=as.numeric(difftime(Sys.time(),t0,units="sec"))/(length(logli)-1)))
}
#' Binned CL-EM algorithm for multivariate binned diagonal Gaussian model
#'
#' Binned CL-EM algorithm for multivariate binned diagonal Gaussian model
#' @param bin Binned data
#' @param mu0 Starting point for means
#' @param sigma0 Starting point for variances
#' @param tol Tolerance of the algorithm
#' @param tim Maximum number of iterations
#' @return List with original binned data, final estimates, log-likelihood and parameters sequences and medium time per iteration
#' @export
em_bin_gaus_comp_nomixt=function(bin,mu0,sigma0,tol,tim)
{
  t0=Sys.time()
  flag=1
  j=1
  R=as.numeric(lapply(bin,function(x) length(x$bin)-1))
  n=sum(bin[[1]]$bin)
  bin=lapply(bin,function(x) nozero(x))
  d=length(mu0)
  lo=numeric(tim+1)
  parm=matrix(0,tim+1,d*2)
  parm[j,]=c(mu0,diag(sigma0))
  lo[j]=logli_marg(bin,1,list(mu0),list(sigma0))
  while(flag==1)
  {
    dim_nomixt=sapply(1:d,function(x) em_dim_nomixt(bin,mu0,sigma0,x),simplify=F)
    mu1=as.numeric(lapply(dim_nomixt,function(x) x$mu))
    sigma1=diag(as.numeric(lapply(dim_nomixt,function(x) x$sigma)))
    #print("ok")
    flag=0
    #print(paste(paste("Iter ",j),paste(c(mu1, sigma1),collapse=" "),sep=": "))
    j=j+1
    parm[j,]=c(mu1,diag(sigma1))
    lo[j]=logli_marg(bin,1,list(mu1),list(sigma1))
    if(abs((lo[j]-lo[j-1])/lo[j-1])>tol & j<=tim)
    {
      flag=1
      mu0=mu1
      sigma0=sigma1
    }
  }
  conv=TRUE
  if(j>tim) conv=FALSE
  lo=lo[1:min(j,tim+1)]
  return(list(bin=bin,density=R,mu=mu1,v=sigma1,convergence=conv,logli=lo,parm=parm[1:min(j,tim+1),],timem=as.numeric(difftime(Sys.time(),t0,units="sec"))/(length(lo)-1)))
}
simulazione=function(dati,maxit,nrep,R,tol,class)
{
  d=ncol(dati)
  R=as.matrix(R)
  n=nrow(dati)
  ret=list()
  for(cl in class)
  {
    if(cl==1)
    {
      temp=list()
      for(i in 1:nrep)
      {
        mu0=sapply(1:d,function(x) runif(1,min(dati[,x]),max(dati[,x])))
        sigma0=diag(d)
        temp[[i]]=list()
        for(j in 1:ncol(R))
        {
          bin=sapply(1:d,function(x) build_bin_mext(dati[,x],R[x,j]),simplify=F)
          temp[[i]][[j]]=try(em_bin_gaus_comp_nomixt(bin,mu0,sigma0,tol,maxit),silent=T)
        }
      }
    } else {
      temp=list()
      for(i in 1:nrep)
      {
        pi0=runif(cl)
        pi0=pi0/sum(pi0)
        mu0=list()
        for(j in 1:cl)
          mu0[[j]]=sapply(1:d,function(x) runif(1,min(dati[,x]),max(dati[,x])))
        sigma0=list()
        for(j in 1:cl)
          sigma0[[j]]=diag(d)
        temp[[i]]=list()
        for(j in 1:ncol(R))
        {
          bin=sapply(1:d,function(x) build_bin_mext(dati[,x],R[x,j]),simplify=F)
          temp[[i]][[j]]=try(em_bin_gaus_comp(bin,pi0,mu0,sigma0,tol,maxit),silent=T)
        }
      }}
    maxi=matrix(-Inf,nrep,ncol(R))
    for(i in 1:nrep)
      for(j in 1:ncol(R))
      {
        if(class(temp[[i]][[j]])=="list") maxi[i,j]=max(temp[[i]][[j]]$logli)
      }
    max_temp=list()
    for(i in 1:ncol(R))
    {
      max_temp[[i]]=temp[[which.max(maxi[,i])]][[i]]
      if(class(max_temp[[i]])=="list")
        max_temp[[i]]$bic=-2*max(max_temp[[i]]$logli)+(ncol(max_temp[[i]]$parm)-1)*log(n)
    }
    ret[[which(class==cl)]]=max_temp
  }
  names(ret)=paste("class",class,sep="")
  return(ret)
}
#' Build multivariate binned data using a square grid
#'
#' Build multivariate binned data using a square grid
#' @param num_sim Number of simulation
#' @param maxit Maximum number of iterations EM
#' @param tol Tolerance of EM
#' @param class Number of groups of the mixture models
#' @param nrep Number of different initialization of EM
#' @param R Needs to be a matrix or a numerical vector. The element (i,j) indicates the density of the j-th grid along the i-th dimension
#' @param n number of instances simulated
#' @param pi Proportions of the true mixture (vector)
#' @param mu Means of the true mixture (list of vectors)
#' @param sigma Variances of the true mixture (list of diagonal matrices)
#' @return Binned CL-EM estimation
#' @export
sim_mult=function(num_sim,maxit,tol,class,nrep,R,n,pi,mu,sigma)
{
  set.seed(123)
  simul=list()
  pb <- txtProgressBar(min = 0, max = num_sim, style = 3)
  for(i in 1:num_sim)
  {
    x=gen_mixt_m(n,pi,mu,sigma)
    dati=x$x
    simul[[i]]=simulazione(dati,maxit,nrep,R,tol,class)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(list(n=n,simul=simul,density=R,class=class,num_sim=num_sim,tol=tol,maxit=maxit,pi_true=pi,mu_true=mu,sigma_true=sigma))
}

#' Binned CL-EM algorithm
#'
#' Binned CL-EM algorithm
#' @param data Needs to a matrix of data
#' @param R Needs to be a matrix or a numerical vector. The element (i,j) indicates the density of the j-th grid along the i-th dimension
#' @param class Number of groups of the mixture models
#' @param maxit Maximum number of iterations EM
#' @param tol Tolerance of EM
#' @param nrep Number of different initialization of EM
#' @return Binned CL-EM algorithm output
#' @export
binmixt=function (data, class, R, maxit, tol, nrep,print=T)
{
  if(missing(print)) print=T
  n = NROW(data)
  R = as.matrix(R)
  d=NCOL(data)
  set.seed(123)
  if(print==T)
  {
  pb <- txtProgressBar(min = 0, max = length(class)*nrep,
                       style = 3)
  if(d==1)
  {
    bin=list()
    for(j in 1:nrow(R))
      bin[[j]] = build_bin_mext(data, R[j,])
    ret = list()
    for (cl in class) {
      if (cl == 1) {
        temp = list()
        for (i in 1:nrep) {
          mu0 = sapply(1:d, function(x) runif(1, min(data[,
                                                          x]), max(data[, x])))
          sigma0 = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                    x])))))
          temp[[i]] = list()
          for (j in 1:nrow(R)) {
            temp[[i]][[j]] = try(em_univ(bin[[j]],m0=mu0, s0=sigma0,tim= maxit,tol=tol), silent = T)
          }
          value=(cl-min(class))*nrep+i
          setTxtProgressBar(pb,value)
        }
      }
      else {
        temp = list()
        for (i in 1:nrep) {
          pi0 = runif(cl)
          pi0 = pi0/sum(pi0)
          mu0 = runif(cl,min(data), max(data))
          sigma0 =runif(cl,0,(var(data)))
          temp[[i]] = list()
          for (j in 1:nrow(R)) {
            temp[[i]][[j]] = try(em_univ(bin[[j]],pi0, mu0, sigma0, maxit,tol), silent = T)
          }
          value=(cl-min(class))*nrep+i
          setTxtProgressBar(pb,value)
        }

      }
      maxi = matrix(-Inf, nrep, nrow(R))
      for (i in 1:nrep) for (j in 1:nrow(R)) {
        if (class(temp[[i]][[j]]) == "list")
          maxi[i, j] = max(temp[[i]][[j]]$logli)
      }
      max_temp = list()
      for (i in 1:nrow(R)) {
        max_temp[[i]] = temp[[which.max(maxi[, i])]][[i]]
        if (class(max_temp[[i]]) == "list")
          max_temp[[i]]$bic = -2 * max(max_temp[[i]]$logli) +
            (ncol(max_temp[[i]]$parm) - 1) * log(n)
      }
      ret[[which(class == cl)]] = list(best=max_temp,all=temp)
    }
  } else {
    bin=list()
    for(j in 1:ncol(R))
      bin[[j]] = sapply(1:d, function(x) build_bin_mext(data[,
                                                             x], R[x, j]), simplify = F)


    ret = list()
    for (cl in class) {
      if (cl == 1) {
        temp = list()
        for (i in 1:nrep) {
          mu0 = sapply(1:d, function(x) runif(1, min(data[,
                                                          x]), max(data[, x])))
          sigma0 = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                    x])))))
          temp[[i]] = list()
          for (j in 1:ncol(R)) {
            temp[[i]][[j]] = try(em_bin_gaus_comp_nomixt(bin[[j]],
                                                         mu0, sigma0, tol, maxit), silent = T)
          }
          value=(cl-min(class))*nrep+i
          setTxtProgressBar(pb,value)
        }
      }
      else {
        temp = list()
        for (i in 1:nrep) {
          pi0 = runif(cl)
          pi0 = pi0/sum(pi0)
          mu0 = list()
          for (j in 1:cl) mu0[[j]] = sapply(1:d, function(x) runif(1,
                                                                   min(data[, x]), max(data[, x])))
          sigma0 = list()
          for (j in 1:cl) sigma0[[j]] = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                                         x])))))
          temp[[i]] = list()
          for (j in 1:ncol(R)) {
            temp[[i]][[j]] = try(em_bin_gaus_comp(bin[[j]],
                                                  pi0, mu0, sigma0, tol, maxit), silent = T)
          }
          value=(cl-min(class))*nrep+i
          setTxtProgressBar(pb,value)
        }

      }
      maxi = matrix(-Inf, nrep, ncol(R))
      for (i in 1:nrep) for (j in 1:ncol(R)) {
        if (class(temp[[i]][[j]]) == "list")
          maxi[i, j] = max(temp[[i]][[j]]$logli)
      }
      max_temp = list()
      for (i in 1:ncol(R)) {
        max_temp[[i]] = temp[[which.max(maxi[, i])]][[i]]
        if (class(max_temp[[i]]) == "list")
          max_temp[[i]]$bic = -2 * max(max_temp[[i]]$logli) +
            (ncol(max_temp[[i]]$parm) - 1) * log(n)
      }
      ret[[which(class == cl)]] = list(best=max_temp,all=temp)
    } }
  close(pb)
  } else {
  if(d==1)
  {
    bin=list()
    for(j in 1:nrow(R))
      bin[[j]] = build_bin_mext(data, R[j,])
    ret = list()
    for (cl in class) {
      if (cl == 1) {
        temp = list()
        for (i in 1:nrep) {
          mu0 = sapply(1:d, function(x) runif(1, min(data[,
                                                          x]), max(data[, x])))
          sigma0 = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                    x])))))
          temp[[i]] = list()
          for (j in 1:nrow(R)) {
            temp[[i]][[j]] = try(em_univ(bin[[j]],m0=mu0, s0=sigma0,tim= maxit,tol=tol), silent = T)
          }

        }
      }
      else {
        temp = list()
        for (i in 1:nrep) {
          pi0 = runif(cl)
          pi0 = pi0/sum(pi0)
          mu0 = runif(cl,min(data), max(data))
          sigma0 =runif(cl,0,(var(data)))
          temp[[i]] = list()
          for (j in 1:nrow(R)) {
            temp[[i]][[j]] = try(em_univ(bin[[j]],pi0, mu0, sigma0, maxit,tol), silent = T)
          }
          value=(cl-min(class))*nrep+i
        }

      }
      maxi = matrix(-Inf, nrep, nrow(R))
      for (i in 1:nrep) for (j in 1:nrow(R)) {
        if (class(temp[[i]][[j]]) == "list")
          maxi[i, j] = max(temp[[i]][[j]]$logli)
      }
      max_temp = list()
      for (i in 1:nrow(R)) {
        max_temp[[i]] = temp[[which.max(maxi[, i])]][[i]]
        if (class(max_temp[[i]]) == "list")
          max_temp[[i]]$bic = -2 * max(max_temp[[i]]$logli) +
            (ncol(max_temp[[i]]$parm) - 1) * log(n)
      }
      ret[[which(class == cl)]] = list(best=max_temp,all=temp)
    }
  } else {
  bin=list()
  for(j in 1:ncol(R))
    bin[[j]] = sapply(1:d, function(x) build_bin_mext(data[,
                                                           x], R[x, j]), simplify = F)


  ret = list()
  for (cl in class) {
    if (cl == 1) {
      temp = list()
      for (i in 1:nrep) {
        mu0 = sapply(1:d, function(x) runif(1, min(data[,
                                                        x]), max(data[, x])))
        sigma0 = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                  x])))))
        temp[[i]] = list()
        for (j in 1:ncol(R)) {
          temp[[i]][[j]] = try(em_bin_gaus_comp_nomixt(bin[[j]],
                                                       mu0, sigma0, tol, maxit), silent = T)
        }
      }
    }
    else {
      temp = list()
      for (i in 1:nrep) {
        pi0 = runif(cl)
        pi0 = pi0/sum(pi0)
        mu0 = list()
        for (j in 1:cl) mu0[[j]] = sapply(1:d, function(x) runif(1,
                                                                 min(data[, x]), max(data[, x])))
        sigma0 = list()
        for (j in 1:cl) sigma0[[j]] = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                                       x])))))
        temp[[i]] = list()
        for (j in 1:ncol(R)) {
          temp[[i]][[j]] = try(em_bin_gaus_comp(bin[[j]],
                                                pi0, mu0, sigma0, tol, maxit), silent = T)
        }
      }

    }
    maxi = matrix(-Inf, nrep, ncol(R))
    for (i in 1:nrep) for (j in 1:ncol(R)) {
      if (class(temp[[i]][[j]]) == "list")
        maxi[i, j] = max(temp[[i]][[j]]$logli)
    }
    max_temp = list()
    for (i in 1:ncol(R)) {
      max_temp[[i]] = temp[[which.max(maxi[, i])]][[i]]
      if (class(max_temp[[i]]) == "list")
        max_temp[[i]]$bic = -2 * max(max_temp[[i]]$logli) +
          (ncol(max_temp[[i]]$parm) - 1) * log(n)
    }
    ret[[which(class == cl)]] = list(best=max_temp,all=temp)
  } } }
  names(ret) = paste("class", class, sep = "")
  return(ret)
}

#' BIC computation for a sim_mult output
#'
#' BIC computation for a sim_mult output
#' @param sim A sim_mult output
#' @return Return an array containing BIC value for each simulation
#' @export
bic_sim=function(sim)
{
  len=length(sim$simul)
  nam=length(names(sim$simul[[1]]))
  bic=array(0,dim=c(nam,length(sim$simul[[1]][[1]]),len))
  for(i in 1:len)
  {
    bic1=NULL
    for(j in 1:nam)
    {
      bic1=rbind(bic1,as.numeric(lapply(sim$simul[[i]][[j]],function(x) x$bic)))
    }
    bic[,,i]=bic1
  }
  return(bic)
}

hist_bic=function(bic,class,R)
{
  len=dim(bic)[3]
  er=dim(bic)[2]
  if(er==1)
  {
    mat_bic=numeric(len)
    for(i in 1:len)
      mat_bic[i]=class[which.min(bic[,,i])]
    par(mfrow=c(2,2))
    barplot(prop.table(table(mat_bic)),main=paste("R=",R))
  } else {
    mat_bic=matrix(0,len,er)
    for(i in 1:len)
      mat_bic[i,]=class[(apply(bic[,,i],2,which.min))]
    par(mfrow=c(2,2))
    for(i in 1:ncol(mat_bic))
      barplot(prop.table(table(mat_bic[,i])),main=paste("R=",R[i]))}
}
#' @export
mediumtime=function(sim,index)
{
  m=matrix(unlist((lapply(sim$simul,function(x) (as.numeric(lapply(x[[index]],function(y) y$timem)))))),length(sim$simul),length(sim$simul[[1]][[1]]),byrow=T)
  apply(m,2,mean)
}
#' @export
totaltime=function(sim,index)
{
  m=matrix(unlist((lapply(sim$simul,function(x) (as.numeric(lapply(x[[index]],function(y) y$timem*length(y$logli))))))),length(sim$simul),length(sim$simul[[1]][[1]]),byrow=T)
  apply(m,2,mean)
}
#' @export
klmed=function(sim,pi,mu,sigma,index)
{
  m=matrix(unlist((lapply(sim$simul,function(x) (as.numeric(lapply(x[[index]],function(y) {adaptIntegrate(function(z) d_mixt_m(z,y$pi,y$mu,y$v)*log(d_mixt_m(z,y$pi,y$mu,y$v)/d_mixt_m(z,pi,mu,sigma)),lower=c(-10,-10),upper=c(10,10))$integral
  })))))),length(sim$simul),length(sim$simul[[1]][[1]]),byrow=T)
  apply(m,2,mean)
}
#' @export
lenlikeli=function(sim,index)
{
  m=matrix(unlist((lapply(sim$simul,function(x) (as.numeric(lapply(x[[index]],function(y) length(y$logli))))))),length(sim$simul),length(sim$simul[[1]][[1]]),byrow=T)
  apply(m,2,mean)
}
clc_dim=function(a,d)
{
  mu0=a$mu
  pi0=a$pi
  sigma0=a$v
  bin1=a$bin[[d]]
  mat=cbind(as.numeric(lapply(mu0,function(x) x[d])),as.numeric(lapply(sigma0,function(x) x[d,d])))
  s=apply(bin1$grid,1,function(x) apply(mat,1,function(y) (pnorm(x[2],y[1],sqrt(y[2]))-pnorm(x[1],y[1],sqrt(y[2])))))
  s1=as.numeric(pi0%*%s)
  int=t(t(s*pi0)/s1)
  class=matrix(0,nrow(int),ncol(int))
  for(i in 1:ncol(int))
    class[which.max(int[,i]),i]=1
  int[int==0]=0
  int=int*log(int)
  int[is.nan(int)]=0
  sum1=sum(bin1$bin*colSums(int))
  return(sum1)
}
#' @export
clc=function(a)
{
  if(class(a$mu)!="list")
    return(logli_marg(a$bin,1,list(a$mu),list(a$v))) else {
      d=length(a$mu[[1]])
      a$bin=lapply(a$bin,function(x) nozero(x))
      logl=logli_marg(a$bin,a$pi,a$mu,a$v)
      s=sum(sapply(1:d,function(x) clc_dim(a,x)))
      return(logl+s)}
}
#' @export
clc_sim=function(sim)
{
  len=length(sim$simul)
  nam=length(sim$simul[[1]])
  clc=array(0,dim=c(nam,length(sim$simul[[1]][[1]]),len))
  for(i in 1:len)
  {
    clc1=NULL
    for(j in 1:nam)
    {
      if(is.null(sim$simul[[i]][[j]][[1]]$pi))
      {clc1=rbind(clc1,as.numeric(lapply(sim$simul[[i]][[j]],function(x) max(x$logli))))
      }else
        clc1=rbind(clc1,as.numeric(lapply(sim$simul[[i]][[j]],function(x) clc(x))))
    }
    clc[,,i]=clc1
  }
  return(clc)
}
#' @export
hist_clc=function(clc,class,R)
{
  len=dim(clc)[3]
  er=dim(clc)[2]
  if(er==1)
  {
    mat_clc=numeric(len)
    for(i in 1:len)
      mat_clc[i]=class[which.max(clc[,,i])]
    par(mfrow=c(2,2))
    barplot(prop.table(table(mat_clc)),main=paste("R=",R))
  } else {
    mat_clc=matrix(0,len,er)
    for(i in 1:len)
      mat_clc[i,]=class[(apply(clc[,,i],2,which.max))]
    par(mfrow=c(2,2))
    for(i in 1:ncol(mat_clc))
      barplot(prop.table(table(mat_clc[,i])),main=paste("R=",R[i]))}
}

#' Modified-BIC computation for a sim_mult output
#'
#' Modified-BIC computation for a sim_mult output
#' @param sim A sim_mult output
#' @return Return an array containing Modified-BIC  value for each simulation
#' @export
bic2_sim=function(sim)
{
  len=length(sim$simul)
  nam=length(names(sim$simul[[1]]))
  bic=array(0,dim=c(nam,length(sim$simul[[1]][[1]]),len))
  for(i in 1:len)
  {
    bic1=NULL
    for(j in 1:nam)
    {
      bic1=rbind(bic1,as.numeric(lapply(sim$simul[[i]][[j]],function(x) x$bic+max(x$logli))))
    }
    bic[,,i]=bic1
  }
  return(bic)
}

#' Transform mclust output into a list of parameters compatible with binMixt package
#'
#' Transform mclust output into a list of parameters compatible with binMixt package
#' @param model A mclust model
#' @return List of parameters compatible with binMixt package
#' @export
param_init_mclust=function(model)
{
  G=model$G
  pi_init=model$parameters$pro
  mu_init=list()
  for(i in 1:G)
    mu_init[[i]]=model$parameters$mean[,i]
  s_init=list()
  for(i in 1:G)
    s_init[[i]]=model$parameters$variance$sigma[,,i]
  return(list(pi=pi_init,mu=mu_init,sigma=s_init))
}

#' @export
logbox=function(binmo)
{
  m=matrix(0,length(binmo$all),length(binmo$best))
  for(i in 1:ncol(m))
    m[,i]=as.numeric(lapply(binmo$all,function(x) if(class(x[[i]])=="list") max(x[[i]]$logli) else -Inf))
  return(m)
}

#' build_bin_marg
#'
#' Build marginal binned data
#' @param x Needs to a maxtrix of data
#' @param R Vector. The i-th element indicates the density of grid along the i-th dimension
#' @export
build_bin_marg=function(x,R)
{
  bin=list()
  for(i in 1:length(R))
    bin[[i]]=build_bin_mext(x[,i],R[i])
  return(bin)
}



#' Tool for drawing confidence ellipses for components of a bivariate Gaussian mixture
#'
#' Tool for drawing confidence ellipses for components of a bivariate Gaussian mixture
#' @param s1 Points to draw in the first dimension
#' @param s2 Points to draw in the second dimension
#' @param mu List of means
#' @param v List of variances
#' @param alpha Level of confidence
#' @param color Colour of the ellipses (default=1)
#' @param add Logical value (default=F). If TRUE, ellpses are added to the existing plot
#' @export
contellipse=function(s1,s2,mu,v,alpha,color=1,lty,add=F)
{
  m1=length(s1)
  m2=length(s2)
  gr=expand.grid(s1,s2)
  len=length(mu)
  m=matrix(apply(gr,1,function(x) ellipse(x,mu[[1]],v[[1]],alpha)),m1,m2)
  contour(s1,s2,m,levels=0,col=color,lty=lty,add=add)
  for(i in 2:len)
  {
    m=matrix(apply(gr,1,function(x) ellipse(x,mu[[i]],v[[i]],alpha)),m1,m2)
    contour(s1,s2,m,levels=0,add=T,lty=lty,col=color)}
}

#' Raw posterior probabilities
#'
#' Raw posterior probabilities
#' @param data Matrix of raw data
#' @param model A set of binned estimates
#' @return Matrix of posterior probabilities
#' @export
posterior=function(data,model)
{
  if(NCOL(data)>1)
  {
    mat = t(t(sapply(1:length(model$pi), function(x) dmvnorm(data,
                                                             model$mu[[x]], model$v[[x]]))) * model$pi)
    return(mat)
  } else
  {
    mat = t(t(sapply(1:length(model$pi1), function(x) dnorm(data,
                                                            model$mu[x], sqrt(model$v[x])))) * model$pi1)
    return(mat)
  }
}

#' Binned posterior probabilities  (bivariate case)
#'
#' Binned posterior probabilities  (bivariate case)
#' @param data Matrix of raw data
#' @param model A set of binned estimates
#' @return Matrix of posterior probabilities
#' @export
posterior.bin_dim2=function(bin, model)
{
  mat4 = sapply(1:length(model$pi), function(x) as.numeric(outer((-pnorm(bin[[1]]$grid[,2], model$mu[[x]][1], sqrt(model$v[[x]][1, 1]),lower.tail = F)
                                                                  + pnorm(bin[[1]]$grid[,1], model$mu[[x]][1], sqrt(model$v[[x]][1, 1]),lower.tail = F)),
                                                                 (-pnorm(bin[[2]]$grid[, 2], model$mu[[x]][2], sqrt(model$v[[x]][2, 2]),lower.tail = F) +
                                                                    pnorm(bin[[2]]$grid[, 1], model$mu[[x]][2], sqrt(model$v[[x]][2, 2]),lower.tail = F)))))
  mat4=t(t(mat4)*model$pi)
  return(mat4)
}


#' Computes at the same time and with the same initial parameters binned CL-EM algorithm and classical EM (using Rmixmod)
#'
#' Computes at the same time and with the same initial parameters binned CL-EM algorithm and classical EM (using Rmixmod)
#' @param data Needs to a matrix of data
#' @param class Number of groups of the mixture models
#' @param R Needs to be a matrix or a numerical vector. The element (i,j) indicates the density of the j-th grid along the i-th dimension
#' @param maxit Maximum number of iterations EM
#' @param tol Tolerance of EM
#' @param nrep Number of different initialization of EM
#' @return Binned CL-EM and classical EM algorithm output
#' @export
binmixmod=function(data,class,R,maxit,tol,nrep)
{
  n = NROW(data)
  R = as.matrix(R)
  d = NCOL(data)
  bin = list()
  set.seed(123)
  pb <- txtProgressBar(min = 0, max = length(class) * nrep,
                       style = 3)
  if(d==1){
    for (j in 1:nrow(R)) bin[[j]] =build_bin_mext(data, R[j,])

    ret = list()
    ret1=list()
    tmix=numeric(nrep)
    for (cl in class)  {
      temp = list()
      mix=list()
      for (i in 1:nrep) {
        pi0 = runif(cl)
        pi0 = pi0/sum(pi0)
        mu0 = runif(cl,min(data), max(data))
        sigma0 =runif(cl,0,(var(data)))
        temp[[i]] = list()
        for (j in 1:nrow(R)) {
          temp[[i]][[j]] = try(em_univ(bin[[j]],pi0, mu0, sigma0, maxit,tol), silent = T)
        }
        a.mix=new("GaussianParameter")
        a.mix@proportions=pi0
        a.mix@mean=as.matrix(mu0)
        a.mix@variance=list()
        for(j in 1:cl)
        a.mix@variance[[j]]=t(sigma0[j])
        t0=Sys.time()
        mixmo=mixmodCluster(data,nbCluster=cl,models=mixmodGaussianModel(listModels = "Gaussian_pk_Lk_C"),
                                 strategy=mixmodStrategy(nbTry=1,initMethod="parameter",epsilonInAlgo=tol,nbIterationInAlgo=maxit,parameter=a.mix))
        tmix[i]=as.numeric(difftime(Sys.time(), t0, units = "sec"))
        mixmo@bestResult@proba=matrix(0,1,1)
        mix[[i]]=mixmo@bestResult
        value = (cl - min(class)) * nrep + i
        setTxtProgressBar(pb, value)

      }
      maxi = matrix(-Inf, nrep, nrow(R))
      for (i in 1:nrep) for (j in 1:nrow(R)) {
        if (class(temp[[i]][[j]]) == "list")
          maxi[i, j] = max(temp[[i]][[j]]$logli)
      }
      max_temp = list()
      for (i in 1:nrow(R)) {
        max_temp[[i]] = temp[[which.max(maxi[, i])]][[i]]
        if (class(max_temp[[i]]) == "list")
          max_temp[[i]]$bic = -2 * max(max_temp[[i]]$logli) +
            (ncol(max_temp[[i]]$parm) - 1) * log(n)
      }
      maxi1 = rep(-Inf, nrep)
      for (i in 1:nrep)  {
        maxi1[i] = mix[[i]]@likelihood
      }
      max_temp1 = mix[[which.max(maxi1)]]
      timix_best=tmix[which.max(maxi1)]
      ret[[which(class == cl)]] = list(best = max_temp, all = temp)
      ret1[[which(class==cl)]]=list(best=max_temp1,all=mix,t_best=timix_best,t=tmix)
    }
  } else {
  for (j in 1:ncol(R)) bin[[j]] = sapply(1:d, function(x) build_bin_mext(data[,
                                                                              x], R[x, j]), simplify = F)

  ret = list()
  ret1=list()
  tmix=numeric(nrep)
  for (cl in class)  {
    temp = list()
    mix=list()
    for (i in 1:nrep) {
      pi0 = runif(cl)
      pi0 = pi0/sum(pi0)
      mu0 = list()
      for (j in 1:cl) mu0[[j]] = sapply(1:d, function(x) runif(1,
                                                               min(data[, x]), max(data[, x])))
      sigma0 = list()
      for (j in 1:cl) sigma0[[j]] = diag(sapply(1:d,
                                                function(x) runif(1, 0, (var(data[, x])))))
      temp[[i]] = list()
      for (j in 1:ncol(R)) {
        temp[[i]][[j]] = try(em_bin_gaus_comp(bin[[j]],
                                              pi0, mu0, sigma0, tol, maxit), silent = T)
      }
      a.mix=new("GaussianParameter")
      a.mix@proportions=pi0
      a.mix@mean=matrix(unlist(mu0),ncol=d,byrow=T)
      a.mix@variance=sigma0
      t0=Sys.time()
      mixmo=mixmodCluster(data,nbCluster=cl,models=mixmodGaussianModel(listModels = "Gaussian_pk_Lk_Bk"),
                          strategy=mixmodStrategy(nbTry=1,initMethod="parameter",epsilonInAlgo=tol,nbIterationInAlgo=maxit,parameter=a.mix))
      tmix[i]=as.numeric(difftime(Sys.time(), t0, units = "sec"))
      mixmo@bestResult@proba=matrix(0,1,1)
      mix[[i]]=mixmo@bestResult
      value = (cl - min(class)) * nrep + i
      setTxtProgressBar(pb, value)

    }
    maxi = matrix(-Inf, nrep, ncol(R))
    for (i in 1:nrep) for (j in 1:ncol(R)) {
      if (class(temp[[i]][[j]]) == "list")
        maxi[i, j] = max(temp[[i]][[j]]$logli)
    }
    max_temp = list()
    for (i in 1:ncol(R)) {
      max_temp[[i]] = temp[[which.max(maxi[, i])]][[i]]
      if (class(max_temp[[i]]) == "list")
        max_temp[[i]]$bic = -2 * max(max_temp[[i]]$logli) +
          (ncol(max_temp[[i]]$parm) - 1) * log(n)
    }
    maxi1 = rep(-Inf, nrep)
    for (i in 1:nrep)  {
      maxi1[i] = mix[[i]]@likelihood
    }
    max_temp1 = mix[[which.max(maxi1)]]
    timix_best=tmix[which.max(maxi1)]
    ret[[which(class == cl)]] = list(best = max_temp, all = temp)
    ret1[[which(class==cl)]]=list(best=max_temp1,all=mix,t_best=timix_best,t=tmix)
  }
  }
  close(pb)
  names(ret) = paste("class", class, sep = "")
  names(ret1) = paste("class", class, sep = "")
  return(list(binned=ret,mixmod=ret1))
}
binmixt.ord=function(data, class, R, maxit, tol, nrep,maxit1,print)
{
  b=binmixt(data,class,R,maxit,tol,nrep,print)
  ord=order(b[[1]]$best[[1]]$pi)
  bin.comp=list()
  bin.comp$pi=sort(b[[1]]$best[[1]]$pi)
  bin.comp$mu=list()
  bin.comp$v=list()
  for(k in 1:class)
  {
    bin.comp$mu[[k]]=b[[1]]$best[[1]]$mu[[order(b[[1]]$best[[1]]$pi)[k]]]
    bin.comp$v[[k]]=b[[1]]$best[[1]]$v[[order(b[[1]]$best[[1]]$pi)[k]]]
  }
  return(list(bin=b[[1]][[1]][[1]]$bin,parm=bin.comp,init=b))

}
#' CL-EM with Initialization axes by axes
#'
#' CL-EM with Initialization axes by axes
#' @param data Needs to a matrix of data
#' @param class Number of groups of the mixture models
#' @param R Needs to be a matrix or a numerical vector (vector at moment only option). The element (i,j) indicates the density of the j-th grid along the i-th dimension
#' @param maxit Maximum number of iterations small-EM
#' @param tol Tolerance of EM
#' @param nrep Number of different initialization of small-EM
#' @param maxit1 Maximum number of iterations final EM
#' @return Binned CL-EM algorithm output
#' @export
binmixt.compinit=function(data, class, R, maxit, tol, nrep,maxit1,print=T,only.init=F)
{
  if(missing(print)) print=T
  ###initialization phase
  t0=Sys.time()
  bin.comp=list()
  bin1=list()
  b=sapply(1:ncol(data),function(x) binmixt.ord(data[,x],class,R[x],maxit,1e-20,nrep,maxit1,print),simplify=F )
  pi0=numeric(class)
  mu0=list()
  for(k in 1:class)
    mu0[[k]]=numeric(ncol(data))
  v0=list()
  for(k in 1:class)
    v0[[k]]=diag(ncol(data))
  for(i in 1:ncol(data))
  {
    pi0=pi0+b[[i]]$parm$pi
    for(k in 1:class)
      mu0[[k]][i]=b[[i]]$parm$mu[[k]]
    for(k in 1:class)
      v0[[k]][i,i]=b[[i]]$parm$v[[k]]
  }
  pi0=pi0/ncol(data)
  bin1=list()
  for(i in 1:ncol(data)) bin1[[i]]=b[[i]]$bin
  t1=Sys.time()
  if(only.init==T)
    return(list(pi.init=pi0,mu.init=mu0,v.init=v0,init=b,t.init=as.numeric(difftime(t1,t0,units="sec"))))
  ##estimation phase
  em=em_bin_gaus_comp(bin1,pi0,mu0,v0,tol,maxit1)
  t2=Sys.time()
  return(list(result=em,init=b,t.init=as.numeric(difftime(t1,t0,units="sec")),t.estim=as.numeric(difftime(t2,t1,units="sec"))))
}

#' CL-EM with classic initialization
#'
#' CL-EM with classic initialization
#' @param data Needs to a matrix of data
#' @param class Number of groups of the mixture models
#' @param R Needs to be a matrix or a numerical vector (vector at moment only option). The element (i,j) indicates the density of the j-th grid along the i-th dimension
#' @param maxit Maximum number of iterations small-EM
#' @param tol Tolerance of EM
#' @param nrep Number of different initialization of small-EM
#' @param maxit1 Maximum number of iterations final EM
#' @return Binned CL-EM algorithm output
#' @export
binmixt.classicinit=function(data, class, R, maxit, tol, nrep,maxit1,print=T)
{
  #initialization phase
  if(missing(print)) print=T
  t0=Sys.time()
  bin=binmixt(data,class,R,maxit,1e-20,nrep,print)
  t1=Sys.time()
  #estimation phase
  if(NCOL(data)==1)
    em=em_univ(bin[[1]][[1]][[1]]$bin,bin[[1]][[1]][[1]]$pi1,
               bin[[1]][[1]][[1]]$mu,bin[[1]][[1]][[1]]$v,
               maxit1,tol) else
  em=em_bin_gaus_comp(bin[[1]][[1]][[1]]$bin,bin[[1]][[1]][[1]]$pi,
                               bin[[1]][[1]][[1]]$mu,bin[[1]][[1]][[1]]$v,
                               tol,maxit1)
  t2=Sys.time()
  return(list(result=em,init=bin,t.init=as.numeric(difftime(t1,t0,units="sec")),t.estim=as.numeric(difftime(t2,t1,units="sec"))))
}


#'Rmixmod parameters to BinMixt format
#'@export
Rmix_to_binM=function(model)
{
  pi0=model@parameters@proportions
  mu0=list()
  for(i in 1:length(pi0))
    mu0[[i]]=model@parameters@mean[i,]
  s0=list()
  for(i in 1:length(pi0))
    s0[[i]]=model@parameters@variance[[i]]
  return(list(pi=pi0,mu=mu0,v=s0))
}

#'Inner function of hybrid binned CL-EM algorithm for raw proportion
#'and means calculation and other useful quantities
em_raw_inner=function(subs,pi0,mu0,sigma0)
{
  class=length(pi0)
  s=t(sapply(1:class,function(x) dmvnorm(subs,mean=mu0[[x]],sigma=((sigma0[[x]])))))
  s1=as.numeric(pi0%*%s)
  tau=t(t(pi0*s)/s1)
  sum1=apply(tau,1,sum)
  int2=apply(tau,1,function(x) x%*%subs)
  return(list(pi.raw=sum1,mu.raw=int2,tau.raw=tau))
}


#'Inner function of hybrid binned CL-EM algorithm for proportion
#'and means calculation and other useful quantities
em_dim_binned_pimu=function(bin,subs,pi0,mu0,sigma0,d,raw.results)
{
  class=length(pi0)
  bin1=bin[[d]]
  ext1=bin1$grid
  sub=subs[,d]
  mat=cbind(as.numeric(lapply(mu0,function(x) x[d])),as.numeric(lapply(sigma0,function(x) x[d,d])))

  ##binned part
  s=apply(ext1, 1, function(x) apply(mat, 1, function(y) {
    max(pnorm(x[2],y[1], sqrt(y[2])) - pnorm(x[1], y[1], sqrt(y[2])),-pnorm(x[2],y[1], sqrt(y[2]),lower.tail=F)+ pnorm(x[1], y[1], sqrt(y[2]),lower.tail=F))
  }))
  s1=as.numeric(pi0%*%s)
  sm=apply(ext1,1,function(x) apply(mat,1,function(y) (dnorm(x[2],y[1],sqrt(y[2]))-dnorm(x[1],y[1],sqrt(y[2])))))
  sv=apply(ext1,1,function(x) apply(mat,1,function(y)
  {
    if(x[1]==-Inf)
      return(x[2]*dnorm(x[2],mean=y[1],sqrt(y[2])))
    if(x[2]==Inf)
      return(-x[1]*dnorm(x[1],mean=y[1],sqrt(y[2])))
    if(x[1]!=-Inf & x[2]!=Inf)
      return(x[2]*dnorm(x[2],mean=y[1],sqrt(y[2]))-x[1]*dnorm(x[1],mean=y[1],sqrt(y[2])))
  }
  ))
  int=t(t(s*pi0)/s1)
  pp=NULL
  for(i in 1:class)
    pp=rbind(pp,s1)
  int3=(1/(pp/pi0))*(s*mat[,1]-sm*mat[,2])
  sum1=apply(int,1,function(x) sum(bin1$bin*x))
  sum2=apply(int3,1,function(x) sum(bin1$bin*x))
  mu1=(sum2+raw.results$mu.raw[d,])/(sum1+raw.results$pi.raw)
  int4=((1/(pp/pi0))*((s+sm*(2*mu1-mat[,1])-sv)*mat[,2]+s*(mu1-mat[,1])^2))
  sum3=apply(int4,1,function(x) sum(bin1$bin*x))

  v1=sum3/sum1
  return(list(sum.binned.prop=sum1,mu1=mu1,sum.binned.var=sum3,mat=mat,sub=sub))
}

#'Inner function of hybrid binned CL-EM algorithm
em_dim_var=function(d,subs,raw.results,binned.results)
{
  class=length(binned.results$sum.binned.prop)
  sum.var.raw=numeric(class)
  for(i in 1:class)
  sum.var.raw[i]=raw.results$tau.raw[i,]%*%(binned.results$sub-binned.results$mu1[i])^2
  v1=(sum.var.raw+binned.results$sum.binned.var)/(binned.results$sum.binned.prop+raw.results$pi.raw)
  return(list(v1=v1,mu1=binned.results$mu1,pi.raw=raw.results$pi.raw,pi.bin=binned.results$sum.binned.prop))
}


#' Hybrid Binned CL-EM algorithm for multivariate binned diagonal Gaussian mixture
#'
#' Hybrid Binned CL-EM algorithm for multivariate binned diagonal Gaussian mixture
#' @param bin  Binned data
#' @param subs Sub-sample
#' @param pi0 Starting point for proportions
#' @param mu0 Starting point for means
#' @param sigma0 Starting point for variances
#' @param tol Tolerance of the algorithm
#' @param tim Maximum number of iterations
#' @return List with original binned data, final estimates, log-likelihood and parameters sequences and medium time per iteration
#' @export
em_bin_gaus_hybrid=function(bin,subs,pi0,mu0,sigma0,tol,tim,print=F)
{
  t0=Sys.time()
  flag=1
  j=1
  R=as.numeric(lapply(bin,function(x) length(x$bin)-1))
  n=sum(bin[[1]]$bin)
  bin=lapply(bin,function(x) nozero(x))
  class=length(pi0)
  d=length(mu0[[1]])
  kl=numeric(tim)
  parm=matrix(0,tim+1,2*d*class+class)
  logli=numeric(tim+1)
  a=logli_marg(bin,pi0,mu0,sigma0)
  b=logli_mixt(subs,pi0,mu0,sigma0)
#  if(is.na(a) | is.nan(a)) return(list(mess="Errore vero binned",pi=pi0,mu=mu0,v=sigma0))
 # if(is.na(b) | is.nan(b)) return(list(mess="Errore vero raw",pi=pi0,mu=mu0,v=sigma0))
  logli[j]=a+b
  parm[j,]=c(pi0,unlist(mu0),unlist(sigma0)[unlist(sigma0)!=0])
  while(flag==1)
  {
    if(print==T) cat(j," ")
    raw.results=em_raw_inner(subs,pi0,mu0,sigma0)
    dim1=sapply(1:d,function(x)
      {
      binned.results=em_dim_binned_pimu(bin,subs,pi0,mu0,sigma0,x,raw.results)
      em_dim_var(x,subs,raw.results,binned.results)
      },simplify = F)
    pi1=(dim1[[1]]$pi.raw+colSums(matrix(unlist(lapply(dim1,function(x) x$pi.bin)),byrow=T,ncol=class)))/(d*n+NROW(subs))
    mu1=list()
    sigma1=list()
    for(i in 1:class)
    {
      mu1[[i]]=as.numeric(lapply(dim1,function(x) x$mu1[[i]]))
      sigma1[[i]]=diag(as.numeric(lapply(dim1,function(x) x$v1[[i]])))
    }
    j=j+1
    a=logli_marg(bin,pi1,mu1,sigma1)
      b=logli_mixt(subs,pi1,mu1,sigma1)
      #if(is.na(a) | is.nan(a)) return(list(mess="Errore vero binned",pi=pi0,mu=mu0,v=sigma0))
      #if(is.na(b) | is.nan(b)) return(list(mess="Errore vero raw",pi=pi0,mu=mu0,v=sigma0))
    logli[j]=a+b
    mu1_ul=unlist(mu1)
    sigma1_ul=unlist(sigma1)
    sigma1_ul=sigma1_ul[sigma1_ul!=0]

    parm[j,]=c(pi1,mu1_ul,sigma1_ul)
    flag=0
    #print(paste(paste("Iter ",j),paste(c(pi1,mu1, sigma1_ul),collapse=" "),sep=": "))
    if(abs((logli[j]-logli[j-1])/logli[j-1])>tol & j<=tim)
    {
      flag=1
      pi0=pi1
      mu0=mu1
      sigma0=sigma1
    }
  }
  conv=TRUE
  if(j>tim) conv=FALSE
  logli=logli[1:min(j,tim+1)]
  return(list(bin=bin,subs=subs,density=R,pi=pi1,mu=mu1,v=sigma1,logli=logli,parm=parm[1:min(j,tim+1),],convergence=conv,timem=as.numeric(difftime(Sys.time(),t0,units="sec"))/(length(logli)-1)))
}

#' Binned hybrid CL-EM algorithm
#'
#' Binned hybrid CL-EM algorithm
#' @param data Needs to a matrix of data
#' @param subind Sample indices
#' @param R Needs to be a matrix or a numerical vector. The element (i,j) indicates the density of the j-th grid along the i-th dimension
#' @param class Number of groups of the mixture models
#' @param maxit Maximum number of iterations EM
#' @param tol Tolerance of EM
#' @param nrep Number of different initialization of EM
#' @return Binned CL-EM algorithm output
#' @export
binmixt.hybrid=function (data,subind, class, R, maxit, tol, nrep,init.type,print=T)
{
  if(missing(print)) print=T
  n = NROW(data)
  R = as.matrix(R)
  d=NCOL(data)
  if(NCOL(data)==1) subs=data[subind] else
  subs=data[subind,]
  set.seed(123)
  if(print==T)
  {
    pb <- txtProgressBar(min = 0, max = length(class)*nrep,
                         style = 3)
    if(d==1)
    {
      bin=list()
      for(j in 1:nrow(R))
        bin[[j]] = build_bin_mext(data, R[j,])
      ret = list()
      for (cl in class) {
        if (cl == 1) {
          temp = list()
          for (i in 1:nrep) {
            mu0 = sapply(1:d, function(x) runif(1, min(data[,
                                                            x]), max(data[, x])))
            sigma0 = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                      x])))))
            temp[[i]] = list()
            for (j in 1:nrow(R)) {
              temp[[i]][[j]] = try(em_univ_hybrid(bin[[j]],subs,m0=mu0, s0=sigma0,tim= maxit,tol=tol), silent = T)
            }
            value=(cl-min(class))*nrep+i
            setTxtProgressBar(pb,value)
          }
        }
        else {
          temp = list()
          for (i in 1:nrep) {
            if(init.type=="casual")
            {
              pi0 = runif(cl)
              pi0 = pi0/sum(pi0)
              mu0=runif(cl,min=min(data),max=max(data))
              sigma0=runif(cl,min=0,max=var(data))
            }
            if(init.type=="VarData")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0 = subs[samp]
              sigma0 =rep(var(data),cl)*runif(cl)
            }
            if(init.type=="VarData.div")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0 = subs[samp]
              sigma0 =rep(var(data)/(cl^2),cl)*runif(cl)
            }
            if(init.type=="VarData.eq")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0 = subs[samp]
              sigma0 =rep(var(data),cl)
            }
            if(init.type=="VarData.eq.div")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0 = subs[samp]
              sigma0 =rep(var(data)/(cl^2),cl)
            }
            if(init.type=="VarGrid")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0 = subs[samp]
              sigma0 =rep(var.grid.pond(bin[[1]]),cl)*runif(cl)
            }
            if(init.type=="VarGrid.div")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0 = subs[samp]
              sigma0 =rep(var.grid.pond(bin[[1]])/(cl^2),cl)*runif(cl)
            }
            if(init.type=="VarGrid.eq")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0 = subs[samp]
              sigma0 =rep(var.grid.pond(bin[[1]]),cl)
            }
            if(init.type=="VarGrid.eq.div")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0 = subs[samp]
              sigma0 =rep(var.grid.pond(bin[[1]])/(cl^2),cl)
            }
            if(init.type=="Rmixmod")
            {
              mixmo=mixmodCluster(as.data.frame(x$x[samp,]),nbCluster=30,models=mixmodGaussianModel(listModels = "Gaussian_pk_Lk_Bk"),strategy=mixmodStrategy(nbTryInInit=1,nbIterationInInit=1,nbIterationInAlgo=1))
             while(length(mixmo@bestResult@parameters@proportions)==0)
               mixmo=mixmodCluster(as.data.frame(x$x[samp,]),nbCluster=30,models=mixmodGaussianModel(listModels = "Gaussian_pk_Lk_Bk"),strategy=mixmodStrategy(nbTryInInit=1,nbIterationInInit=1,nbIterationInAlgo=1))
             param=Rmix_to_binM(mixmo@bestResult)
             pi0=param$pi
             mu0=param$mu
             sigma0=param$v
              }

            temp[[i]] = list()
            for (j in 1:nrow(R)) {
              temp[[i]][[j]] = try(em_univ_hybrid(bin[[j]],subs,pi0, mu0, sigma0, maxit,tol), silent = T)
            }
            value=(cl-min(class))*nrep+i
            setTxtProgressBar(pb,value)
          }

        }
        maxi = matrix(-Inf, nrep, nrow(R))
        for (i in 1:nrep) for (j in 1:nrow(R)) {
          if (class(temp[[i]][[j]]) == "list")
            maxi[i, j] = max(temp[[i]][[j]]$logli)
        }
        max_temp = list()
        for (i in 1:nrow(R)) {
          max_temp[[i]] = temp[[which.max(maxi[, i])]][[i]]
          if (class(max_temp[[i]]) == "list")
            max_temp[[i]]$bic = -2 * max(max_temp[[i]]$logli) +
              (ncol(max_temp[[i]]$parm) - 1) * log(n)
        }
        ret[[which(class == cl)]] = list(best=max_temp,all=temp)
      }
    } else
      {
      bin=list()
      for(j in 1:ncol(R))
        bin[[j]] = sapply(1:d, function(x) build_bin_mext(data[,
                                                               x], R[x, j]), simplify = F)


      ret = list()
      for (cl in class) {
        if (cl == 1) {
          temp = list()
          for (i in 1:nrep) {
            mu0 = sapply(1:d, function(x) runif(1, min(data[,
                                                            x]), max(data[, x])))
            sigma0 = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                      x])))))
            temp[[i]] = list()
            for (j in 1:ncol(R)) {
              temp[[i]][[j]] = try(em_bin_gaus_comp_nomixt(bin[[j]],
                                                           mu0, sigma0, tol, maxit), silent = T)
            }
            value=(cl-min(class))*nrep+i
            setTxtProgressBar(pb,value)
          }
        }
        else {
          temp = list()
          for (i in 1:nrep) {
            if(init.type=="casual")
            {
              pi0 = runif(cl)
              pi0 = pi0/sum(pi0)
              mu0 = list()
              for (j in 1:cl) mu0[[j]] = sapply(1:d, function(x) runif(1,
                min(data[, x]), max(data[, x])))
              sigma0 = list()
              for (j in 1:cl) sigma0[[j]] = diag(sapply(1:d, function(x) runif(1,0,(var(data[,x])))))

            }
            if(init.type=="VarData")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0=list()
              for(j in 1:cl)
                mu0[[j]]=subs[samp[j],]
              sigma0 = list()
              for(j in 1:cl)
              {
                sigma0[[j]]=diag(diag(var(data)))*runif(d)
               }
            }
            if(init.type=="VarData.div")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0=list()
              for(j in 1:cl)
                mu0[[j]]=subs[samp[j],]
              sigma0 = list()
              for(j in 1:cl)
              {
                sigma0[[j]]=diag(diag(var(data)))*runif(d)/(cl^2)
              }
            }
            if(init.type=="VarData.eq")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0=list()
              for(j in 1:cl)
                mu0[[j]]=subs[samp[j],]
              sigma0 = list()
              for(j in 1:cl)
              {
                sigma0[[j]]=diag(diag(var(data)))
              }
            }
            if(init.type=="VarData.eq.div")
            {
              pi0=rep(1/cl,cl)
              samp=sample(1:length(subind),cl,replace=F)
              mu0=list()
              for(j in 1:cl)
                mu0[[j]]=subs[samp[j],]
              sigma0 = list()
              for(j in 1:cl)
              {
                sigma0[[j]]=diag(diag(var(data)))/(cl^2)
              }
            }
            if(init.type=="VarGrid")
            {
              pi0=rep(1/cl,cl)

              samp=sample(1:length(subind),cl,replace=F)
              mu0=list()
              for(j in 1:cl)
                mu0[[j]]=subs[samp[j],]
              sigma0 = list()
              for(j in 1:cl)
              {
                sigma0[[j]]=diag(d)
                for(k in 1:d) sigma0[[j]][k,k]=var(bin[[1]][[k]]$grid[2:(nrow(bin[[1]][[k]]$grid)),1])*runif(1)
              }
              }
            if(init.type=="VarGrid.div")
            {
              pi0=rep(1/cl,cl)

              samp=sample(1:length(subind),cl,replace=F)
              mu0=list()
              for(j in 1:cl)
                mu0[[j]]=subs[samp[j],]
              sigma0 = list()
              for(j in 1:cl)
              {
                sigma0[[j]]=diag(d)
                for(k in 1:d) sigma0[[j]][k,k]=var.grid.pond(bin[[1]][[k]])*runif(1)/(cl^2)
              }
              }
            if(init.type=="VarGrid.eq")
            {
              pi0=rep(1/cl,cl)

              samp=sample(1:length(subind),cl,replace=F)
              mu0=list()
              for(j in 1:cl)
                mu0[[j]]=subs[samp[j],]
              sigma0 = list()
              for(j in 1:cl)
              {
                sigma0[[j]]=diag(d)
                for(k in 1:d) sigma0[[j]][k,k]=var.grid.pond(bin[[1]][[k]])
              }
            }
            if(init.type=="VarGrid.eq.div")
            {
              pi0=rep(1/cl,cl)

              samp=sample(1:length(subind),cl,replace=F)
              mu0=list()
              for(j in 1:cl)
                mu0[[j]]=subs[samp[j],]
              sigma0 = list()
              for(j in 1:cl)
              {
                sigma0[[j]]=diag(d)
                for(k in 1:d) sigma0[[j]][k,k]=var.grid.pond(bin[[1]][[k]])/(cl^2)
              }
            }
            if(init.type=="Rmixmod")
            {
            mixmo=mixmodCluster(as.data.frame(subs),nbCluster=cl,models=mixmodGaussianModel(listModels = "Gaussian_pk_Lk_Bk"),strategy=mixmodStrategy(initMethod="random",nbTryInInit=1,nbIterationInAlgo=1))
            while(length(mixmo@bestResult@parameters@proportions)==0)
              mixmo=mixmodCluster(as.data.frame(subs),nbCluster=cl,models=mixmodGaussianModel(listModels = "Gaussian_pk_Lk_Bk"),strategy=mixmodStrategy(nbTryInInit=1,nbIterationInInit=1,nbIterationInAlgo=1))
            param=Rmix_to_binM(mixmo@bestResult)
            pi0=param$pi
            mu0=param$mu
            sigma0=param$v
            }
            temp[[i]] = list()
            for (j in 1:ncol(R)) {
              temp[[i]][[j]] = try(em_bin_gaus_hybrid(bin[[j]],subs,
                                                    pi0, mu0, sigma0, tol, maxit), silent = T)
            }
            value=(cl-min(class))*nrep+i
            setTxtProgressBar(pb,value)
          }

        }
        maxi = matrix(-Inf, nrep, ncol(R))
        for (i in 1:nrep) for (j in 1:ncol(R)) {
          if (class(temp[[i]][[j]]) == "list")
            maxi[i, j] = max(temp[[i]][[j]]$logli)
        }
        max_temp = list()
        for (i in 1:ncol(R)) {
          max_temp[[i]] = temp[[which.max(maxi[, i])]][[i]]
          if (class(max_temp[[i]]) == "list")
            max_temp[[i]]$bic = -2 * max(max_temp[[i]]$logli) +
              (ncol(max_temp[[i]]$parm) - 1) * log(n)
        }
        ret[[which(class == cl)]] = list(best=max_temp,all=temp)
      } }
    close(pb)
  } else
    {
    if(d==1)
    {
      bin=list()
      for(j in 1:nrow(R))
        bin[[j]] = build_bin_mext(data, R[j,])
      ret = list()
      for (cl in class) {
        if (cl == 1) {
          temp = list()
          for (i in 1:nrep) {
            mu0 = sapply(1:d, function(x) runif(1, min(data[,
                                                            x]), max(data[, x])))
            sigma0 = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                      x])))))
            temp[[i]] = list()
            for (j in 1:nrow(R)) {
              temp[[i]][[j]] = try(em_univ_subs(bin[[j]],subs,m0=mu0, s0=sigma0,tim= maxit,tol=tol), silent = T)
            }

          }
        }
        else {
          temp = list()
          for (i in 1:nrep) {
            pi0 = runif(cl)
            pi0 = pi0/sum(pi0)
            mu0 = runif(cl,min(data), max(data))
            sigma0 =runif(cl,0,(var(data)))/(cl^2)
            temp[[i]] = list()
            for (j in 1:nrow(R)) {
              temp[[i]][[j]] = try(em_univ_hybrid(bin[[j]],subs,pi0, mu0, sigma0, maxit,tol), silent = T)
            }
            value=(cl-min(class))*nrep+i
          }

        }
        maxi = matrix(-Inf, nrep, nrow(R))
        for (i in 1:nrep) for (j in 1:nrow(R)) {
          if (class(temp[[i]][[j]]) == "list")
            maxi[i, j] = max(temp[[i]][[j]]$logli)
        }
        max_temp = list()
        for (i in 1:nrow(R)) {
          max_temp[[i]] = temp[[which.max(maxi[, i])]][[i]]
          if (class(max_temp[[i]]) == "list")
            max_temp[[i]]$bic = -2 * max(max_temp[[i]]$logli) +
              (ncol(max_temp[[i]]$parm) - 1) * log(n)
        }
        ret[[which(class == cl)]] = list(best=max_temp,all=temp)
      }
    } else {
      bin=list()
      for(j in 1:ncol(R))
        bin[[j]] = sapply(1:d, function(x) build_bin_mext(data[,
                                                               x], R[x, j]), simplify = F)


      ret = list()
      for (cl in class) {
        if (cl == 1) {
          temp = list()
          for (i in 1:nrep) {
            mu0 = sapply(1:d, function(x) runif(1, min(data[,
                                                            x]), max(data[, x])))
            sigma0 = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                      x])))))
            temp[[i]] = list()
            for (j in 1:ncol(R)) {
              temp[[i]][[j]] = try(em_bin_gaus_comp_nomixt(bin[[j]],
                                                           mu0, sigma0, tol, maxit), silent = T)
            }
          }
        }
        else {
          temp = list()
          for (i in 1:nrep) {
            pi0 = runif(cl)
            pi0 = pi0/sum(pi0)
            mu0 = list()
            for (j in 1:cl) mu0[[j]] = sapply(1:d, function(x) runif(1,
                                                                     min(data[, x]), max(data[, x])))
            sigma0 = list()
            for (j in 1:cl) sigma0[[j]] = diag(sapply(1:d, function(x) runif(1,0,(var(data[,
                                                                                           x])))))
            temp[[i]] = list()
            for (j in 1:ncol(R)) {
              temp[[i]][[j]] = try(em_bin_gaus_hybrid(bin[[j]],subs,
                                                    pi0, mu0, sigma0, tol, maxit), silent = T)
            }
          }

        }
        maxi = matrix(-Inf, nrep, ncol(R))
        for (i in 1:nrep) for (j in 1:ncol(R)) {
          if (class(temp[[i]][[j]]) == "list")
            maxi[i, j] = max(temp[[i]][[j]]$logli)
        }
        max_temp = list()
        for (i in 1:ncol(R)) {
          max_temp[[i]] = temp[[which.max(maxi[, i])]][[i]]
          if (class(max_temp[[i]]) == "list")
            max_temp[[i]]$bic = -2 * max(max_temp[[i]]$logli) +
              (ncol(max_temp[[i]]$parm) - 1) * log(n)
        }
        ret[[which(class == cl)]] = list(best=max_temp,all=temp)
      } } }
  names(ret) = paste("class", class, sep = "")
  return(ret)
}

binmixt.hybrid.ord=function(data, subind,class, R, maxit, tol, nrep,maxit1,init.type,print)
{
  b=binmixt.hybrid(data,subind,class,R,maxit,tol,nrep,init.type,print)
  ord=order(b[[1]]$best[[1]]$pi)
  bin.comp=list()
  bin.comp$pi=sort(b[[1]]$best[[1]]$pi)
  bin.comp$mu=list()
  bin.comp$v=list()
  for(k in 1:class)
  {
    bin.comp$mu[[k]]=b[[1]]$best[[1]]$mu[[order(b[[1]]$best[[1]]$pi)[k]]]
    bin.comp$v[[k]]=b[[1]]$best[[1]]$v[[order(b[[1]]$best[[1]]$pi)[k]]]
  }
  return(list(bin=b[[1]][[1]][[1]]$bin,subs=b[[1]][[1]][[1]]$subs,parm=bin.comp,init=b))

}
#' CL-EM with Initialization axes by axes
#'
#' CL-EM with Initialization axes by axes
#' @param data Needs to a matrix of data
#' @param class Number of groups of the mixture models
#' @param R Needs to be a matrix or a numerical vector (vector at moment only option). The element (i,j) indicates the density of the j-th grid along the i-th dimension
#' @param maxit Maximum number of iterations small-EM
#' @param tol Tolerance of EM
#' @param nrep Number of different initialization of small-EM
#' @param maxit1 Maximum number of iterations final EM
#' @return Binned CL-EM algorithm output
#' @export
binmixt.hybrid.compinit=function(data, subind, class, R, maxit, tol, nrep,maxit1,init.type,print=T,only.init=F)
{
  if(missing(print)) print=T
  ###initialization phase
  t0=Sys.time()
  bin.comp=list()
  bin1=list()
  b=sapply(1:ncol(data),function(x) binmixt.hybrid.ord(data[,x],subind,class,R[x],maxit,1e-20,nrep,maxit1,init.type,print),simplify=F )
  pi0=numeric(class)
  mu0=list()
  for(k in 1:class)
    mu0[[k]]=numeric(ncol(data))
  v0=list()
  for(k in 1:class)
    v0[[k]]=diag(ncol(data))
  for(i in 1:ncol(data))
  {
    pi0=pi0+b[[i]]$parm$pi
    for(k in 1:class)
      mu0[[k]][i]=b[[i]]$parm$mu[[k]]
    for(k in 1:class)
      v0[[k]][i,i]=b[[i]]$parm$v[[k]]
  }
  pi0=pi0/ncol(data)
  bin1=list()
  for(i in 1:ncol(data)) bin1[[i]]=b[[i]]$bin
  subs=matrix(0,nrow=length(subind),ncol=ncol(data))
  for(i in 1:ncol(data)) subs[,i]=b[[i]]$subs
  t1=Sys.time()
  if(only.init==T)
    return(list(pi.init=pi0,mu.init=mu0,v.init=v0,init=b,t.init=as.numeric(difftime(t1,t0,units="sec"))))
  ##estimation phase
  em=em_bin_gaus_hybrid(bin1,subs,pi0,mu0,v0,tol,maxit1)
  t2=Sys.time()
  return(list(result=em,init=b,t.init=as.numeric(difftime(t1,t0,units="sec")),t.estim=as.numeric(difftime(t2,t1,units="sec"))))
}

#' CL-EM with classic initialization
#'
#' CL-EM with classic initialization
#' @param data Needs to a matrix of data
#' @param subsize subsize
#' @param class Number of groups of the mixture models
#' @param R Needs to be a matrix or a numerical vector (vector at moment only option). The element (i,j) indicates the density of the j-th grid along the i-th dimension
#' @param maxit Maximum number of iterations small-EM
#' @param tol Tolerance of EM
#' @param nrep Number of different initialization of small-EM
#' @param maxit1 Maximum number of iterations final EM
#' @return Binned CL-EM algorithm output
#' @export
binmixt.hybrid.classicinit=function(data, subind, class, R, maxit, tol, nrep,maxit1,init.type,print=T,only.init=F)
{
  #initialization phase
  if(missing(print)) print=T
  t0=Sys.time()
  bin=binmixt.hybrid(data,subind,class,R,maxit,1e-20,nrep,init.type,print)
  t1=Sys.time()
  #estimation phase
  if(only.init==T)   return(list(init=bin,t.init=as.numeric(difftime(t1,t0,units="sec"))))

  if(NCOL(data)==1)
    em=em_univ(bin[[1]][[1]][[1]]$bin,bin[[1]][[1]][[1]]$pi1,
               bin[[1]][[1]][[1]]$mu,bin[[1]][[1]][[1]]$v,
               maxit1,tol) else
                 em=em_bin_gaus_hybrid(bin[[1]][[1]][[1]]$bin,bin[[1]][[1]][[1]]$subs,bin[[1]][[1]][[1]]$pi,
                                     bin[[1]][[1]][[1]]$mu,bin[[1]][[1]][[1]]$v,
                                     tol,maxit1)
  t2=Sys.time()
  return(list(result=em,init=bin,t.init=as.numeric(difftime(t1,t0,units="sec")),t.estim=as.numeric(difftime(t2,t1,units="sec"))))
}


#' Hybrid Binned-EM algorithm for univariate Gaussian mixtures
#'
#' Hybrid Binned-EM algorithm for univariate Gaussian mixtures
#' @param bin Univariate binned data
#' @param subs Subsample
#' @param pi0 Starting point for proportions. If missing, a binned_EM algorithm for a single Gaussian is performed
#' @param mu0 Starting point for means (numerical vector)
#' @param sigma0 Starting point for variances (numerical vector)
#' @param tim Maximum number of iterations
#' @param tol Tolerance of the algorithm
#' @return List with final estimates, log-likelihood sequence, parameters sequence, BIC and medium time per iteration
#' @export
em_univ_hybrid=function (bin, subs,pi0, m0, s0,tim,tol)
{
  if(missing(pi0))
  {
    t0=Sys.time()
    bin1 = nozero(bin)
    n=sum(bin1$bin)
    ext1 = bin1$grid
    logli=numeric(tim)
    parm=matrix(0,tim+1,2)
    parm[1,]=c(m0,s0)
    flag=1
    j=1
    logli[j]=logli_mult(bin1,1,m0,s0)
    while(flag==1)
    {
      rip1=apply(ext1,1,function(x) (pnorm(x[2],mean=m0,sqrt(s0))-pnorm(x[1],mean=m0,sqrt(s0))))
      nor1=apply(ext1,1,function(x) (dnorm(x[2],mean=m0,sqrt(s0))-dnorm(x[1],mean=m0,sqrt(s0))))
      c1=apply(ext1,1,function(x)
      {
        if(x[1]==-Inf)
          return(x[2]*dnorm(x[2],mean=m0,sqrt(s0)))
        if(x[2]==Inf)
          return(-x[1]*dnorm(x[1],mean=m0,sqrt(s0)))
        if(x[1]!=-Inf & x[2]!=Inf)
          return(x[2]*dnorm(x[2],mean=m0,sqrt(s0))-x[1]*dnorm(x[1],mean=m0,sqrt(s0)))
      }
      )
      mu1=sum(bin1$bin*(m0-s0*nor1/rip1))/n
      sigma1=sum(bin1$bin*((s0*(rip1+(2*mu1-m0)*nor1-c1)+rip1*(mu1-m0)^2)/rip1))/n
      j=j+1
      parm[j,]=c(mu1,sigma1)
      logli[j]=logli_mult(bin1,1,mu1,sigma1)
      if(abs((logli[j]-logli[j-1])/logli[j-1])>tol &
         j <= tim)
      {
        s0=v1
        m0=mu1
        flag=1
      }
    }
    conv = TRUE
    if (j > tim)
      conv = FALSE
    logli = logli[1:min(j, tim + 1)]
    parm=parm[1:min(j, tim + 1),]
    t1=Sys.time()
    return(list(mu = mu1, v = sigma1,logli=lo,parm=parm,bic=bic,timem=as.numeric(difftime(t1,t0,units="sec"))/(length(logli)-1)))} else {
      t0=Sys.time()
      class = length(pi0)
      bin1 = nozero(bin)
      ext1 = bin1$grid
      n=sum(bin1$bin)
      logli=numeric(tim)
      parm=matrix(0,tim+1,3*class)
      parm[1,]=c(pi0,m0,s0)
      flag=1
      j=1
      logli[j]=logli_mult(bin1,pi0,m0,s0)+logli_mixt(subs,pi0,m0,s0)

      while(flag==1)
      {
        mat=cbind(m0,s0)
        s = apply(ext1, 1, function(x) apply(mat, 1, function(y) {
          max(pnorm(x[2],y[1], sqrt(y[2])) - pnorm(x[1], y[1], sqrt(y[2])),-pnorm(x[2],y[1], sqrt(y[2]),lower.tail=F)+ pnorm(x[1], y[1], sqrt(y[2]),lower.tail=F))
        }))
        s1 = as.numeric(pi0 %*% s)
        sm = apply(ext1, 1, function(x) apply(mat, 1, function(y) (dnorm(x[2], y[1], sqrt(y[2])) - dnorm(x[1], y[1], sqrt(y[2])))))
        sv = apply(ext1, 1, function(x) apply(mat, 1, function(y) {
          if (x[1] == -Inf)
            return(x[2] * dnorm(x[2], mean = y[1], sqrt(y[2])))
          if (x[2] == Inf)
            return(-x[1] * dnorm(x[1], mean = y[1], sqrt(y[2])))
          if (x[1] != -Inf & x[2] != Inf)
            return(x[2] * dnorm(x[2], mean = y[1], sqrt(y[2])) -
                     x[1] * dnorm(x[1], mean = y[1], sqrt(y[2])))
        }))
        int = t(t(s * pi0)/s1)
        pp = NULL
        for (i in 1:class) pp = rbind(pp, s1)
        int3 = (1/(pp/pi0)) * (s * mat[, 1] - sm * mat[, 2])
        sum1 = apply(int, 1, function(x) sum(bin1$bin * x))
        s.raw=t(sapply(1:class,function(x) dnorm(subs,m0[x],sqrt((s0[x])))))
        s1.raw=as.numeric(pi0%*%s.raw)
        tau=t(t(pi0*s.raw)/s1.raw)
        sum1.raw=apply(tau,1,sum)
        int2.raw=apply(tau,1,function(x) x%*%subs)
        pi1=(sum1+sum1.raw)/(sum(bin1$bin)+NROW(subs))
        sum2 = apply(int3, 1, function(x) sum(bin1$bin * x))
        mu1 = (sum2+int2.raw)/(sum1+sum1.raw)
        int4 = ((1/(pp/pi0)) * ((s + sm * (2 * mu1 - mat[, 1]) -
                                   sv) * mat[, 2] + s * (mu1 - mat[, 1])^2))
        sum3 = apply(int4, 1, function(x) sum(bin1$bin * x))
        sum.var.raw=numeric(class)
        for(i in 1:class)
          sum.var.raw[i]=tau[i,]%*%(subs-mu1[i])^2
        v1=(sum.var.raw+sum3)/(sum1+sum1.raw)
        j=j+1
        parm[j,]=c(pi1,mu1,v1)
        logli[j]=logli_mult(bin1,pi1,mu1,v1)+logli_mixt(subs,pi1,mu1,v1)
        flag=0
        if(abs((logli[j]-logli[j-1])/logli[j-1])>tol &
           j <= tim)
        {
          s0=v1
          m0=mu1
          pi0=pi1
          flag=1
        }
      }
      conv = TRUE
      if (j > tim)
        conv = FALSE
      logli = logli[1:min(j, tim + 1)]
      parm=parm[1:min(j, tim + 1),]
      t1=Sys.time()
      return(list(bin=bin1,subs=subs,mu = mu1, v = v1, pi1 = pi1,logli=logli,parm=parm,timem=as.numeric(difftime(t1,t0,units="sec"))/(length(logli)-1)))}
}

#' Binmixt paramter to Rmixmod parameter
#' @param binmixt parameter
#' Binmixt paramter to Rmixmod parameter
#' @export
binmixt_to_Rmix=function(binmixt)
{
  d=length(binmixt$mu[[1]])
  a.mix = new("GaussianParameter")
  a.mix@proportions = binmixt$pi
  a.mix@mean = matrix(unlist(binmixt$mu), ncol = d, byrow = T)
  a.mix@variance = binmixt$v
  a.mix
}
