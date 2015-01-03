takeStep = function(i1, j1, k,y,  b,alphas,eps,C,error) {
  if (i1 == j1) {
    return(list(val=0))
  }
  alphai = alphas[i1]
  alphaj = alphas[j1]
  yi = y[i1]
  yj = y[j1]
  ui = sum(alphas*y*k[i1,]) - b
  uj = sum(alphas*y*k[j1,]) - b
  #   ui = x[i1,]%*%w - b
  #   uj = x[j1,]%*%w - b
  Ei = ui - yi
  Ej = uj - yj
  
  #compute H and L use old alpha i and j
  if(y[i1] == y[j1]) {
    L = max(0, alphas[j1] + alphas[i1] - C)
    H = min(C,alphas[i1] + alphas[j1])
  } else {
    L = max(0, alphas[j1] - alphas[i1])
    H = min(C, C + alphas[j1] - alphas[i1])
  }
  if (L == H) {
    return(list(val=0))
  }
  eta = 2*k[i1,j1] - k[i1,i1] - k[j1,j1]
  
  
  #compute new alpha j and i as alphaj and alphai
  if (eta < 0) {
    alphaj = alphaj - yj*(Ei - Ej)/eta
    if (alphaj < L) {
      alphaj = L
    } else if (alphaj > H) {
      alphaj = H
    }
  } else {
    Lobj = objfun(alphas,j1,L,y,k)
    Hobj = objfun(alphas,j1,H,y,k)
    if (Lobj > Hobj + eps){
      alphaj = L
    } else if (Lobj < Hobj - eps) {
      alphaj = H
    } 
  }
  if (alphaj < 1e-8) {
    alphaj = 0
  } else if (alphaj > C-1e-8) {
    alphaj = C
  }
  if (abs(alphaj - alphas[j1]) < eps*(alphaj+alphas[j1]+eps)) {
    return(list(val=0))
  }
  alphai = alphas[i1]+y[i1]*y[j1]*(alphas[j1]-alphaj)
  
  
  #compute new b and assign
  b1 = b + Ei + y[i1]*(alphai-alphas[i1])*k[i1,i1] + y[j1]*(alphaj-alphas[j1])*k[i1,j1]
  b2 = b + Ej + y[i1]*(alphai-alphas[i1])*k[i1,j1] + y[j1]*(alphaj-alphas[j1])*k[j1,j1]
  # update b
  if (alphas[i1]>0 & alphas[i1]<C) {
    b = b1
  } else if (alphas[j1]>0 & alphas[j1]<C) {
    b = b2
  } else {
    b = (b1+b2)/2
  }
  
  #assign new alpha i and j to alphas
  alphas[i1] = alphai
  alphas[j1] = alphaj
  
  #update new weight vector
  #w = colSums(apply(x, 2, function(z) y*z*alphas))
  #update error
  ui = sum(alphas*y*k[i1,]) - b
  uj = sum(alphas*y*k[j1,]) - b
  #ui = x[i1,]%*%w - b
  #uj = x[j1,]%*%w - b
  error[i1] = ui - y[i1]
  error[j1] = uj - y[j1]
  result = list(val=1,alphas=alphas,error=error,   b=b)
  return(result)
}


examineExample = function(j,k,y,  b,alphas,eps,tol,C,error,n){
  yj = y[j]
  alphaj = alphas[j]
  #w = colSums(apply(x, 2, function(z) y*z*alphas))
  #uj = x[j,]%*%w - b
  uj = sum(alphas*y*k[j,])
  Ej = uj - yj
  rj = Ej*yj
  if ((rj < -tol & alphaj < C) | (rj > tol & alphaj > 0)) {
    # first heuristic: choose i that maximize abs(Ei -Ej), only select from suport vector
    #pos = which(alphas != 0 | alphas != C)
    pos = which(alphas>eps & alphas<(C-eps))
    
    #rest = which(!(alphas != 0 & alphas != C))
    rest = which(!(alphas>eps & alphas<(C-eps)))
    if (length(pos)>1) {
      # select i from unbound max abs(Ej - Ei)
      temp = matrix(0,nrow(alphas),1)
      temp[pos] = abs(as.numeric(Ej) - error[pos])
      i1 = which.max(temp)
      heu1 = takeStep(i1, j, k=k,y=y,   b=b,alphas=alphas,eps=eps,C=C,error)
      if (heu1$val) {
        #w = heu1$w
        alphas = heu1$alphas
        error = heu1$error
        b = heu1$b
        n = n+1
#         print(cat('Iteration',n,'  heu1','  b',b))
        result = list(val=1,alphas=alphas,error=error,b=b,n=n)
        return(result)
      }
    }
    
    # second heuristic: loop over i from all SV
    for (i1 in sample(pos,length(pos))) {
      if (alphas[i1]>eps & alphas[i1]<(C-eps)) {
        heu2 = takeStep(i1, j, k=k,y=y,b=b,alphas=alphas,eps=eps,C=C,error)
        if (heu2$val){
          #w = heu2$w
          alphas = heu2$alphas
          error = heu2$error
          b = heu2$b
          n = n+1
#           print(cat('Iteration',n,'  heu 2','  b',b))
          result = list(val=1,alphas=alphas,error=error    ,b=b,n=n)
          return(result)
        }
      }
      
    }
    # third heuristic: loop through the rest of data set
    #rest = which(!(alphas != 0 & alphas != C))
    rest = which(!(alphas>eps & alphas<(C-eps)))
    #print(cat('rest',length(rest)))
    for (i1 in sample(rest,length(rest))) {
      heu3 = takeStep(i1, j, k=k,y=y,    b=b,alphas=alphas,eps=eps,C=C,error)
      if (heu3$val) {
        #w = heu3$w
        alphas = heu3$alphas
        error = heu3$error
        b = heu3$b
        n = n+1
#         print(cat('Iteration',n,'  heu  3','  b',b))
        result = list(val=1,alphas=alphas,error=error,    b=b,n=n)
        return(result)
      }
    }
  }
  return(list(val=0,alphas=alphas,error=error,    b=b,n=n))
}


main = function(k,y,C,tol) {
  #w = as.numeric(0)
  n=0
  alphas = matrix(0,nrow(k),1) #initial alphas
  #w = colSums(apply(x, 2, function(z) y*z*alphas)) #in order to initial b
  b = matrix(0,1,1) #initial b
  error = matrix(0,nrow(k),1)
  eps = 1e-3
  numChanged = 0
  examineAll = 1
  while (numChanged > 0 | examineAll) {
    numChanged = 0
    if (examineAll) {
      for (i0 in 1:nrow(k)) {
        # loop through all data set
        examA = examineExample(i0,k=k,y=y,   b=b,alphas=alphas,eps=eps,tol=tol,C=C,error=error,n=n)
        numChanged = numChanged + examA$val
        #w = examA$w
        alphas = examA$alphas
        error = examA$error
        b = examA$b
        n = examA$n
      } 
    } else {
      unbound = which(alphas!=0 & alphas!=C)
      #unbound = which(alphas>eps & alphas<(C-eps))
      for (i0 in unbound) {
        # loop through only the unbound data
        examB = examineExample(i0,k=k,y=y,   b=b,alphas=alphas,eps=eps,tol=tol,C=C,error=error,n=n)
        numChanged = numChanged + examB$val
        #w = examB$w
        alphas = examB$alphas
        error = examB$error
        b = examB$b
        n = examB$n
      }
    }
    
    if (examineAll == 1) {
      examineAll = 0
    } else if (numChanged == 0) {
      examineAll = 1
    }
  }
  result = list(alphas=alphas,b=b)
  return(result)
}


objfun = function(a,ix,ai,y,k) {
  a[ix] = ai
  object = sum(a%*%t(a)*y%*%t(y)*k)/2-sum(a)
  return(object)
}



##############################################################################################################
##################################      Kernal function      #################################################
##############################################################################################################


# pre processing data
kernal_linear = function(x) {
  K = x%*%t(x)
  return(K)
}


kernal_RBF = function(x,sigma) {
  x2 = rowSums(x^2)
  temp = -2*(x%*%t(x))
  kkk = t(apply(temp, 1, function(z) z+x2))
  kk = apply(kkk, 2, function(z) z+x2)
  K = exp(-kk/(2*sigma^2))
  return(K)
}

##############################################################################################################
##################################      linear     plot      #################################################
##############################################################################################################


linearplot = function(X,Y,w,b,main) {
  x1 = X[Y==-1,]
  x2 = X[Y==1,]
  plot(x1, pch=1,col="blue",xlim=range(X[,1]),ylim=range(X[,2]),main=main,xlab='',ylab='')
  points(x2, pch=1,col="red")
  
  xplot = seq(min(X[,1]),max(X[,1]),(max(X[,1])-min(X[,1]))/100)
  yplot = (b - w[1]*xplot)/w[2]
  lines(xplot,yplot)
}

##############################################################################################################
##################################      RBF Kernal plot      #################################################
##############################################################################################################

Predict = function(x,svalpha,svX,svY,sigma=0.1) {
  #   indx = alpha>0
  x1 = rowSums(x^2)
  sv1 = rowSums(svX^2)
  temp = -2*x%*%t(svX)
  k1 = apply(temp,1,function(z) z+sv1)
  k = apply(t(k1),2,function(z) z+x1)
  k = exp(-k/(2*sigma^2))
  k1 = apply(k,1,function(z) z*svY)
  k2 = apply(t(k1),1,function(z) z*svalpha)
  K = t(k2)
  K = rowSums(K)
  K[K>=0] = 1
  K[K<0] =-1
  return(K)
}

gen_space = function(data,Y,alpha,l,sigma=0.1){
  idx = alpha>0
  x1plot = seq(min(data[,1]),max(data[,1]),(max(data[,1])-min(data[,1]))/(l-1))
  x2plot = seq(min(data[,2]),max(data[,2]),(max(data[,2])-min(data[,2]))/(l-1))
  xx1 = t(matrix(rep(x1plot,l),l,l))
  xx2 = matrix(rep(x2plot,l),l,l)
  
  svX2 = data[idx,]
  svY2 = Y[idx,]
  svalphas2 = alpha[idx,]
  vals = matrix(0,nrow(xx1),ncol(xx2))
  for(i in 1:ncol(xx1)) {
    xx = cbind(xx1[,i],xx2[,i])
    vals[i,] = Predict(xx,svalphas2,svX2,svY2,sigma=sigma)
  }
  return(list(x=x1plot,y=x2plot,z=vals))
  #   contour(x1plot,x2plot,vals,nlevels=nlevel,lty=lty)
}









