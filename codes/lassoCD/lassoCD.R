######################################################################
## Coordinate Descent Algorithm for lasso 
######################################################################

calSmatavec = function(y, x, w, y0 = NULL, x0 = NULL)
{
  if(is.null(y0) | is.null(x0))
  {
    wsum = sum(w)
    y0 = sum(y * w) / wsum
    x0 = colSums(x * drop(w)) / wsum
  }

  xctd = scale(x, center = x0, scale = F)
  xw = xctd * drop(w)
  Smat = t(xw) %*% xctd
  avec = t(xw) %*% (y - y0)
  list(x0 = x0, y0 = y0, Smat = Smat, avec = avec)
}


lassoCD = function(y, x, w = rep(1, length(y)), lambda = 0,
                   maxiter = 100, tol = 1e-7)
{
  temp = calSmatavec(y, x, w, NULL, NULL)
  ans = lassoCD.core(temp$Smat, temp$avec, lambda, maxiter, tol)
  list(lambda = lambda, beta = ans$x, iter = ans$iter, conv = ans$conv,
       Smat = temp$Smat, avec = temp$avec)
}


C.lassoCD = function(y, x, w = rep(1, length(y)), lambda = 0,
                   maxiter = 100, tol = 1e-7)
{
  temp = calSmatavec(y, x, w, NULL, NULL)
  ans = C.lassoCD.core(temp$Smat, temp$avec, lambda, maxiter, tol)
  list(lambda = lambda, beta = ans$x, iter = ans$iter, conv = ans$conv,
       Smat = temp$Smat, avec = temp$avec)
}

C.lassoCD2 = function(y, x, w = rep(1, length(y)), lambda = 0,
                   maxiter = 100, tol = 1e-7)
{
  temp = calSmatavec(y, x, w, NULL, NULL)
  ans = C.lassoCD.core2(temp$Smat, temp$avec, lambda, maxiter, tol)
  list(lambda = lambda, beta = ans$x, iter = ans$iter, conv = ans$conv,
       Smat = temp$Smat, avec = temp$avec)
}

######################################################################
## find minimizer of
## Q(x) = x^T S x /2 - a^T x + lambda * sum(abs(x))
######################################################################

lassoCD.core = function(Smat, avec, lambda, maxiter = 100, tol = 1e-7)
{
  p = length(avec)
  if((NROW(Smat) != p) | (NCOL(Smat) != p))
    stop("NROW(Smat) != p or NCOL(Smat) != p.\n")
  if(any(diag(Smat) < 0))
    stop("Smat is not positive definite.\n")
  if(lambda < 0)
    stop("lambda is negative.\n")

  iter = 0; conv = tol + 1; newx = rep(0, p)
  while((iter < maxiter) & (conv > tol))
  {
    iter = iter + 1
    oldx = newx
    for(k in 1:p) newx = lassoCD.k(newx, Smat, avec, lambda, k)
    conv = max(abs(newx - oldx))
  }

  list(lambda = lambda, x = newx, iter = iter, conv = conv)
}


C.lassoCD.core = function(Smat, avec, lambda, maxiter = 100, tol = 1e-7)
{
  p = length(avec)
  if((NROW(Smat) != p) | (NCOL(Smat) != p))
    stop("NROW(Smat) != p or NCOL(Smat) != p.\n")
  if(any(diag(Smat) < 0))
    stop("Smat is not positive definite.\n")
  if(lambda < 0)
    stop("lambda is negative.\n")

  ans = .C("lasso_CD_core", as.integer(p), 
     as.double(Smat), as.double(avec), as.double(lambda), 
     as.integer(maxiter), as.double(tol),
     x = double(p), iter = integer(1), conv = double(1))

  list(lambda = lambda, x = ans$x, iter = ans$iter, conv = ans$conv)
}

C.lassoCD.core2 = function(Smat, avec, lambda, maxiter = 100, tol = 1e-7)
{
  p = length(avec)
  if((NROW(Smat) != p) | (NCOL(Smat) != p))
    stop("NROW(Smat) != p or NCOL(Smat) != p.\n")
  if(any(diag(Smat) < 0))
    stop("Smat is not positive definite.\n")
  if(any(lambda < 0))
    stop("lambda is negative.\n")

  ans = .C("lasso_CD_core2", as.integer(p), 
     as.double(Smat), as.double(avec), as.double(lambda), 
     as.integer(maxiter), as.double(tol),
     x = double(p), iter = integer(1), conv = double(1))

  list(lambda = lambda, x = ans$x, iter = ans$iter, conv = ans$conv)
}


######################################################################
## calculate sign(a) * (abs(a) - lambda)_+ / s
######################################################################

unilassoCD = function(s, a, lambda)
{
  if(a > lambda)
    ans = (a - lambda) / s
  else if(a < -lambda)
    ans = (a + lambda) / s
  else
    ans = 0.0  

  ans 
}

lassoCD.k = function(x, Smat, avec, lambda, k)
{
  ans = x; ans[k] = 0;
  s = Smat[k, k]
  a = avec[k] - sum(Smat[, k] * ans)
  ans[k] = unilassoCD(s, a, lambda)
  ans
}

######################################################################
## END OF THE FILE 
######################################################################
