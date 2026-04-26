##########

rm(list = ls())
library("lars")
source("lassoCD.R")
dyn.load("lassoCD-clib.dll")

#################################################################
#################################################################

lars.wlasso = function(y, x, lambda = 0,
                       w = rep(1, length(y)), intercept = T)
{
  n = NROW(x); p = NCOL(x);
  if(n != length(y)) stop("length(y) != NROW(x).\n")

  w = drop(w); w2 = sqrt(w)
  if(intercept == F)
  {
    yt = y * w2; xt = x * w2; int = NULL
  }
  else  ## intercept == T
  {
    wsum = sum(w);
    ymu = sum(y * w) / wsum; yt = (y - ymu) * w2;
    if(p > 1)
    {
      xmu = apply(x * w, 2, sum) / wsum
      xt = sweep(x, 2, xmu, "-") * w2
    }
    else { xmu = sum(x * w) / wsum; xt = (x - xmu) * w2 }
  }

  fit = lars(xt, yt, type = "lasso", normalize = F, intercept = F)
  b = coef.lars(fit, s = lambda, mode = "lambda")

  if(intercept == T) int = ymu - sum(xmu * b)
  ans = c(int, b)
  list(beta = ans, lambda = lambda, fit = fit)
}

#################################################################
#################################################################

x = matrix(runif(8 * 100, -10, 10), ncol = 8)
y = runif(100, -10, 10)
w = runif(length(y), 1, 3)
lda = runif(1, 0, 10)

ansc = C.lassoCD(y, x, w = w, lambda = lda, maxiter = 2000)
ansr = lassoCD(y, x, w = w, lambda = lda, maxiter = 2000)
ans0 = lars.wlasso(y, x, lambda = lda, w = w)
c(sum(abs(ansc$b - ansr$b)), sum(abs(ansc$b - ans0$beta[-1])))

#################################################################

rm(list = ls())
source("lassoCD.R")
dyn.load("lassoCD-clib.dll")

x = matrix(runif(8 * 100, -10, 10), ncol = 8)
y = runif(100, -10, 10)
w = runif(length(y), 1, 3)
lda = runif(8, 0, 10)

W = diag(lda)
temp = calSmatavec(y, x, w, NULL, NULL)
ans1 = C.lassoCD.core(solve(W) %*% temp$Smat %*% solve(W), solve(W) %*% temp$avec, 1)
ans2 = C.lassoCD.core2(temp$Smat, temp$avec, lda)
max(abs(solve(W) %*% ans1$x - ans2$x))

#################################################################


##########
##########

data(diabetes, package = "lars")
y = diabetes$y
x = diabetes$x


lambda = runif(1, 0, 1000)
a = wlasso.coef(ans1, lambda)
b = coef(ans2, s = lambda, mode = "lambda")
sum(abs(a - b))

lambda = runif(1, 0, 1000)
a0 = wlasso.pred(y, x, ans1, lambda)
b0 = predict(ans2, x, s = lambda, mode = "lambda")
sum(abs(a0 - b0$fit))
 
##########
