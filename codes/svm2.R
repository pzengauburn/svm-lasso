######################################################################
# support vector machine 
# 2025-03-25 
# 2025-06-04 (rewrite svm.hinge.lasso.gurobi() to include sqrt(p))
######################################################################

######################################################################
# simulate samples from multivariate normal distribution 
# mean = mu and Sigma^{1/2} = Sigma.half
######################################################################

rmvnorm = function(n, mu, Sigma.half)
{
    stopifnot(length(mu) == nrow(Sigma.half), length(mu) == ncol(Sigma.half));
    p = length(mu); 
    sweep(matrix(rnorm(n * p), nrow = n) %*% Sigma.half, 2, mu, "+"); 
}

simulate.data = function(n, mu, Sigma)
{
    y = c(rep(1, n/2), rep(-1, n/2));  
    x = matrix(0, nrow = length(y), ncol = length(mu)); 
    Schol = chol(Sigma); 
    x[y ==  1, ] = rmvnorm(sum(y ==  1),  mu, Schol); 
    x[y == -1, ] = rmvnorm(sum(y == -1), -mu, Schol); 

    list(x = x, y = y); 
}

######################################################################
# hinge loss with lasso penalty using Gurobi 
######################################################################

svm.hinge.lasso.gurobi2 = function(y, x, lambda, ...)
{
    require(gurobi); 

    n = nrow(x); p = ncol(x); 
    stopifnot(length(y) == n, all(abs(y) == 1)); 
    stopifnot(length(lambda) == 1, lambda >= 0); 

    ##################################################################
    # prepare for linear programming 
    ##################################################################

    model = list(modelsense = "min"); 

    ##################################################################
    # min sum(xi) + lambda * sum(abs(w))
    # parameters are, w+, w-, xi 
    ##################################################################

    model$obj = c(rep(lambda, p + p), rep(1, n)); 
    model$lb = rep(0.0, p + p + n); 
    model$ub = rep(Inf, p + p + n); 

    ##################################################################
    # linear constraints  
    # yx'w >= 1 - xi 
    ##################################################################

    model$A = cbind(x * y, -x * y, diag(rep(sqrt(p), n))); 
    model$sense = rep(">=", n);
    model$rhs = rep(sqrt(p), n); 

    ##################################################################
    # linear constraints  
    ##################################################################

    params = list(OutputFlag = 0, ...); 
    LPfit = gurobi(model, params); 
    if(LPfit$status != "OPTIMAL")
        warning("\n\nThe status of gurobi(model, params) is ", LPfit$status, 
                "\n===== Something is wrong! =====\n"); 
    
    b = LPfit$x[1:p] - LPfit$x[(p+1):(p+p)];
    list(method = "gurobi", b = b, Qvalue = LPfit$objval, gurobi = LPfit);
}

######################################################################
# THE END 
######################################################################
