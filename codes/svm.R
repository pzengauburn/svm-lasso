######################################################################
# support vector machine 
# 2025-03-25 
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

svm.hinge.lasso.gurobi = function(y, x, lambda, ...)
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

    model$A = cbind(x * y, -x * y, diag(rep(1, n))); 
    model$sense = rep(">=", n);
    model$rhs = rep(1, n); 

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
# hinge loss with lasso penalty using CVXR 
######################################################################

svm.hinge.lasso.cvx = function(y, x, lambda, ...)
{
    require(CVXR); 

    n = nrow(x); p = ncol(x); 
    stopifnot(length(y) == n, all(abs(y) == 1)); 
    stopifnot(length(lambda) == 1, lambda >= 0); 

    ############################################################
    # define the optimization problem  
    ############################################################

    bvar = Variable(p); 
    obj = sum(pos(1 - y * (x %*% bvar))) + lambda * p_norm(bvar, 1); 

    prob = Problem(Minimize(obj)); 
    result = solve(prob, ...);
    b = drop( result$getValue(bvar) ); 

    ############################################################
    # output 
    ############################################################

    list(method = "CVXR", b = b, Qvalue = result$value, cvx = result); 
}

######################################################################
# hinge loss with ridge penalty using Gurobi 
######################################################################

svm.hinge.ridge.gurobi = function(y, x, lambda, ...)
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
    # min sum(xi) + lambda * sum(w^2)
    # parameters are, w+, w-, xi 
    ##################################################################

    model$Q = diag(c(rep(lambda, p + p), rep(0, n))); 
    model$obj = c(rep(0, p + p), rep(1, n)); 
    model$lb = rep(0.0, p + p + n); 
    model$ub = rep(Inf, p + p + n); 

    ##################################################################
    # linear constraints  
    # yx'w >= 1 - xi 
    ##################################################################

    model$A = cbind(x * y, -x * y, diag(rep(1, n))); 
    model$sense = rep(">=", n);
    model$rhs = rep(1, n); 

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
# hinge loss with Ridge penalty using CVXR 
######################################################################

svm.hinge.ridge.cvx = function(y, x, lambda, ...)
{
    require(CVXR); 

    n = nrow(x); p = ncol(x); 
    stopifnot(length(y) == n, all(abs(y) == 1)); 
    stopifnot(length(lambda) == 1, lambda >= 0); 

    ############################################################
    # define the optimization problem  
    ############################################################

    bvar = Variable(p); 
    obj = sum(pos(1 - y * (x %*% bvar))) + lambda * p_norm(bvar, 2)^2; 

    prob = Problem(Minimize(obj)); 
    result = solve(prob, ...);
    b = drop( result$getValue(bvar) ); 

    ############################################################
    # output 
    ############################################################

    list(method = "CVXR", b = b, Qvalue = result$value, cvx = result); 
}

######################################################################
# THE END 
######################################################################

fit.svm.e1071 = function(x, y, lambda) 
{
    require(e1071); 
    stopifnot(all(y %in% c(-1, 1))); 
    stopifnot(length(y) == NROW(x), lambda >= 0); 

    n = NROW(x); p = NCOL(x); 

    ############################################################
    # fit model using e1071 
    ############################################################

    fitsvm = svm(x, y, scale = FALSE, type = "C-classification",
                  cost = 0.5 / lambda, kernel = "linear", 
                  cross = 0, probability = FALSE, tolerance = 1e-5); 
    b = coef(fitsvm); 

    ############################################################
    # output 
    ############################################################

    list(method = "kernlab", lambda = lambda, a = b[1], b = b[-1],
         kernlab = fitsvm); 
}


fit.svm.kernlab = function(x, y, lambda, ...) 
{
    require(kernlab); 
    stopifnot(all(y %in% c(-1, 1))); 
    stopifnot(length(y) == NROW(x), lambda >= 0); 

    n = NROW(x); p = NCOL(x); 

    ############################################################
    # fit model using kernlab 
    # scaled = FALSE --- do not center and scale predictors 
    # C              --- 
    # C-svc          --- C classification 
    # vanilladot     --- linear kernal 
    # cross = 0      --- no cross-validation 
    # prob.model = F --- do not compute class probabilities 
    ############################################################

    fitsvm = ksvm(x, y, scaled = FALSE, type = "C-svc",
                  C = lambda, kernel = "vanilladot", 
                  cross = 0, prob.model = FALSE); 
    b = coef(fitsvm); 

    ############################################################
    # output 
    ############################################################

    list(method = "kernlab", lambda = lambda, b = b,
         kernlab = fitsvm); 
}

