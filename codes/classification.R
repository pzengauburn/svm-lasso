######################################################################
# classification 
# 2025-09-24 
######################################################################

twoclass = function(x1, x2, lambda,
    ncores = 1, maxiter = 50, tol = 1e-3, mc.sample = 1e+5, ratio = 0.5)
{
    p = ncol(x1); n1 = nrow(x1); n2 = nrow(x2); 
    stopifnot(ncol(x2) == p, lambda > 0);

    ##################################################################
    # prepare data for analysis 
    ##################################################################

    x1mean = colMeans(x1);    x1sd = apply(x1, 2, sd);
    x2mean = colMeans(x2);    x2sd = apply(x2, 2, sd); 
    x_center = (x1mean + x2mean) / 2; 
    x1c = sweep(x1, 2, x_center, "-");  # mean = (x1mean - x2mean)/2
    x2c = sweep(x2, 2, x_center, "-");  # mean = (x2mean - x1mean)/2
    mu = (x1mean - x2mean) / 2; 

    # estimate sigma 
    S1 = compute_cov(x1);
    S2 = compute_cov(x2); 
    Sigma = ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2); 

    # compute the inverse of Sigma and get its diagonal 
    Sigma.inv = solve(Sigma);
    Sigma.jj = diag(Sigma.inv); 

    # compute saddle points 
    saddle = saddle.point(
                alpha = (n1+n2)/p, mu = mu, Sigma = Sigma, lambda = lambda,  
                ncores = ncores, loss = "SVM", penalty = "lasso",
                maxiter = maxiter, tol = tol, mc.sample = mc.sample,
                ratio = ratio,  trace.params = TRUE, details = FALSE); 
    
    # fit SVM
    y0 = c(rep(-1, nrow(x1)), rep(1, nrow(x2)));
    x0 = rbind(x1, x2); 
    fitsvm = svm.hinge.lasso.gurobi2(y0, x0, lambda);
    b.hat = fitsvm$b; 

    # compute the bias term 
    Sigma.inv.x = x0 %*% Sigma.inv; 
    u = y0 * (x0 %*% b.hat) / sqrt(p);
    Vprim = -fitsvm$gurobi$pi * sqrt(p);   
    b.adjust = colSums( Sigma.inv.x * drop(y0 * Vprim) );

    # compute the debiased b and its standard error 
    zeta = saddle$zeta;
    tau = sqrt(saddle$zeta.0) / saddle$zeta;
    b.debiased = b.hat - b.adjust / sqrt(p) / zeta;  
    b.se = tau * sqrt(Sigma.jj); 

    pvalue = 2 * pnorm(abs(b.debiased) / b.se, lower.tail = FALSE); 

    list(hat = b.hat, debiased = b.debiased, se = b.se, pvalue = pvalue, 
         saddle = saddle); 
}


compute_cov = function(x)
{
    cov(x) + diag(1e-5, nrow = ncol(x)); 
}


######################################################################
# THE END 
######################################################################
