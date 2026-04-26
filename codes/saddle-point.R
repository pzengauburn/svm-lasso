######################################################################
# compute saddle-point parameters 
# 2023-04-28 (add option penalty = "exact-lasso")
# 2023-05-05 (fix a bug for ini)
# 2023-07-06 (add option penalty = "phi-lasso") --- version 4
# 2025-02-19 (add argument ratio = 1.0) --- version 5
# 2025-12-11 (add argument mu.length, function adjusted_mu) -- version 6 
######################################################################

library(parallel); 

# dyn.load(paste0("lassoCD", .Platform$dynlib.ext));

saddle.point = function(alpha, mu, Sigma, lambda, mu.length = NULL, 
    ini = list(), ncores = 1, 
    loss = c("logistic", "SVM"), 
    penalty = c("lasso", "ridge", "exact-ridge", "exact-lasso", "phi-lasso"), 
    maxiter = 100, tol = 1e-3, mc.sample = 10000, int.limit = Inf,
    ratio = 1.0, 
    trace.params = FALSE, details = FALSE)
{
    params = list(); 
    params$q  = ifelse(is.null(ini$q),  1.0, ini$q); 
    params$R  = ifelse(is.null(ini$R),  1.0, ini$R); 
    params$q0 = ifelse(is.null(ini$q0), 1.0, ini$q0); 
    params$zeta  = ifelse(is.null(ini$zeta),  1.0, ini$zeta); 
    params$zeta.0  = ifelse(is.null(ini$zeta.0),  1.0, ini$zeta.0); 
    params$R0 = ifelse(is.null(ini$R0), 1.0, ini$R0); 

    loss = match.arg(loss); 
    penalty = match.arg(penalty); 

    p = length(mu);                             # number of variables 

    if(is.null(mu.length))
        mu.length = sqrt(sum(mu * mu));         # length of mu

    mu.direction = mu / sqrt(sum(mu * mu));     # direction of mu

    Sigma.inv.mu = solve(Sigma, mu.direction);  # solve(Sigma) %*% mu.hat
    muSmu = sum(mu.direction * Sigma.inv.mu);

    Sigma.eig = eigen(Sigma); 
    Sigma.half = (Sigma.eig$vectors %*% 
                  diag( sqrt(Sigma.eig$values) ) %*% 
                  t(Sigma.eig$vectors));        # Sigma^{1/2}

    if(details)
    {
        cat("The number of cores used is ", ncores, ".\n", sep = "");
        cat("alpha =", alpha, "lambda =", lambda, "\n"); 
    }

    if((penalty == "exact-lasso") || (penalty == "phi-lasso"))
    {
        sigma2 = Sigma[1];
        if( max(abs(Sigma - diag(rep(sigma2, p)))) > 1e-10 )
            stop(paste0('penalty = "', penalty,'" requires Sigma is diagonal.\n'));
    }

    if(trace.params) params.trace = list(); 

    iter = 0; conv = tol + 1.0; 
    while((iter < maxiter) && (conv > tol))
    {
        iter = iter + 1; 
        params.old = params;

        ##############################################################
        # compute zeta.0, zeta, R0
        ##############################################################

        params.tmp = params; 
        params = update.zeta0zetaR0(params, alpha, mu.length, 
                                    loss = loss, int.limit = int.limit); 

        # damping 
        params$zeta.0 = ratio * params$zeta.0 + (1.0 - ratio) * params.tmp$zeta.0;
        params$zeta = ratio * params$zeta + (1.0 - ratio) * params.tmp$zeta;
        params$R0 = ratio * params$R0 + (1.0 - ratio) * params.tmp$R0;

        if(details)
        {
            cat(">>>>> iter =", iter, "<<<<<\n");
            with(params, cat("zeta.0 =", zeta.0, "zeta =", zeta, "R0 =", R0, "\n")); 
        }

        ##############################################################
        # compute q0, q, R
        ##############################################################

        params.tmp = params; 
        if(penalty == "exact-ridge")
        {
            params = update.q0qR.ridge(params, mu.direction, Sigma, lambda,
                            muSmu = muSmu);
        }
        else if(penalty == "exact-lasso")
        {
            params = update.q0qR.lasso(params, mu.direction, sigma2, lambda);
        }
        else if(penalty == "phi-lasso")
        {
            params = update.q0qR.phi(params, mu.direction, sigma2, lambda);
        }
        else
        { 
            params = update.q0qR(params, mu.direction, Sigma, lambda,  
                             muSmu = muSmu, Sigma.half = Sigma.half, 
                             penalty = penalty, mc.sample = mc.sample, 
                             ncores = ncores);
        }

        # damping 
        params$q0 = ratio * params$q0 + (1.0 - ratio) * params.tmp$q0;
        params$q = ratio * params$q + (1.0 - ratio) * params.tmp$q;
        params$R = ratio * params$R + (1.0 - ratio) * params.tmp$R;

        if(details)
            with(params, cat("q0 =", q0, "q =", q, "R =", R, "\n")); 

        if(trace.params) params.trace[[iter]] = params; 

        ##############################################################
        # check convergence 
        ##############################################################

        conv = max(abs(unlist(params[c("q", "q0", "R")]) - unlist(params.old[c("q", "q0", "R")]))); 
    }

    trace = NULL; 
    if(trace.params) 
    {
        trace = do.call(rbind, lapply(params.trace, data.frame)); 
    }

    c(list(alpha = alpha, lambda = lambda, 
            iter = iter, conv = conv, trace = trace), 
      params);
}

update.zeta0zetaR0 = function(params, alpha, mu.length, 
    loss = c("logistic", "SVM"), int.limit = Inf)
{
    loss = match.arg(loss); 

    ##################################################################
    # compute E(u.hat^2), E(u.hat * e), E(u.hat) 
    ##################################################################

    Rmu = params$R * mu.length; 
    sqrt.q0 = sqrt(params$q0); 

    if(loss == "logistic")
    {
        u.hat.fun = compute.u.hat.logistic; 
    }
    else if(loss == "SVM")
    {
        u.hat.fun = compute.u.hat.SVM; 
    }

    E.uhat.2 = integrate(function(e) 
                        { 
                            u.hat.fun(Rmu + sqrt.q0 * e, params$q)^2 * dnorm(e)
                        }, 
                        -int.limit, int.limit)$value; 

    E.uhat.e = integrate(function(e) 
                        { 
                            u.hat.fun(Rmu + sqrt.q0 * e, params$q) * e * dnorm(e)
                        }, 
                        -int.limit, int.limit)$value; 

    E.uhat = integrate(function(e) 
                        { 
                            u.hat.fun(Rmu + sqrt.q0 * e, params$q) * dnorm(e)
                        }, 
                        -int.limit, int.limit)$value; 

    ##################################################################
    # compute zeta_0 
    ##################################################################

    params$zeta.0 = alpha / params$q^2 * 
                (E.uhat.2 + Rmu * Rmu + params$q0 
                - 2 * Rmu * E.uhat - 2 * sqrt.q0 * E.uhat.e); 

    ##################################################################
    # compute zeta 
    ##################################################################

    params$zeta = -alpha / params$q / sqrt.q0 * (E.uhat.e - sqrt.q0); 

    ##################################################################
    # compute R0 
    #################################################################

    params$R0 = alpha * mu.length / params$q * (E.uhat - Rmu); 

    params; 
}

update.zeta0zetaR0.v0 = function(params, alpha, mu.length, 
    int.limit = 100)
{
    ##################################################################
    # compute zeta_0 
    ##################################################################

    exp.uhat.2.fun = function(eps, Rmu, sqrt.q0, q)
    {
        ans = numeric(length(eps)); 
        for(i in 1:length(eps))
        {
            u.hat = compute.u.hat.logistic(Rmu + sqrt.q0 * eps[i], q)
            ans[i] = (u.hat - Rmu - sqrt.q0 * eps[i])^2 * dnorm(eps[i]);
        }
             
        ans; 
    } 

    exp.uhat.2 = integrate(exp.uhat.2.fun, -int.limit, int.limit, 
                           Rmu = params$R * mu.length, 
                           sqrt.q0 = sqrt(params$q0), 
                           q = params$q); 
    params$zeta.0 = alpha / params$q^2 * exp.uhat.2$value; 

    ##################################################################
    # compute zeta 
    ##################################################################

    exp.uhat.1.fun = function(eps, Rmu, sqrt.q0, q)
    {
        ans = numeric(length(eps)); 
        for(i in 1:length(eps))
        {
            u.hat = compute.u.hat.logistic(Rmu + sqrt.q0 * eps[i], q)
            ans[i] = (u.hat - Rmu - sqrt.q0 * eps[i]) * eps[i] * dnorm(eps[i]);
        }
             
        ans; 
    } 

    exp.uhat.1 = integrate(exp.uhat.1.fun, -int.limit, int.limit,
                           Rmu = params$R * mu.length, 
                           sqrt.q0 = sqrt(params$q0),
                           q = params$q); 
    params$zeta = -alpha / params$q / sqrt(params$q0) * exp.uhat.1$value; 

    ##################################################################
    # compute R0 
    ##################################################################

    exp.uhat.0.fun = function(eps, Rmu, sqrt.q0, q)
    {
        ans = numeric(length(eps)); 
        for(i in 1:length(eps))
        {
            u.hat = compute.u.hat.logistic(Rmu + sqrt.q0 * eps[i], q)
            ans[i] = (u.hat - Rmu) * dnorm(eps[i]);
        }
             
        ans; 
    } 

    exp.uhat.0 = integrate(exp.uhat.0.fun, -int.limit, int.limit,
                           Rmu = params$R * mu.length, 
                           sqrt.q0 = sqrt(params$q0),
                           q = params$q); 
    params$R0 = alpha * mu.length / params$q * exp.uhat.0$value; 

    params; 
}

update.q0qR.ridge = function(params, mu.direction, Sigma, lambda,
    muSmu = NULL)
{
    p = length(mu.direction);                   # number of variables 

    if(is.null(muSmu))
    {
        muSmu = sum(mu.direction * solve(Sigma, mu.direction)); 
    }

    SigmaI.inv = solve(Sigma + diag(rep(2 * lambda / params$zeta, p))); 

    ##################################################################
    # compute w.hat for one z 
    ##################################################################

    tau = sqrt(params$zeta.0) / params$zeta;
    R0zeta = params$R0 / params$zeta; 

    SSS = SigmaI.inv %*% Sigma %*% SigmaI.inv; 
    E.wSw = tau^2 * sum(diag( SSS %*% Sigma )) + 
            p * R0zeta^2 * sum(mu.direction * (SSS %*% mu.direction));
    E.wSz = tau * sum(diag(SigmaI.inv %*% Sigma));
    E.wmu = sqrt(p) * R0zeta * sum(mu.direction * (SigmaI.inv %*% mu.direction));

    ##################################################################
    # compute q_0, q, R as sample average 
    ##################################################################

    params$q0 = E.wSw / p;
    params$q = E.wSz / p / sqrt(params$zeta.0);
    params$R = E.wmu / sqrt(p);

    # c(params, list(wSw = E.wSw, wSz = E.wSz, wmu = E.wmu)); 
    params; 
}


update.q0qR.lasso = function(params, mu.direction, sigma2, lambda)
{
    p = length(mu.direction);                   # number of variables 
    # Sigma = diag(rep(sigma2, p));

    ##################################################################
    # compute w.hat for one z 
    ##################################################################

    tau = sqrt(params$zeta.0) / params$zeta;
    lambda.zeta.sigma2 = lambda / params$zeta / sigma2; 
    tau.sigma = tau / sqrt(sigma2); 
    const = sqrt(p) * params$R0 / params$zeta / sigma2; 

    E.wSw = 0; 
    E.wSz = 0; 
    E.wmu = 0; 

    muk = unique(mu.direction);

    for(k in 1:length(muk))
    {
        off = const * muk[k]; 
        neg = -(lambda.zeta.sigma2 + const * muk[k]) / tau.sigma;
        pos =  (lambda.zeta.sigma2 - const * muk[k]) / tau.sigma;
        E.wSw = E.wSw + sum(muk[k] == mu.direction) * (
                    integrate(function(z)
                        { (off + tau.sigma * z + lambda.zeta.sigma2)^2 * dnorm(z) }, 
                        -Inf, neg)$value 
                    +
                    integrate(function(z)
                        { (off + tau.sigma * z - lambda.zeta.sigma2)^2 * dnorm(z) },
                        pos, Inf)$value); 
        E.wSz = E.wSz + sum(muk[k] == mu.direction) * (
                    integrate(function(z)
                        { (off + tau.sigma * z + lambda.zeta.sigma2) * z * dnorm(z) }, 
                        -Inf, neg)$value 
                    + 
                    integrate(function(z)
                        { (off + tau.sigma * z - lambda.zeta.sigma2) * z * dnorm(z) }, 
                        pos, Inf)$value); 
        E.wmu = E.wmu +  sum(muk[k] == mu.direction) * muk[k] * (
                    integrate(function(z)
                        { (off + tau.sigma * z + lambda.zeta.sigma2) * dnorm(z) }, 
                        -Inf, neg)$value 
                    +
                    integrate(function(z)
                        { (off + tau.sigma * z - lambda.zeta.sigma2) * dnorm(z) },
                        pos, Inf)$value); 
    }

    E.wSw = E.wSw * sigma2;
    E.wSz = E.wSz * sqrt(sigma2);   

    ##################################################################
    # compute q_0, q, R as sample average 
    ##################################################################

    params$q0 = E.wSw / p;
    params$q = E.wSz / p / sqrt(params$zeta.0);
    params$R = E.wmu / sqrt(p);

    # c(params, list(wSw = mean.wSw, wSz = mean.wSz, wmu = mean.wmu)); 
    params; 
}


update.q0qR.phi = function(params, mu.direction, sigma2, lambda)
{
    p = length(mu.direction);                   # number of variables 
    # Sigma = diag(rep(sigma2, p));

    ##################################################################
    # compute w.hat for one z 
    ##################################################################

    tau = sqrt(params$zeta.0) / params$zeta;
    lambda.zeta.sigma2 = lambda / params$zeta / sigma2; 
    tau.sigma = tau / sqrt(sigma2); 
    const = sqrt(p) * params$R0 / params$zeta / sigma2; 

    E.wSw = 0; 
    E.wSz = 0; 
    E.wmu = 0; 

    muk = unique(mu.direction);

    for(k in 1:length(muk))
    {
        off = const * muk[k]; 
        ck = off - lambda.zeta.sigma2;
        dk = off + lambda.zeta.sigma2;
        Lk = -(lambda.zeta.sigma2 + const * muk[k]) / tau.sigma;
        Uk =  (lambda.zeta.sigma2 - const * muk[k]) / tau.sigma;

        PhiLk = pnorm(Lk);     phi_Lk = dnorm(Lk);
        PhiUk = 1 - pnorm(Uk); phi_Uk = dnorm(Uk); 

        E.wSw = E.wSw + sum(muk[k] == mu.direction) * (
                    (dk * dk * sigma2 + tau * tau) * PhiLk 
                    - dk * tau * sqrt(sigma2) * phi_Lk 
                    + (ck * ck * sigma2 + tau * tau) * PhiUk 
                    + ck * tau * sqrt(sigma2) * phi_Uk); 
        E.wSz = E.wSz + sum(muk[k] == mu.direction) * (
                    PhiLk + PhiUk); 
        E.wmu = E.wmu +  sum(muk[k] == mu.direction) * muk[k] * (
                    dk * PhiLk - tau.sigma * phi_Lk  
                    + ck * PhiUk + tau.sigma * phi_Uk); 
    }

    E.wSz = E.wSz * tau;   

    ##################################################################
    # compute q_0, q, R as sample average 
    ##################################################################

    params$q0 = E.wSw / p;
    params$q = E.wSz / p / sqrt(params$zeta.0);
    params$R = E.wmu / sqrt(p);

    # c(params, list(wSw = mean.wSw, wSz = mean.wSz, wmu = mean.wmu)); 
    params; 
}

update.q0qR = function(params, mu.direction, Sigma, lambda,  
    muSmu = NULL, Sigma.half = NULL, 
    penalty = c("lasso", "ridge"), mc.sample = 1000, ncores = 2)
{
    penalty = match.arg(penalty); 

    p = length(mu.direction);                   # number of variables 

    if(is.null(muSmu)) 
    {
        Sigma.inv.mu = solve(Sigma, mu.direction);  # solve(Sigma) %*% mu.hat
        muSmu = sum(mu.direction * Sigma.inv.mu);   # t(mu.hat) %*% Solve(Sigma) %*% mu.hat 
    }

    if(is.null(Sigma.half))
    {
        Sigma.eig = eigen(Sigma); 
        Sigma.half = (Sigma.eig$vectors %*% 
                      diag( sqrt(Sigma.eig$values) ) %*% 
                      t(Sigma.eig$vectors));        # Sigma^{1/2}

    }

    ##################################################################
    # compute w.hat for one z 
    ##################################################################

    tau = sqrt(params$zeta.0) / params$zeta;

    if(penalty == "lasso")
    {
        w.hat.fun = compute.w.hat.lasso; 
    }
    else if(penalty == "ridge")
    {
        w.hat.fun = compute.w.hat.ridge; 
    }

    one.case = function(index)
    {
        z = rnorm(p); 
        Sigma.half.z = Sigma.half %*% z;
        w.hat = w.hat.fun(Sigma, 
                    sqrt(p) * params$R0 / params$zeta * mu.direction + tau * Sigma.half.z, 
                    lambda / params$zeta); 
        c("index" = index,
          "wSw" = sum(w.hat * (Sigma %*% w.hat)), 
          "wSz" = sum(w.hat * Sigma.half.z), 
          "wmu" = sum(w.hat * mu.direction));
    }

    results = mclapply(1:mc.sample, one.case, mc.cores = ncores); 
    exps = do.call(rbind, results);
    mean.wSw = mean(exps[, "wSw"]);
    mean.wSz = mean(exps[, "wSz"]);
    mean.wmu = mean(exps[, "wmu"]);

    ##################################################################
    # compute q_0, q, R as sample average 
    ##################################################################

    params$q0 = mean.wSw / p;
    params$q = mean.wSz / p / sqrt(params$zeta.0);
    params$R = mean.wmu / sqrt(p);

    # c(params, list(wSw = mean.wSw, wSz = mean.wSz, wmu = mean.wmu)); 
    params; 
}

######################################################################
# compute u.hat
#     arg.min V(u) + (u - off.set)^2 / (2 * q)
# where Rmuq0e = R * mu + sqrt(q0) * e
# Note that R, mu, q0 are scalar, while e can be a vector  
#
# for logistic regression: V(u) = log(1 + exp(-u)) 
# for SVM: V(u) = max(0, 1 - u) 
######################################################################

compute.u.hat.logistic.V0 = function(Rmuq0e, q)
{
    ans = numeric(length(Rmuq0e)); 
    for(i in 1:length(Rmuq0e))
    {
        ans[i] = uniroot(function(u) 
                         {
                            (u - Rmuq0e[i]) * (exp(u) + 1) / q - 1 
                         }, 
                         c(Rmuq0e[i], q + Rmuq0e[i]))$root; 
    }

    ans; 
}

compute.u.hat.logistic = function(Rmuq0e, q)
{
    ans = numeric(length(Rmuq0e)); 
    for(i in 1:length(Rmuq0e))
    {
        ans[i] = ifelse(Rmuq0e[i] <= -0.5 * q, 
                        uniroot(function(u) 
                        {
                            (u - Rmuq0e[i]) * (exp(u) + 1) - q 
                        }, 
                        c(Rmuq0e[i] - 1, min(0, q + Rmuq0e[i]) + 1))$root,
                        uniroot(function(u)
                        {
                            (u - Rmuq0e[i]) * (exp(-u) + 1) - q * exp(-u)
                        },
                        c(max(0, Rmuq0e[i]) - 1, q + Rmuq0e[i] + 1))$root);
    }

    ans; 
}

compute.u.hat.SVM = function(Rmuq0e, q)
{
    ans = rep(1.0, length(Rmuq0e));      # default value is 1

    index = (Rmuq0e > 1); 
    ans[index] = Rmuq0e[index]; 
    index = (Rmuq0e < 1 - q); 
    ans[index] = q + Rmuq0e[index]; 

    ans;
}

######################################################################
# compute w.hat
#     arg.min 0.5 * zeta * |w - R0Smuz|_Sigma + lambda * sum(abs(w))
# where 
#     R0zetamuSz = sqrt(p) * R0 / zeta * mu.hat + tau * Sigma^{1/2} * z
# which is equivalent to 
#     arg.min 0.5 * w^T Sigma w - w^T R0zetamuSz + sum(abs(w)) * lambda.zeta 
# where lambda.zeta = lambda / zeta 
######################################################################

compute.w.hat.lasso = function(Sigma, R0zetamuSz, lambda.zeta, 
    maxiter = 100, tol = 1e-7)
{
    p = length(R0zetamuSz); 

    # it is a lasso problem 
    ans = .C("lasso_CD_core", 
            as.integer(p), 
            as.double(Sigma), 
            as.double(R0zetamuSz), 
            as.double(lambda.zeta), 
            as.integer(maxiter), 
            as.double(tol),
            x = double(p), 
            iter = integer(1), 
            conv = double(1));
    ans$x; 
}

######################################################################
# compute w.hat
#     arg.min 0.5 * zeta * |w - R0Smuz|_Sigma + lambda * sum(w^2)
# where 
#     R0zetamuSz = sqrt(p) * R0 /zeta * mu.hat + tau * Sigma^{1/2} * z
# which is equivalent to 
#     arg.min 0.5 * w^T Sigma w - w^T R0zetamuSz + sum(w^2) * lambda.zeta 
# where lambda.zeta = lambda / zeta 
######################################################################

compute.w.hat.ridge = function(Sigma, R0zetamuSz, lambda.zeta, 
    maxiter = 100, tol = 1e-7)
{
    mat = Sigma; 
    diag(mat) = diag(mat) + 2 * lambda.zeta; 
    solve(mat, R0zetamuSz);
}

######################################################################
## find minimizer of
## Q(x) = 0.5 * x^T S x - a^T x + lambda * sum(abs(x))
######################################################################

C.lassoCD.core = function(Smat, avec, lambda, maxiter = 100, tol = 1e-7)
{
    p = length(avec);
    stopifnot(nrow(Smat) == p, ncol(Smat) == p); 
    stopifnot(all(diag(Smat)) > 0); 
    stopifnot(lambda >= 0); 

    ans = .C("lasso_CD_core", 
            as.integer(p), 
            as.double(Smat), 
            as.double(avec), 
            as.double(lambda), 
            as.integer(maxiter), 
            as.double(tol),
            x = double(p), 
            iter = integer(1), 
            conv = double(1))

  list(lambda = lambda, x = ans$x, iter = ans$iter, conv = ans$conv); 
}

######################################################################
# adjusted length of mu
######################################################################

adjusted_mu_length = function(xplus, xminus, size = 10000)
{
    stopifnot(ncol(xplus) == ncol(xminus))

    n_plus = nrow(xplus); n_minus = nrow(xminus); 
    r_plus  = n_plus  / (n_plus + n_minus);
    r_minus = n_minus / (n_plus + n_minus); 

    p = ncol(xplus); alpha = (n_plus + n_minus) / p; 

    mean_plus  = colMeans(xplus); 
    mean_minus = colMeans(xminus); 
    mu_c = mean_plus - mean_minus; 

    z_plus  = scale(xplus,  center = apply(xplus,  2, median));
    z_minus = scale(xminus, center = apply(xminus, 2, median));
    mad_x = median(abs(c(z_plus, z_minus)));
    z = rnorm(size);  
    mad_norm = median(abs(z - median(z)));
    sigma = mad_x / mad_norm; 

    0.5 * sqrt(sum(mu_c * mu_c) - sigma^2 / (alpha * r_plus * r_minus)); 
}

######################################################################
# THE END
######################################################################
