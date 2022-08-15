function (y, burn = 300, MHsim = 1000, m_max = 5, p_max = 5, 
    p_fix = NULL, m_fix = NULL, y0 = NULL, p_init = NULL, m_init = NULL, 
    MMp = 10^(-2), MMm = 10^(-2), chm = 1.1, chp = 1.1, Jm_scale = 5, 
    Jp_scale = 5, determc = c("both", "drift"), initial.breaks = NULL) 
{
    require(mvtnorm)
    require(coda)
    require(VGAM)
    require(MASS)
    indBrkL = dget("indBrkL.R")
    gibbsdraw_k = dget("gibbsdraw_k.R")
    wetbreak = dget("wetbreak.R")
    lagw = dget("lagw.R")
    inv = dget("inv.R")
    dual2dec = dget("dual2dec.R")
    mlogpost = dget("mlogpost.R")
    Xupdate = function(y, p, k = NULL, type, y0 = NULL, sbreaks = m > 
        0) {
        n = length(y)
        stopifnot(k[length(k)] - 1 == n)
        if (is.null(y0)) {
            y0 = rep(0, p)
        }
        (ylag1 = matrix(c(y0[1], y[1:(n - 1)]), ncol = 1))
        if (p > 1) {
            (d1.y = matrix(c(y0[1], diff(y)), ncol = 1))
            (d1.ylagk = lagw(d1.y, k = p - 1, initial = y0[-p])[, 
                -1])
            (RiteX = cbind(ylag1, d1.ylagk))
            colnames(RiteX) = c("y1", paste("dy", 1:(p - 1), 
                sep = ""))
        }
        else if (p == 1) {
            (RiteX = cbind(ylag1))
            colnames(RiteX) = "y1"
        }
        else {
            stop("AR order 'p' must be at least 1!\n\n")
        }
        if (sbreaks && !is.null(k)) {
            Dt = switch(type[1], drift = indBrkL(k)$ei1, both = indBrkL(k)$ei)
            X = cbind(Dt, RiteX)
        }
        else {
            X = switch(type[1], drift = matrix(cbind(1, RiteX), 
                ncol = ncol(RiteX) + 1, dimnames = list(NULL, 
                  c("drift", colnames(RiteX)))), both = matrix(cbind(1, 
                1:n, RiteX), ncol = ncol(RiteX) + 2, dimnames = list(NULL, 
                c("drift", "trend", colnames(RiteX)))))
            cat("No structural breaks model!\n")
        }
        return(list(X = X, RiteX = RiteX))
    }
    modelID = function(p, p_max, type, m_max, m) {
        if (p > p_max || m > m_max) {
            stop("Wrong number of breaks or number of lags specified!\n")
        }
        if (m > 0) {
            bb = numeric(m_max + 1)
            bb[1:(m + 1)] = rep(1, m + 1)
            dt = switch(type[1], drift = rep(bb, 1), both = rep(bb, 
                2))
        }
        else {
            dt = switch(type[1], drift = 1, both = rep(1, 2))
        }
        st = numeric(p_max)
        st[1:p] = rep(1, p)
        return(c(dt, st))
    }
    moments12 = function(y, X, s2eps) {
        Sigma2 = ginv(as.double(s2eps)^(-1) * t(X) %*% X + as.double(s2eps)^(-1) * 
            diag(ncol(X)))
        mu = (as.double(s2eps)^(-1)) * Sigma2 %*% t(X) %*% y
        rownames(mu) = rownames(Sigma2)
        return(list(Sigma2 = Sigma2, mu = mu))
    }
    start_bpoints = function(m_max, nobs, enter.dates = NULL) {
        save.bpoints = vector("list", m_max)
        for (m in 1:m_max) {
            if (is.null(enter.dates)) {
                k = c(1, (1:m) * floor(nobs/(m + 1)), nobs + 
                  1)
            }
            else {
                if (length(enter.dates) != m_max) {
                  stop("'enter.dates' must be a numeric vector of length equal to 'm_{max}'")
                }
                else {
                  k = c(1, enter.dates[1:m], nobs + 1)
                }
            }
            cat("\nBreak dates:", k, "\n\n")
            save.bpoints[[m]] = rbind(save.bpoints[[m]], k)
        }
        return(save.bpoints = save.bpoints)
    }
    message("Note: Zeros are used as initial values!\n")
    if (any(is.na(y)) == TRUE) {
        y = ts(na.omit(y))
    }
    (n = length(y))
    (type = match.arg(arg = determc, choices = c("drift", "both"), 
        several.ok = TRUE))
    save.bpoints = brfreq = vector("list", m_max)
    save.lags = nofb = mdec = NULL
    if (is.null(p_init)) {
        p = sample(seq(1, p_max), 1)
    }
    else {
        p = p_init
    }
    if (!is.null(p_fix)) {
        p = p_fix
    }
    save.lags = c(save.lags, p)
    if (m_max > 0) {
        if (!is.null(initial.breaks)) {
            stopifnot(length(initial.breaks) == m_max && is.list(initial.breaks))
            save.bpoints = initial.breaks
        }
        else {
            for (m in 1:m_max) {
                k = c(1, (1:m) * floor(n/(m + 1)), n + 1)
                cat("\nBreak dates:", k, "\n\n")
                save.bpoints[[m]] = rbind(save.bpoints[[m]], 
                  k)
            }
        }
        if (is.null(m_init)) {
            m = sample(1:m_max, 1)
            k = as.numeric(save.bpoints[[m]])
        }
        else {
            m = ifelse(m_max > 0, yes = m_init, no = 0)
        }
    }
    else {
        m_max = 0
        k = NULL
        m = 0
    }
    if (!is.null(m_fix)) {
        m = m_max = m_fix
        k = c(1, (1:m) * floor(n/(m + 1)), n + 1)
        save.bpoints = brfreq = NULL
        save.bpoints = rbind(save.bpoints, k)
    }
    nofb = c(m, nofb)
    xy = Xupdate(y = y, p = p, k = k, type = type)
    X = xy$X
    RiteX = xy$RiteX
    npar = switch(type[1], drift = (2 * m_max + p_max + 1), both = (3 * 
        m_max + p_max + 2))
    if (npar >= n) {
        stop("\nThere are more parameters to estimate than observations available!\n")
    }
    accrateA = accrateB = numeric(1)
    s2eps = numeric(MHsim + burn)
    (theta = solve(t(X) %*% X) %*% t(X) %*% y)
    id = which(colnames(X) == "y1")
    pAR = theta[id]
    err = y - X %*% theta
    d = length(theta)
    s2eps[1] = t(err) %*% err/(n - d)
    v0 = 1.001
    lambda0 = 2.0001
    mom.theta = moments12(y = y, X = X, s2eps = s2eps[1])
    mu_theta = mom.theta$mu
    var_theta = mom.theta$Sigma2
    logpost0p = logpost0m = mlogpost(yy = y, Xgam = X, betatil.gam = theta) - 
        d * log(n, base = 10)
    ii = 2
    hh = 1
    repeat {
        cat("\nDraw ", ii, " (", length(save.lags), "/", length(nofb), 
            "/", MHsim + burn, "):\n", sep = "")
        (pst = ceiling(rlaplace(1, location = p, scale = Jp_scale)))
        (pstar = ifelse(pst <= 0 || pst >= p_max, yes = sample(1:p_max, 
            1), no = pst))
        if (runif(1) < chp && is.null(p_fix)) {
            xy.new = Xupdate(y = y, p = pstar, k = k, type = type)
            Xstar = xy.new$X
            RiteX.star = xy$RiteX.star
            mom.theta_star = moments12(y = y, X = Xstar, s2eps = s2eps[ii - 
                1])
            mu_theta_star = mom.theta_star$mu
            var_theta_star = mom.theta_star$Sigma2
            dstar = length(mu_theta_star)
            JA10 = dlaplace(p, location = pstar, scale = Jp_scale)
            JA01 = dlaplace(pstar, location = p, scale = Jp_scale)
            theta.star = try(ginv(t(Xstar) %*% Xstar) %*% t(Xstar) %*% 
                y, silent = T)
            if (p != pstar) {
                logpost1p = mlogpost(yy = y, Xgam = Xstar, betatil.gam = theta.star) - 
                  log(n, base = 10) * dstar
            }
            else {
                logpost1p = logpost0p
            }
            cat("proposal 'p*'(lpost):", logpost1p, "\n")
            cat("current 'p'(lpost):", logpost0p, "\n")
            (log_u = log(runif(1) * MMp))
            accA = 0
            cat("proposal 'p*':", pstar, "\n")
            cat("current 'p':", p, "\n")
            if (logpost1p - logpost0p + log(JA10) - log(JA01) > 
                log_u) {
                X = Xstar
                p = pstar
                var_theta = var_theta_star
                mu_theta = mu_theta_star
                d = length(mu_theta)
                RiteX = RiteX.star
                logpost0p = logpost1p
                save.lags = c(save.lags, p)
                accA = 1
                cat("'p*' accepted!\n")
            }
            else {
                if (runif(1) < 0.1) {
                  logpost0p = -10^10
                }
            }
            mbits = modelID(p = p, p_max = p_max, type = type, 
                m_max = m_max, m = m)
            mdec = c(mdec, dual2dec(mbits))
            accrateA[ii] = (accrateA[ii - 1] * (ii - 1) + accA)/ii
            if (ii > 1) {
                cat("accept.rate 'p':", accrateA[ii], "\n\n")
            }
        }
        mst = ceiling(rlaplace(1, location = m, scale = Jm_scale))
        mstar = ifelse(mst < 0 || mst > m_max, yes = sample(0:m_max, 
            1), no = mst)
        if (runif(1) < chm && m_max > 0 && is.null(m_fix)) {
            if (mstar > 0) {
                kstar = gibbsdraw_k(y = y, RiteX = RiteX, m = mstar, 
                  type = type, s2eps = s2eps[ii - 1], start.breaks = save.bpoints[[mstar]][nrow(save.bpoints[[mstar]]), 
                    ])
                kstar.cor = c(kstar[1], kstar[-c(1, length(kstar))] - 
                  1, kstar[length(kstar)])
            }
            else {
                kstar = kstar.cor = NULL
            }
            print(kstar)
            print(kstar.cor)
            xy.new = Xupdate(y = y, p = p, k = kstar.cor, type = type)
            Xstar = xy.new$X
            RiteX.star = xy.new$RiteX
            theta.star = try(ginv(t(Xstar) %*% Xstar) %*% t(Xstar) %*% 
                y, silent = T)
            mom.theta_star = moments12(y = y, X = Xstar, s2eps = s2eps[ii - 
                1])
            (mu_theta_star = mom.theta_star$mu)
            (var_theta_star = mom.theta_star$Sigma2)
            dstar = length(mu_theta_star)
            JB10 = dlaplace(m, location = mstar, scale = Jm_scale)
            JB01 = dlaplace(mstar, location = m, scale = Jm_scale)
            if (m != mstar) {
                logpost1m = mlogpost(yy = y, Xgam = Xstar, betatil.gam = theta.star) - 
                  log(n, base = 10) * dstar
            }
            else {
                logpost1m = logpost0m
            }
            cat("proposal 'm*'(lpost):", logpost1m, "\n")
            cat("current 'm'(lpost):", logpost0m, "\n")
            (log_v = log(runif(1) * MMm))
            accB = 0
            cat("proposal 'm*':", mstar, "\n")
            cat("current 'm':", m, "\n")
            if (logpost1m - logpost0m + log(JB10) - log(JB01) > 
                log_v) {
                m = mstar
                X = Xstar
                k = kstar
                var_theta = var_theta_star
                mu_theta = mu_theta_star
                d = length(mu_theta)
                RiteX = RiteX.star
                logpost0m = logpost1m
                if (m > 0) {
                  save.bpoints[[m]] = rbind(save.bpoints[[m]], 
                    k)
                }
                nofb = c(m, nofb)
                accB = 1
                cat("'m*' accepted!\n")
            }
            else {
                if (runif(1) < 0.2) {
                  logpost0 = -10^10
                }
            }
            mbits = modelID(p = p, p_max = p_max, type = type, 
                m_max = m_max, m = m)
            mdec = c(mdec, dual2dec(mbits))
            accrateB[ii] = (accrateB[ii - 1] * (ii - 1) + accB)/ii
            if (ii > 1) {
                cat("accept.rate 'm':", accrateB[ii], "\n")
            }
        }
        else {
            if (m > 0) {
                (k = gibbsdraw_k(y = y, RiteX = RiteX, m = m, 
                  type = type, s2eps = s2eps[ii - 1], start.breaks = save.bpoints[nrow(save.bpoints), 
                    ]))
                save.bpoints = rbind(save.bpoints, k)
                xy.new = Xupdate(y = y, p = p, k = k, type = type)
                X = xy.new$X
                RiteX = xy.new$RiteX
                mom.theta = moments12(y = y, X = X, s2eps = s2eps[ii - 
                  1])
                (mu_theta = mom.theta$mu)
                (var_theta = mom.theta$Sigma2)
                d = length(mu_theta)
            }
            hh = hh + 1
        }
        theta = t(rmvnorm(1, mean = mu_theta, sigma = var_theta))
        id = which(colnames(X) == "y1")
        err = y - X %*% theta
        if (!is.null(p_fix) && !is.null(m_fix)) {
            pAR = cbind(pAR, theta)
            rownames(pAR) = colnames(X)
        }
        else {
            pAR = cbind(pAR, theta[id])
        }
        v.e = v0 + n/2
        lambda.e = lambda0 + t(err) %*% err/2
        (s2eps[ii] = 1/rgamma(1, shape = v.e, rate = lambda.e))
        if (length(nofb) < (MHsim + burn) && length(save.lags) < 
            (MHsim + burn) && hh < (MHsim + burn)) {
            ii = ii + 1
        }
        else {
            break
        }
        cat("\nCurrent (p/m): (", p, ",", m, ")\n", sep = "")
    }
    if (m_max > 0) {
        if (is.null(m_fix)) {
            for (dd in 1:m_max) {
                brp = save.bpoints[[dd]]
                brp[, ncol(brp)] = brp[, ncol(brp)] - 1
                brfreq[[dd]] = round(wetbreak(kMat = brp, n = n), 
                  3)
            }
            nofbs = mcmc(nofb[-c(1:burn)])
        }
        else {
            brp = save.bpoints
            brp[, ncol(brp)] = brp[, ncol(brp)] - 1
            brfreq = round(wetbreak(kMat = brp, n = n), 3)
            nofbs = NULL
        }
    }
    else {
        nofbs = brfreq = NULL
    }
    return(list(bp = brfreq, lags = mcmc(save.lags[-c(1:burn)]), 
        nofb = nofbs, accrateA = accrateA, mdec = mdec, accrateB = accrateB, 
        s2eps = mcmc(s2eps[-c(1:burn)]), pAR = mcmc(t(pAR[, -c(1:burn)]))))
}
