function (data, burn = 300, MHsim = 1000, m_max = 5, p_max = 5, 
    p_fix = NULL, m_fix = NULL, p_init = NULL, m_init = NULL, 
    MMp = 10^(-2), MMm = 10^(-2), chm = 1.1, chp = 1.1, Jm_scale = 5, 
    Jp_scale = 5, v0 = 2.001, lambda0 = 0.001, thinPAR = 2, thinS2 = 1, 
    thinM = 1, thinP = 1, determc = c("both", "drift"), initial.breaks = NULL, 
    Mlife = 0.1, Plife = 0.1, CC = 100, ppL = 10^(-2), ppU = 10^(-10), 
    tune = 50, tria.m = c(A = 0, C = 1), tria.p = c(A = 1, C = 2)) 
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
    mlogpost = dget("mlogpost2.R")
    Xupdate = function(y, p, k = NULL, type, y0 = F, sbreaks = m > 
        0) {
        lagw = dget("lagw.R")
        indBrkL = dget("indBrkL.R")
        n = length(y)
        if (y0) {
            y0 = rep(0, p)
            (ylag1 = matrix(c(y0[1], y[1:(n - 1)]), ncol = 1))
            y01 = y0[1]
            y0p = y0[-p]
        }
        else {
            y0 = y01 = y0p = NULL
            ylag1 = matrix(y[p:(n - 1)], ncol = 1)
        }
        if (p > 1) {
            (d1.y = matrix(c(y01, diff(y)), ncol = 1))
            (d1.ylagk = lagw(d1.y, k = p - 1, initial = y0p)[, 
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
            (Dt = switch(type[1], drift = indBrkL(k)$ei1, both = indBrkL(k)$ei))
            X = cbind(Dt, RiteX)
        }
        else {
            X = switch(type[1], drift = matrix(cbind(1, RiteX), 
                ncol = ncol(RiteX) + 1, dimnames = list(NULL, 
                  c("drift", colnames(RiteX)))), both = matrix(cbind(1, 
                (p + 1):n, RiteX), ncol = ncol(RiteX) + 2, dimnames = list(NULL, 
                c("drift", "trend", colnames(RiteX)))))
            cat("No structural breaks model!\n")
        }
        return(list(ynew = y[-c(1:p)], X = X, RiteX = RiteX))
    }
    modelID = function(p, p_max, type, m_max, m) {
        if (p > p_max || m > m_max) {
            stop("Wrong number of breaks or number of lags specified!\n")
        }
        if (type[1] == "drift") {
            dt = rep(0, m_max + 1)
            dt[1:(m + 1)] = rep(1, m + 1)
        }
        if (type[1] == "both") {
            dt = rep(0, 2 * (m_max + 1))
            dt[1:(2 * (m + 1))] = rep(1, 2 * (m + 1))
        }
        st = numeric(p_max)
        st[1:p] = rep(1, p)
        return(c(dt, st))
    }
    moments12 = function(y, X, s2eps) {
        Sigma2 = inv(as.double(s2eps)^(-1) * t(X) %*% X + as.double(s2eps)^(-1) * 
            diag(ncol(X)))
        mu = (as.double(s2eps)^(-1)) * Sigma2 %*% t(X) %*% y
        rownames(mu) = rownames(Sigma2)
        return(list(Sigma2 = Sigma2, mu = mu))
    }
    initial_pm = function(p_fix, m_fix, m_max, p_max) {
        save.bpoints = NULL
        if (!is.null(p_fix) && !is.null(m_fix) && m_fix > 0) {
            p_max = p_fix
            m_max = m_fix
            save.bpoints = matrix(, ncol = m_max + 2)
            nNew = n - p_fix
            k = matrix(c(1, (1:m_fix) * floor(nNew/(m_fix + 1)), 
                nNew + 1), nrow = 1)
            save.bpoints[1, ] = k
        }
        if (!is.null(p_fix) && is.null(m_fix) && m_max > 0) {
            save.bpoints = vector("list", m_max)
            for (m in 1:m_max) {
                datearr = array(, dim = c(1, m + 2, p_fix))
                nNew = n - p_fix
                k = matrix(c(1, (1:m) * floor(nNew/(m + 1)), 
                  nNew + 1), nrow = 1)
                save.bpoints[[m]] = k
            }
        }
        if (is.null(p_fix) && !is.null(m_fix) && m_fix > 0) {
            save.bpoints = vector("list", p_max)
            for (p in 1:p_max) {
                nNew = n - p
                k = matrix(c(1, (1:m_fix) * floor(nNew/(m_fix + 
                  1)), nNew + 1), nrow = 1)
                save.bpoints[[p]] = k
            }
        }
        if (is.null(p_fix) && is.null(m_fix) && m_max > 0) {
            save.bpoints = vector("list", m_max)
            if (m_max > 0) {
                for (m in 1:m_max) {
                  datearr = array(, dim = c(1, m + 2, p_max))
                  for (p in 1:p_max) {
                    nNew = n - p
                    k = matrix(c(1, (1:m) * floor(nNew/(m + 1)), 
                      nNew + 1), nrow = 1)
                    save.bpoints[[m]][[p]] = k
                  }
                }
            }
        }
        return(save.bpoints)
    }
    correct_datesfreq = function(p_fix, m_fix, m_max, p_max, 
        save.bpoints) {
        brfreq = NULL
        if (!is.null(p_fix) && !is.null(m_fix) && m_fix > 0) {
            (save.bpoints[, -1] = save.bpoints[, -1] + (p_fix - 
                1))
            brfreq = matrix(, ncol = ncol(save.bpoints), nrow = nrow(save.bpoints))
            brfreq = round(wetbreak(kMat = save.bpoints, n = n), 
                6)
        }
        if (!is.null(p_fix) && is.null(m_fix) && m_max > 0) {
            brfreq = vector(class(save.bpoints), length = length(save.bpoints))
            for (m in 1:m_max) {
                (save.bpoints[[m]][, -1] = save.bpoints[[m]][, 
                  -1] + (p_fix - 1))
                brfreq[[m]] = round(wetbreak(kMat = save.bpoints[[m]], 
                  n = n), 6)
            }
        }
        if (is.null(p_fix) && !is.null(m_fix) && m_fix > 0) {
            brfreq = vector(class(save.bpoints), length = length(save.bpoints))
            for (p in 1:p_max) {
                (save.bpoints[[p]][, -1] = save.bpoints[[p]][, 
                  -1] + (p - 1))
                brfreq[[p]] = round(wetbreak(kMat = save.bpoints[[p]], 
                  n = n), 6)
            }
        }
        if (is.null(p_fix) && is.null(m_fix) && m_max > 0) {
            brfreq = vector(class(save.bpoints), length = length(save.bpoints))
            for (m in 1:m_max) {
                for (p in 1:p_max) {
                  (save.bpoints[[m]][[p]][, -1] = save.bpoints[[m]][[p]][, 
                    -1] + (p - 1))
                  brfreq[[m]][[p]] = round(wetbreak(kMat = save.bpoints[[m]][[p]], 
                    n = n), 6)
                }
            }
        }
        return(brfreq)
    }
    inout_bdates = function(p_fix, m_fix, m, p, save.bpoints, 
        k_in = NULL) {
        if (is.null(k_in)) {
            if (!is.null(p_fix) && !is.null(m_fix) && m_fix > 
                0) {
                (k = as.numeric(save.bpoints[nrow(save.bpoints), 
                  ]))
            }
            if (is.null(p_fix) && !is.null(m_fix) && m_fix > 
                0) {
                (k = as.numeric(save.bpoints[[p]][nrow(save.bpoints[[p]]), 
                  ]))
            }
            if (!is.null(p_fix) && is.null(m_fix) && m > 0) {
                (k = as.numeric(save.bpoints[[m]][nrow(save.bpoints[[m]]), 
                  ]))
            }
            if (is.null(p_fix) && is.null(m_fix) && m > 0) {
                (k = as.numeric(save.bpoints[[m]][[p]][nrow(save.bpoints[[m]][[p]]), 
                  ]))
            }
            return(k)
        }
        else {
            if (!is.null(p_fix) && !is.null(m_fix) && m_fix > 
                0) {
                (save.bpoints = rbind(save.bpoints, k_in))
            }
            if (!is.null(p_fix) && is.null(m_fix) && m > 0) {
                (save.bpoints[[m]] = rbind(save.bpoints[[m]], 
                  k_in))
            }
            if (is.null(p_fix) && !is.null(m_fix) && m_fix > 
                0) {
                (save.bpoints[[p]] = rbind(save.bpoints[[p]], 
                  k_in))
            }
            if (is.null(p_fix) && is.null(m_fix) && m > 0) {
                (save.bpoints[[m]][[p]] = rbind(save.bpoints[[m]][[p]], 
                  k_in))
            }
            return(save.bpoints)
        }
    }
    initialise_savebetas = function(p_max, m_max, p_fix, m_fix) {
        if (is.null(m_fix) && is.null(p_fix)) {
            (savebetas = vector(length = p_max * (m_max + 1), 
                "list"))
            (indi_ii = matrix(0, ncol = p_max * (m_max + 1), 
                nrow = 1))
            hh = 1
            nam = NULL
            for (mm in 0:m_max) {
                for (pp in 1:p_max) {
                  (mbits = modelID(p = pp, p_max = p_max, type = type, 
                    m_max = m_max, m = mm))
                  (nam = c(nam, paste(mbits, collapse = "")))
                  savebetas[[hh]] = 0
                  hh = hh + 1
                }
            }
        }
        if (!is.null(m_fix) && is.null(p_fix)) {
            (savebetas = vector(length = p_max, "list"))
            (indi_ii = matrix(0, ncol = p_max, nrow = 1))
            hh = 1
            nam = NULL
            for (pp in 1:p_max) {
                (mbits = modelID(p = pp, p_max = p_max, type = type, 
                  m_max = m_max, m = m_fix))
                (nam = c(nam, paste(mbits, collapse = "")))
                savebetas[[hh]] = 0
                hh = hh + 1
            }
        }
        if (is.null(m_fix) && !is.null(p_fix)) {
            (savebetas = vector(length = m_max + 1, "list"))
            (indi_ii = matrix(0, ncol = m_max, nrow = 1))
            hh = 1
            nam = NULL
            for (mm in 0:m_max) {
                (mbits = modelID(p = p_fix, p_max = p_max, type = type, 
                  m_max = m_max, m = mm))
                (nam = c(nam, paste(mbits, collapse = "")))
                savebetas[[hh]] = 0
                hh = hh + 1
            }
        }
        if (!is.null(m_fix) && !is.null(p_fix)) {
            (savebetas = vector(length = 1, "list"))
            indi_ii = 0
            (mbits = modelID(p = p_fix, p_max = p_max, type = type, 
                m_max = m_max, m = m_fix))
            (nam = paste(mbits, collapse = ""))
            savebetas[[1]] = 0
        }
        names(savebetas) = nam
        return(list(savebetas = savebetas, indi_ii = indi_ii))
    }
    triangular = function(x, A, B, C, g = F, logS = F, tune = 50) {
        ff = numeric(length(x))
        for (ii in 1:length(x)) {
            if (x[ii] < A) {
                ff[ii] = 0
            }
            if (A <= x[ii] && x[ii] <= C) {
                ff[ii] = 2 * (x[ii] - A)/((B - A) * (C - A))
            }
            if (C < x[ii] && x[ii] <= B) {
                ff[ii] = 2 * (B - x[ii])/((B - A) * (B - C))
            }
            if (B < x[ii]) {
                ff[ii] = 0
            }
        }
        if (logS) {
            ff[which(ff == 0)] <- 10^(-tune)
            ff = log(ff)
        }
        if (g) {
            plot(x, ff, type = "b", xlab = "", ylab = "", main = "Triangular prior", 
                font.main = 11)
            abline(v = c(A, C, B), lty = 3, col = "lightblue")
        }
        return(ff)
    }
    if (any(is.na(data)) == TRUE) {
        data = ts(na.omit(data))
    }
    (n = length(data))
    (type = match.arg(arg = determc, choices = c("drift", "both"), 
        several.ok = TRUE))
    save.lags = nofb = mdec = save_theta = NULL
    if (is.null(m_init)) {
        m = sample(1:m_max, 1)
    }
    else {
        m = ifelse(m_max > 0, yes = m_init, no = 0)
    }
    if (is.null(p_init)) {
        p = sample(1:p_max, 1)
    }
    else {
        p = ifelse(p_max > 1, yes = p_init, no = 1)
    }
    if (!is.null(p_fix)) {
        p = p_max = p_fix
    }
    if (!is.null(m_fix)) {
        m = m_max = m_fix
    }
    save.lags = c(save.lags, p)
    nofb = c(m, nofb)
    postdraws_mp = bics = matrix(0, ncol = p_max, nrow = m_max + 
        1, dimnames = list(paste("m", 0:m_max, sep = ""), paste("p", 
        1:p_max, sep = "")))
    postdraws_mp[m + 1, p] = 1
    if (m > 0) {
        (save.bpoints = initial_pm(p_fix = p_fix, m_fix = m_fix, 
            m_max = m_max, p_max = p_max))
        (k = inout_bdates(p_fix = p_fix, m_fix = m_fix, m = m, 
            p = p, save.bpoints = save.bpoints))
    }
    (nNew = n - p)
    xy = Xupdate(y = data, p = p, k = k, type = type)
    X = xy$X
    RiteX = xy$RiteX
    y = xy$y
    npar = switch(type[1], drift = (2 * m_max + p_max + 1), both = (3 * 
        m_max + p_max + 2))
    if (npar >= n) {
        stop("\nThere are more parameters to estimate than observations available!\n")
    }
    accrateA = accrateB = numeric(1)
    s2eps = numeric(MHsim + burn)
    (sb = initialise_savebetas(p_max = p_max, m_max = m_max, 
        p_fix = p_fix, m_fix = m_fix))
    indi_ii = sb$indi
    (savebetas = sb$save)
    (theta = solve(t(X) %*% X) %*% t(X) %*% y)
    id = which(colnames(X) == "y1")
    pAR = theta[id]
    err = y - X %*% theta
    d = length(theta)
    s2eps[1] = t(err) %*% err/(nNew - d)
    mom.theta = moments12(y = y, X = X, s2eps = s2eps[1])
    (mu_theta = mom.theta$mu)
    var_theta = mom.theta$Sigma2
    npar = switch(type[1], drift = (2 * m + p + 1), both = (3 * 
        m + p + 2))
    (logpost0p = logpost0m = mlogpost(yy = y, Xgam = X, a = v0, 
        b = lambda0, cc = CC) - npar * log(n, base = exp(1)))
    ii = jj = hh = 2
    repeat {
        cat("\nDraw ", ii, " (", length(save.lags), "/", length(nofb), 
            "/", MHsim + burn, "):\n", sep = "")
        (pst = ceiling(rlaplace(1, location = p, scale = Jp_scale)))
        pstar = pst
        if (pst < 1) {
            (pstar = 1)
        }
        if (pst > p_max) {
            (pstar = p_max)
        }
        (nNew.star = n - pstar)
        npar.star = switch(type[1], drift = (m + pstar + 1), 
            both = (2 * m + pstar + 2))
        if (runif(1) < chp && is.null(p_fix)) {
            if (m > 0) {
                (kstar = inout_bdates(p_fix = p_fix, m_fix = m_fix, 
                  m = m, p = pstar, save.bpoints = save.bpoints))
            }
            else {
                kstar = NULL
            }
            xy.new = Xupdate(y = data, p = pstar, k = kstar, 
                type = type)
            Xstar = xy.new$X
            RiteX.star = xy.new$RiteX.star
            ystar = xy.new$y
            mom.theta_star = moments12(y = ystar, X = Xstar, 
                s2eps = s2eps[ii - 1])
            mu_theta_star = mom.theta_star$mu
            var_theta_star = mom.theta_star$Sigma2
            JA10 = dlaplace(p, location = pstar, scale = Jp_scale)
            JA01 = dlaplace(pstar, location = p, scale = Jp_scale)
            if (p != pstar) {
                if (is.null(ppL) && is.null(ppU)) {
                  (prob = rep(1/p_max, times = p_max))
                }
                if (!is.null(ppL) && is.null(ppU)) {
                  (prob = c(ppL, rep((1 - (ppL))/(p_max - 1), 
                    p_max - 1)))
                }
                if (is.null(ppL) && !is.null(ppU)) {
                  (prob = c(rep((1 - (ppU))/(p_max - 1), p_max - 
                    1), ppU))
                }
                if (!is.null(ppL) && !is.null(ppU)) {
                  (prob = c(ppL, rep((1 - (ppL + ppU))/(p_max), 
                    p_max), ppU))
                }
                logpost1p = mlogpost(yy = ystar, Xgam = Xstar, 
                  a = v0, b = lambda0, cc = CC)
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
                RiteX = RiteX.star
                y = ystar
                k = kstar
                nNew = nNew.star
                logpost0p = logpost1p
                save.lags = c(save.lags, p)
                accA = 1
                cat("'p*' accepted!\n")
                postdraws_mp[m + 1, pstar] = postdraws_mp[m + 
                  1, pstar] + 1
            }
            else {
                if (runif(1) < 0.05) {
                  logpost0p = -10^10
                }
            }
            accrateA[jj] = (accrateA[jj - 1] * (jj - 1) + accA)/jj
            if (jj > 1) {
                cat("accept.rate 'p':", accrateA[jj], "\n\n")
            }
            jj = jj + 1
        }
        mst = ceiling(rlaplace(1, location = m, scale = Jm_scale))
        mstar = mst
        if (mst < 0) {
            (mstar = 0)
        }
        if (mst > m_max) {
            (mstar = m_max)
        }
        npar.star = switch(type[1], drift = (mstar + p + 1), 
            both = (2 * mstar + p + 2))
        if (runif(1) < chm && m_max > 0 && is.null(m_fix)) {
            if (mstar > 0) {
                (start.breaks = inout_bdates(p_fix = p_fix, m_fix = m_fix, 
                  m = mstar, p = p, save.bpoints = save.bpoints))
                kstar = gibbsdraw_k(y = y, RiteX = RiteX, m = mstar, 
                  type = type, s2eps = s2eps[ii - 1], start.breaks = start.breaks)
                (save.bpoints = inout_bdates(p_fix = p_fix, m_fix = m_fix, 
                  m = mstar, p = p, save.bpoints = save.bpoints, 
                  k_in = kstar))
            }
            else {
                kstar = NULL
            }
            xy.new = Xupdate(y = data, p = p, k = kstar, type = type)
            Xstar = xy.new$X
            RiteX.star = xy.new$RiteX
            mom.theta_star = moments12(y = y, X = Xstar, s2eps = s2eps[ii - 
                1])
            (mu_theta_star = mom.theta_star$mu)
            (var_theta_star = mom.theta_star$Sigma2)
            JB10 = dlaplace(m, location = mstar, scale = Jm_scale)
            JB01 = dlaplace(mstar, location = m, scale = Jm_scale)
            if (m != mstar) {
                if (is.null(ppL) && is.null(ppU)) {
                  (prob = rep(1/(m_max + 1), times = m_max + 
                    1))
                }
                if (!is.null(ppL) && is.null(ppU)) {
                  (prob = c(ppL, rep((1 - (ppL))/(m_max), m_max)))
                }
                if (is.null(ppL) && !is.null(ppU)) {
                  (prob = c(rep((1 - (ppU))/(m_max), m_max), 
                    ppU))
                }
                if (!is.null(ppL) && !is.null(ppU)) {
                  (prob = c(ppL, rep((1 - (ppL + ppU))/(m_max - 
                    1), m_max - 1), ppU))
                }
                logpost1m = mlogpost(yy = y, Xgam = Xstar, a = v0, 
                  b = lambda0, cc = CC) + log(prob)[mstar + 1]
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
                RiteX = RiteX.star
                logpost0m = logpost1m
                nofb = c(m, nofb)
                accB = 1
                postdraws_mp[mstar + 1, p] = postdraws_mp[mstar + 
                  1, p] + 1
                cat("'m*' accepted!\n")
            }
            else {
                if (runif(1) < 0.05) {
                  logpost0m = -10^10
                }
            }
            if (ii > burn) {
                mbits = modelID(p = p, p_max = p_max, type = type, 
                  m_max = m_max, m = m)
                mdec = c(mdec, dual2dec(mbits))
            }
            accrateB[hh] = (accrateB[hh - 1] * (hh - 1) + accB)/hh
            if (hh > 1) {
                cat("accept.rate 'm':", accrateB[hh], "\n")
            }
            hh = hh + 1
        }
        if (m > 0) {
            (start.breaks = inout_bdates(p_fix = p_fix, m_fix = m_fix, 
                m = m, p = p, save.bpoints = save.bpoints))
            k = gibbsdraw_k(y = y, RiteX = RiteX, m = m, type = type, 
                s2eps = s2eps[ii - 1], start.breaks = start.breaks)
            (save.bpoints = inout_bdates(p_fix = p_fix, m_fix = m_fix, 
                m = m, p = p, save.bpoints = save.bpoints, k_in = k))
            xy.new = Xupdate(y = data, p = p, k = k, type = type)
            X = xy.new$X
            RiteX = xy.new$RiteX
            mom.theta = moments12(y = y, X = X, s2eps = s2eps[ii - 
                1])
            (mu_theta = mom.theta$mu)
            (var_theta = mom.theta$Sigma2)
        }
        (theta = t(rmvnorm(1, mean = mu_theta, sigma = var_theta)))
        rownames(theta) = colnames(X)
        id = which(colnames(X) == "y1")
        err = y - X %*% theta
        if (!is.null(p_fix) && !is.null(m_fix)) {
            (save_theta = mcmc(cbind(save_theta, theta)))
        }
        if (ii > burn) {
            if (is.null(m_fix) && is.null(p_fix)) {
                combi = m * p_max + p
            }
            if (!is.null(m_fix) && is.null(p_fix)) {
                combi = p
            }
            if (is.null(m_fix) && !is.null(p_fix)) {
                combi = m + 1
            }
            if (!is.null(m_fix) && !is.null(p_fix)) {
                combi = 1
            }
            zz = indi_ii[combi] + 1
            savebetas[[combi]] = (savebetas[[combi]] * (zz - 
                1) + theta)/zz
            indi_ii[combi] = indi_ii[combi] + 1
        }
        if (!is.null(p_fix) && !is.null(m_fix)) {
            pAR = cbind(pAR, theta)
            rownames(pAR) = colnames(X)
        }
        else {
            pAR = cbind(pAR, theta[id])
        }
        v.e = v0 + nNew/2
        lambda.e = lambda0 + t(err) %*% err/2
        (s2eps[ii] = 1/rgamma(1, shape = v.e, rate = lambda.e))
        if (is.null(m_fix) || is.null(p_fix)) {
            if (length(nofb) < (MHsim + burn) && length(save.lags) < 
                (MHsim + burn)) {
                ii = ii + 1
            }
            else {
                break
            }
        }
        else {
            if (ii < (MHsim + burn)) {
                ii = ii + 1
            }
            else {
                break
            }
        }
        cat("\nCurrent (p/m): (", p, ",", m, ")\n", sep = "")
    }
    brfreq = correct_datesfreq(p_fix = p_fix, m_fix = m_fix, 
        m_max = m_max, p_max = p_max, save.bpoints = save.bpoints)
    if (m_max > 0) {
        if (is.null(m_fix)) {
            nofbs = window(mcmc(nofb[-c(1:burn)]), thin = thinM)
        }
        else {
            nofbs = NULL
        }
    }
    else {
        nofbs = brfreq = NULL
    }
    if (is.null(p_fix)) {
        lags = window(mcmc(save.lags[-c(1:burn)]), thin = thinP)
    }
    else {
        lags = NULL
    }
    bics = (1/postdraws_mp) * bics
    return(invisible(list(bp = brfreq, lags = lags, postdraws_mp = postdraws_mp, 
        bics = bics, thetas = save_theta, nofb = nofbs, accrateA = accrateA, 
        mdec = mdec, accrateB = accrateB, s2eps = window(mcmc(s2eps[-c(1:burn)]), 
            thin = thinS2), indi_ii = indi_ii, savebetas = savebetas, 
        pAR = window(mcmc(t(pAR[, -c(1:burn)])), thin = thinPAR))))
}
