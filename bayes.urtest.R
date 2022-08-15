function (data, lags = NULL, breaknr = NULL, determc = "both") 
{
    indBrkL = dget("indBrkL.R")
    gibbsdraw_k = dget("gibbsdraw_k.R")
    wetbreak = dget("wetbreak.R")
    lagw = dget("lagw.R")
    inv = dget("inv.R")
    dual2dec = dget("dual2dec.R")
    rjmcmc = dget("rjbreak_ADF5.R")
    Xupdate = dget("Xupdate.R")
    marglike_BFbreak = dget("marglike_BFbreak.R")
    bayes.test = dget("bayes.test.R")
    if (is.null(lags) || is.null(breaknr)) {
        out = rjmcmc(y = data, burn = 500, MHsim = 2000, m_max = 7, 
            p_max = 10, MMm = 10^(-2), MMp = 10^(-1), p_fix = lags, 
            m_fix = breaknr, chm = 1, chp = 1, Jm_scale = 5, 
            determc = determc)
        (lags = out$lags)
        (sigma2 = out$s2eps)
        (bpoints = out$bp)
        (nofb = out$nofb)
        (k = apply(bpoints[[3]], 2, which.max))
        (k = c(k[-length(k)], k[length(k)] + 1))
        (modelfreq = round(table(nofb)/length(nofb), 4))
        (m = as.numeric(names(which.max(modelfreq))))
        (lagsfreq = round(table(lags)/length(lags), 4))
        (p = as.numeric(names(which.max(lagsfreq))))
    }
    else {
        m = breaknr
        p = lags
    }
    xy = Xupdate(y = data, p = p, k = k, type = determc)
    data_new = data.frame(y = data, X = xy$X)
    grid = round(seq(from = 0.05, to = 1.25, by = 0.001), digit = 5)
    outmarg = marglike_BFbreak(rhos = grid, data = data_new, 
        TT.orig = length(data), norm = F, k = p, type = determc, 
        nofb = m)
    if (is.null(sigma2)) {
        ols = summary(lm(outs$lhs ~ outs$X - 1))
        (sigma = ols$sigma)
    }
    sigma = sqrt(as.numeric(summary(sigma2)$stat["Mean"]))
    bayes = bayes.test(grid, outmarg, res.se = sigma, n = length(data))
    postp = round(bayes$test, 5)[1:3, 1]
    point = round(bayes$test, 5)[4:6, 1]
    cat("\n", paste(rep("#", 30)), "\n", "# Bayesian Unit Root Test with Multiple Structural Breaks #\n", 
        paste(rep("#", 30)), "\n", paste(rep("-", 40)), "\n", 
        paste(rep("-", 40)), "\n", sprintf(" Autoregressive Lag Order (p) = %i\n", 
            p), sprintf(" Used Time Observations (T) = %i\n", 
            n), sprintf(" Number of Breaks (m) = %i\n", m), paste(rep("-", 
            40)), "\n", paste(rep("-", 40)), "\n")
    cat("\n  OUTPUTS:\n\n")
    return(list(bpoints = bpoints, lags = lags, grid = grid, 
        margll = outmarg, post_Jeff = bayes$Jeff.mpost, post_Norm = bayes$Normal.mpost, 
        post_BY = bayes$BY.mpost, post_prob = postp, point_est = point, 
        modusMLL = grid[which.max(outmarg)]))
}
