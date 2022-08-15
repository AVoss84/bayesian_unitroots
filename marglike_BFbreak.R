function (rhos, data, sigma = NULL, nofb = NULL, TT.orig = NULL, 
    k = NULL, type = c("drift", "both"), norm = TRUE, ...) 
{
    require(Bolstad)
    require(MASS)
    like = loglike = NULL
    cons = 1
    y0 = 0
    qq = 1
    if (is.null(colnames(data))) {
        stop("colnames() missing!\n")
    }
    (data = as.matrix(data))
    (y = data[, 1])
    (ylag = data[, "X.y1"])
    (x = data[, -1])
    (type = match.arg(arg = type, choices = c("drift", "both"), 
        several.ok = TRUE))
    if (!is.null(k) && k > 1) {
        (lagdiff = data[, (ncol(data) - k + 2):ncol(data)])
        (delete.ids = c(which(colnames(data) == "y"), which(colnames(data) == 
            "X.y1"), (ncol(data) - k + 2):ncol(data)))
        (xt = data[, -delete.ids])
        (lagdiff.star = lagdiff - rep(mean(y - ylag), nrow(data)))
        (xrho = switch(type[1], drift = cbind(xt[, 1:(nofb + 
            1)], lagdiff.star), both = cbind(xt, lagdiff.star)))
    }
    else {
        delete.ids = c(which(colnames(data) == "y"), which(colnames(data) == 
            "X.y1"))
        xt = data[, -delete.ids]
        xrho = switch(type[1], drift = cbind(xt[, 1:(nofb + 1)]), 
            both = cbind(xt))
    }
    (m.star = t(xrho) %*% xrho)
    (det.m.star = det(m.star))
    (inv.m.star = ginv(X = m.star))
    for (ii in 1:length(rhos)) {
        (yrho = y - rhos[ii] * ylag)
        (beta.star = inv.m.star %*% (t(xrho) %*% yrho))
        (s.star = t(yrho) %*% yrho + (y0^2)/qq - t(beta.star) %*% 
            m.star %*% beta.star)
        (loglikekern = (-(TT.orig - k + 1)/2) * log(s.star) - 
            0.5 * log(det.m.star) - 0.5 * log(qq))
        loglike = c(loglike, loglikekern)
    }
    like = exp(loglike - max(loglike))
    if (norm) {
        (cons = sintegral(rhos, like)$value)
    }
    return(like = like/cons)
}
