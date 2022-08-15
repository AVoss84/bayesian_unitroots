function (rhos, data, sigma = NULL, k = NULL, type = c("drift", 
    "both"), norm = TRUE, ...) 
{
    require(Bolstad)
    require(MASS)
    data = as.matrix(data)
    if (ncol(data) > 1) {
        stop("Sorry, 'marglike_NoBreak.R' was written for vector arguments 'data' only not for matrices.\n")
    }
    like = loglike = lagdiff.star = NULL
    y0 = 0
    cons = qq = 1
    TT = nrow(matrix(data, ncol = 1))
    (type = match.arg(arg = type, choices = c("drift", "both"), 
        several.ok = TRUE))
    if (!is.null(k) && k > 1) {
        (d1_y = diff(data))
        (XX = embed(d1_y, dimension = k))
        (xt = embed(cbind(1, 2:TT), dimension = k))
        (lagdiff = as.matrix(XX[, -1]))
        (lagdiff.star = lagdiff - rep(mean(lagdiff[, 1]), nrow(lagdiff)))
    }
    else {
        xt = cbind(1, 2:TT)
    }
    (y = matrix(data[-c(1:k)], ncol = 1))
    ylag = matrix(data[k:(TT - 1)], ncol = 1)
    (xrho = switch(type[1], drift = cbind(xt[, (ncol(xt) - 1)], 
        lagdiff.star), both = cbind(xt[, (ncol(xt) - 1):ncol(xt)], 
        lagdiff.star)))
    (m.star = t(xrho) %*% xrho)
    (inv.m.star = ginv(m.star))
    (det.m.star = det(m.star))
    if ((-0.01 < det.m.star) && (det.m.star < 0.01)) {
        warning("det.m.star==0!! det.m.star<-0.1 has been assigned instead.")
        det.m.star = runif(1, 0, 0.1)
    }
    for (ii in 1:length(rhos)) {
        (yrho = y - rhos[ii] * ylag)
        (beta.star = inv.m.star %*% (t(xrho) %*% yrho))
        (s.star = t(yrho) %*% yrho + (y0^2)/qq - t(beta.star) %*% 
            m.star %*% beta.star)
        (loglikekern = (-(TT - k + 1)/2) * log(s.star) - 0.5 * 
            log(det.m.star) - 0.5 * log(qq))
        loglike = c(loglike, loglikekern)
    }
    (like = exp(loglike - max(loglike)))
    if (norm) {
        (cons = sintegral(rhos, like)$value)
    }
    return(like = like/cons)
}
