function (yy, Xgam, betatil.gam = NULL, cc = 100, base = exp(1), 
    a, b, kernel.only = F) 
{
    require(MASS)
    if (is.null(betatil.gam)) {
        betatil.gam = matrix(rep(0, ncol(Xgam)), ncol = 1)
    }
    (n = nrow(Xgam))
    Xgam = as.matrix(Xgam)
    betatil.gam = as.matrix(betatil.gam)
    Mgam = diag(1/cc, nrow = ncol(Xgam), ncol = ncol(Xgam))
    (S = diag(n) - Xgam %*% ginv(Mgam + t(Xgam) %*% Xgam) %*% 
        t(Xgam))
    (QF = t(yy - Xgam %*% betatil.gam) %*% S %*% (yy - Xgam %*% 
        betatil.gam))
    if (kernel.only) {
        l_nconst = 0
    }
    else {
        (Det = det(S)^{
            1/2
        })
        (l_nconst = log(Det))
    }
    (logpost_ni = l_nconst - 0.5 * (n + a) * log(b + QF))
    return(as.numeric(logpost_ni))
}
