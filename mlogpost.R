function (yy, Xgam, betatil.gam) 
{
    require(MASS)
    (n = nrow(Xgam))
    Xgam = as.matrix(Xgam)
    betatil.gam = as.matrix(betatil.gam)
    Mgam = diag(100, nrow = ncol(Xgam), ncol = ncol(Xgam))
    (S = diag(n) + Xgam %*% ginv(Mgam) %*% t(Xgam))
    (logpost_ni = -0.5 * log(det(S)) - (n/2) * log(t(yy - Xgam %*% 
        betatil.gam) %*% ginv(S) %*% (yy - Xgam %*% betatil.gam)))
    return(as.double(logpost_ni))
}
