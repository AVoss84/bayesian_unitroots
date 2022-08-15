function (y, p, k = NULL, type, y0 = NULL, sbreaks = m > 0) 
{
    if (is.null(y0)) {
        y0 = rep(0, p)
    }
    (ylag1 = matrix(c(y0[1], y[1:(n - 1)]), ncol = 1))
    if (p > 1) {
        (d1.y = matrix(c(y0[1], diff(y)), ncol = 1))
        (d1.ylagk = lagw(d1.y, k = p - 1, initial = y0[-p])[, 
            -1])
        (RiteX = cbind(ylag1, d1.ylagk))
        colnames(RiteX) = c("y1", paste("dy", 1:(p - 1), sep = ""))
    }
    else if (p == 1) {
        (RiteX = cbind(ylag1))
        colnames(RiteX) = "y1"
    }
    else {
        stop("AR order 'p' must be at least 1!\n\n")
    }
    if (sbreaks && !is.null(k)) {
        ei = indBrkL(k)$ei
        ei1 = indBrkL(k)$ei1
        X = switch(type[1], drift = cbind(ei1, RiteX), both = cbind(ei, 
            RiteX))
    }
    else {
        X = switch(type[1], drift = cbind(1, RiteX), both = cbind(1, 
            1:nrow(RiteX), RiteX))
        colnames(X)[1:2] = c("drift", "trend")
        cat("No structural breaks model!\n")
    }
    return(list(X = X, RiteX = RiteX))
}
