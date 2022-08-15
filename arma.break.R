function (n = 150, k = c(50, 100), coef = list(ar = c(0, 0, 0), 
    ma = c(0, 0, 0), inter = c(0, 0, 0), slope = c(0, 0, 0), 
    errvar = c(1, 1, 1)), init = NULL, start.innov = NULL, innov = NULL, 
    graph = FALSE, ...) 
{
    k = kp = sort(k)
    time = 1:n
    stopifnot(n > k[length(k)])
    dd = 1
    if (!is.null(innov)) {
        stopifnot(length(innov) == n)
    }
    else {
        innov = rnorm(n)
    }
    if (is.null(start.innov)) {
        start.innov = rnorm(10)
    }
    innov = c(start.innov, innov)
    yt = ut = numeric(length(innov))
    if (!is.null(init)) {
        yt[length(start.innov) + 1] = init
    }
    if (sum(c(length(coef$inter), length(coef$slope), length(coef$errvar))%%(length(k) + 
        1)) != 0) {
        stop("Dimension of coefficient vectors must be 'k+1' (with 'k' indicating the number of break dates)!")
    }
    coefmat = matrix(, ncol = length(coef$inter), nrow = 3)
    rownames(coefmat) = c("alpha", "beta", "sigma")
    colnames(coefmat) = paste("Regime", 1:(length(k) + 1), sep = "")
    if (k[1] > 1) {
        k = c(1, k)
    }
    if (k[length(k)] < length(time)) {
        k = c(k, length(time) + 1)
    }
    (phis = switch(length(coef$ar), phi = c(coef$ar, 0, 0), phi = c(coef$ar[1:2], 
        0), phi = coef$ar[1:3]))
    (thetas = switch(length(coef$ma), theta = c(coef$ma, 0, 0), 
        theta = c(coef$ma[1:2], 0), theta = coef$ma[1:3]))
    if (is.null(phis)) {
        phis = rep(0, 3)
    }
    if (is.null(thetas)) {
        thetas = rep(0, 3)
    }
    coefmat[1, ] <- coef$inter
    coefmat[2, ] <- coef$slope
    coefmat[3, ] <- coef$errvar
    regime = matrix(, ncol = (length(k) - 1), nrow = length(time))
    for (i in 2:length(k)) {
        regime[, i - 1] = as.numeric((k[i - 1] <= time) & (time < 
            k[i]))
    }
    for (cc in 4:length(innov)) {
        X = c(1, cc)
        if (cc <= length(start.innov)) {
            (Dt = t(coefmat[1:2, ])[1, ] %*% X)
            (ut[cc] = coefmat[3, 1] * innov[cc])
        }
        else {
            Dt = regime[dd, ] %*% t(coefmat[1:2, ]) %*% X
            (ut[cc] = regime[dd, ] %*% coefmat[3, ] * innov[cc])
            dd = dd + 1
        }
        (AR = phis[1] * yt[cc - 1] + phis[2] * yt[cc - 2] + phis[3] * 
            yt[cc - 3])
        (MA = ut[cc] + thetas[1] * ut[cc - 1] + thetas[2] * ut[cc - 
            2] + thetas[3] * ut[cc - 3])
        (yt[cc] = as.double(Dt + AR + MA))
    }
    if (graph) {
        par(cex.main = 1.2, las = 1, font.main = 11)
        plot(time, yt[(length(yt) - n + 1):length(yt)], type = "l", 
            ylab = "")
        abline(v = kp, col = "gray", lty = 2, lwd = 1)
    }
    return(list(yt = yt[(length(yt) - n + 1):length(yt)], regime = regime, 
        coef = coef))
}
