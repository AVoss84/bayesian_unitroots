function (n = 150, k = c(50, 100), coef = list(ar1 = 0, inter = c(0, 
    0, 0), slope = c(0, 0, 0), errvar = c(1, 1, 1)), start.innov = NULL, 
    innov = NULL, init = NULL, graph = FALSE, ...) 
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
    innov = c(start.innov, innov)
    yt = numeric(length(innov))
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
    phi1 <- coef$ar1
    coefmat[1, ] <- coef$inter
    coefmat[2, ] <- coef$slope
    coefmat[3, ] <- coef$errvar
    regime = matrix(, ncol = (length(k) - 1), nrow = length(time))
    for (i in 2:length(k)) {
        regime[, i - 1] = as.numeric((k[i - 1] <= time) & (time < 
            k[i]))
    }
    for (cc in 2:length(innov)) {
        (X = c(1, cc, innov[cc]))
        if (cc <= length(start.innov)) {
            yt[cc] = t(coefmat)[1, ] %*% X + phi1 * yt[cc - 1]
        }
        else {
            yt[cc] = regime[dd, ] %*% t(coefmat) %*% X + phi1 * 
                yt[cc - 1]
            dd = dd + 1
        }
    }
    if (graph) {
        par(cex.main = 1.2, las = 1, font.main = 11)
        plot(time, yt[(length(yt) - n + 1):length(yt)], type = "l", 
            ylab = "", ...)
        abline(v = kp, col = "gray", lty = 2, lwd = 1)
    }
    return(list(yt = yt[(length(yt) - n + 1):length(yt)], regime = regime, 
        coef = coef))
}
