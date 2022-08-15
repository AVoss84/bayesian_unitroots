function (x, k, initial = NULL) 
{
    n = length(x)
    lagk = NULL
    k = as.integer(k)
    if (k < 0) {
        stop("'k' must be a non negative integer value!\n")
    }
    else {
        if (is.null(initial)) {
            for (i in 0:k) {
                lagk = cbind(lagk, x[(k + 1 - i):(n - i)])
            }
        }
        else {
            stopifnot(length(initial) == k)
            i = 1
            lagk = as.matrix(x)
            while (i <= k) {
                lagk = cbind(lagk, c(initial[1:i], x[1:(n - i)]))
                i = i + 1
            }
        }
    }
    colnames(lagk) = paste("lag", 0:k, sep = "")
    return(as.matrix(lagk))
}
