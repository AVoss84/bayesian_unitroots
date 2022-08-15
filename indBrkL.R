function (k) 
{
    if (!is.vector(k)) {
        stop("k must be a numeric vector!\n")
    }
    k[1] = 1
    nofb = length(k) - 2
    tim = 1:(k[length(k)] - 1)
    ei1 <- ei2 <- matrix(0, nrow = k[length(k)] - 1, ncol = length(k) - 
        1)
    colnames(ei1) = paste("Drift", 1:ncol(ei1), sep = "")
    colnames(ei2) = paste("Trend", 1:ncol(ei1), sep = "")
    i = 2
    while (i <= length(k)) {
        ei1[, i - 1] = (k[i - 1] <= tim) & (tim < k[i])
        ei2[, i - 1] = ei1[, i - 1] * tim
        i = i + 1
    }
    return(list(ei = cbind(ei1 = ei1, ei2 = ei2), ei1 = ei1))
}
