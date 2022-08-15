function (y, RiteX, m, type, s2eps, start.breaks = NULL) 
{
    indBrkL = dget("indBrkL.R")
    inv = dget("inv.R")
    require(MASS)
    k = start.breaks
    cat("\nStart break dates:", k, "\n")
    i = 2
    while (i <= (m + 1)) {
        kBeg = ifelse(i == 2, yes = k[i - 1] + 2, no = k[i - 
            1] + 1)
        kEnd = k[i + 1] - 1
        prob = numeric(kEnd - kBeg + 1)
        j = kBeg
        while (j <= kEnd) {
            k[i] = j
            ind = indBrkL(k)
            ei = ind$ei
            ei1 = ind$ei1
            X = switch(type[1], drift = cbind(ei1, RiteX), both = cbind(ei, 
                RiteX))
            b = ginv(t(X) %*% X) %*% t(X) %*% y
            err = y[kBeg:kEnd] - X[kBeg:kEnd, ] %*% b
            if (i == 2) {
                prob[j - k[i - 1] - 1] = -length(prob) * log(sqrt(s2eps)) - 
                  sum(err^2)/(2 * s2eps)
            }
            else {
                prob[j - k[i - 1]] = -length(prob) * log(sqrt(s2eps)) - 
                  sum(err^2)/(2 * s2eps)
            }
            j = j + 1
        }
        prob = exp(prob - max(prob))
        prob = prob/sum(prob)
        if (length(prob) > 0) {
            mn = try(rmultinom(n = 1, size = length(prob), prob = prob), 
                silent = T)
            postki = try(which.max(mn), silent = T)
            k[i] = ifelse(i == 2, yes = postki + k[i - 1] + 1, 
                no = postki + k[i - 1])
        }
        i = i + 1
    }
    cat("\nDrawn break dates:", k, "\n\n")
    return(k)
}
