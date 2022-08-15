function (kMat, n) 
{
    nofb = ncol(kMat)
    nofs = nrow(kMat)
    kWet = matrix(0, nrow = n, ncol = nofb)
    colnames(kWet) = paste("Break", 1:nofb, sep = "")
    tIndex = seq(1, n, by = 1)
    i = 1
    while (i <= nofb) {
        k = 1
        while (k <= n) {
            kWet[k, i] = sum(kMat[, i] == tIndex[k])/nofs
            k = k + 1
        }
        i = i + 1
    }
    return(kWet)
}
