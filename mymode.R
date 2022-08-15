function (mat, initial = -100) 
{
    stopifnot(is.matrix(mat))
    xy = numeric(2)
    for (ii in 1:nrow(mat)) {
        for (jj in 1:ncol(mat)) {
            if (mat[ii, jj] > initial) {
                xy = c(ii - 1, jj)
                initial = mat[ii, jj]
            }
        }
    }
    return(xy)
}
