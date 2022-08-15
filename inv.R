function (X) 
{
    EV = eigen(X)
    EV$vector %*% diag(1/EV$values) %*% t(EV$vector)
}
