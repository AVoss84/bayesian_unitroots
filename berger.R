function (rho, singular = 10) 
{
    outs = numeric(length(rho))
    ii = 1
    while (ii <= length(rho)) {
        if (abs(rho[ii]) < 1) {
            outs[ii] = 1/(2 * pi * (1 - rho[ii]^2)^(0.5))
        }
        else {
            outs[ii] = 1/(2 * pi * abs(rho[ii]) * (rho[ii]^(2) - 
                1)^(0.5))
        }
        ii = ii + 1
    }
    if (any(outs == Inf)) {
        outs[which(outs == Inf)] = singular
    }
    return(outs)
}
