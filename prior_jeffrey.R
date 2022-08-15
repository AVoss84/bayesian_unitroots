function (rhos, TT, sigma, norm = F) 
{
    require(Bolstad)
    prior = NULL
    for (rr in rhos) {
        if (rr != 1) {
            alpha0 = (TT - (1 - rr^(2 * TT))/(1 - rr^2))/(1 - 
                rr^2)
        }
        else {
            alpha0 = TT * (TT - 1)/2
        }
        prior = c(prior, alpha0)
    }
    prior = sqrt(prior)/(sigma^3)
    if (norm) {
        cons = sintegral(rhos, prior)$value
        prior = prior/cons
    }
    return(prior)
}
