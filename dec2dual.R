function (dec) 
{
    bins = numeric()
    ii = 1
    repeat {
        (integ = dec%/%2)
        (modulo = dec%%2)
        (dec = integ)
        (bins[ii] = modulo)
        if (integ == 0) 
            break
        ii = ii + 1
    }
    rev(bins)
}
