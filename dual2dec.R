function (bits) 
{
    lga = length(bits)
    dez = 0
    m = lga
    for (i in 1:lga) {
        dez = dez + bits[i] * 2^(m - 1)
        m = m - 1
    }
    dez
}
