function (p_max = 5, m_max = 5, Mdec, hh, determc, data) 
{
    identify_modelprob = function(p, m, p_max, m_max, Mdec, modelfr, 
        determc) {
        modelID = function(p, p_max, type, m_max, m) {
            if (p > p_max || m > m_max) {
                stop("Wrong number of breaks or number of lags specified!\n")
            }
            if (type[1] == "drift") {
                dt = rep(0, m_max + 1)
                dt[1:(m + 1)] = rep(1, m + 1)
            }
            if (type[1] == "both") {
                dt = rep(0, 2 * (m_max + 1))
                dt[1:(2 * (m + 1))] = rep(1, 2 * (m + 1))
            }
            st = numeric(p_max)
            st[1:p] = rep(1, p)
            return(c(dt, st))
        }
        (model = modelID(p = p, p_max = p_max, type = determc, 
            m_max = m_max, m = m))
        (mdec_model = dual2dec(model))
        mat = cbind(matrix(Mdec, ncol = 1), hh)
        (model_prob = mat[which(mat[, 1] %in% mdec_model), 2])
        return(model_prob)
    }
    pointestimates = function(grid, n, determc, m, p, data) {
        if (m > 0) {
            (dating = results$bpoints[[m]][[p]])
            cump = apply(dating, 2, cumsum)
            (k_map = apply(cump >= 0.5, 2, which.max))
            sigma2 = NULL
            kk_map = k_map
            (kk_map[-1] = kk_map[-1] - (p - 1))
            xy = Xupdate(y = data, p = p, k = kk_map, type = determc, 
                sbreaks = m > 0)
        }
        else {
            xy = Xupdate(y = data, p = p, k = NULL, type = determc, 
                sbreaks = FALSE)
        }
        (X = xy$X)
        (y = xy$ynew)
        (data_new <- data.frame(y = y, X = X))
        if (is.null(sigma2)) {
            ols = summary(lm(y ~ X - 1))
            sigma = ols$sigma
        }
        else {
            sigma = sqrt(length(data))
        }
        if (m > 0) {
            outmarg <- marglike_Break(rhos = grid, data = data_new, 
                TT.orig = n, norm = T, k = p, type = determc, 
                nofb = m)
        }
        else {
            outmarg <- marglike_NoBreak(rhos = grid, data = data, 
                sigma = sigma, k = p, type = determc, norm = T)
        }
        bayes <- bayes.test(grid, outmarg, res.se = sqrt(n), 
            n = n, singular = 3)
        point = round(bayes$test, 5)[5:8, ]
        return(point)
    }
    n = length(data)
    mpr = bma = 0
    pts.est = mfr = NULL
    for (ii in 1:p_max) {
        for (jj in 0:m_max) {
            print(ii)
            print(jj)
            mprob = identify_modelprob(p = ii, m = jj, p_max = p_max, 
                m_max = m_max, Mdec = Mdec, modelfr = hh, determc = determc)
            pest = pointestimates(grid = grid, n = n, determc = determc, 
                m = jj, p = ii, data = data)
            pts.est = rbind(pts.est, pest[1:2])
            mfr = rbind(mfr, mprob)
            mpr = mpr + mprob
            bma = bma + pest[1:2] * mprob
        }
    }
    print(cbind(ptest = pts.est, modelf = mfr))
    return(bma)
}
