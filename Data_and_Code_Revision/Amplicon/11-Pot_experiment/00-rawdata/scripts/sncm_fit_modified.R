# 修改版 sncm.fit (修复了 confint 导致的无穷大报错问题)
sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
    require(minpack.lm)
    require(Hmisc)
    require(stats4)
    
    options(warn=-1)
    
    #Calculate the number of individuals per community
    N <- mean(apply(spp, 1, sum))
    
    #Calculate the average relative abundance of each taxa across communities
    if(is.null(pool)){
        p.m <- apply(spp, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
    } else {
        p.m <- apply(pool, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
    }
    
    #Calculate the occurrence frequency of each taxa across communities
    spp.bi <- 1*(spp>0)
    freq <- apply(spp.bi, 2, mean)
    freq <- freq[freq != 0]
    
    #Combine
    C <- merge(p, freq, by=0)
    C <- C[order(C[,2]),]
    C <- as.data.frame(C)
    C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
    p <- C.0[,2]
    freq <- C.0[,3]
    names(p) <- C.0[,1]
    names(freq) <- C.0[,1]
    
    #Calculate the limit of detection
    d = 1/N
    
    ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
    # 使用 nlsLM 拟合
    m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
    
    # 【修改核心 1】：用 tryCatch 包裹 confint，如果算不出来置信区间，就优雅地返回 NA，而不是直接让整个程序崩溃闪退
    m.ci <- tryCatch({
        confint(m.fit, 'm', level=0.95)
    }, error = function(e) {
        message("Warning: Cannot compute confidence interval for m. Returning NA.")
        return(c(NA, NA))
    })
    
    ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
    sncm.LL <- function(m, sigma){
        R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
        R = dnorm(R, 0, sigma)
        -sum(log(R))
    }
    
    # 【修改核心 2】：同样用 tryCatch 保护 mle 估计，防止 vmmin 报错
    m.mle <- tryCatch({
        mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
    }, error = function(e) {
        message("Warning: MLE optimization failed. Returning pseudo-object.")
        # 如果失败，返回一个假的结构体让后面的代码能跑下去
        fake_mle <- list(coef = c(m = NA, sigma = NA), details = list(value = NA))
        class(fake_mle) <- "mle_failed"
        return(fake_mle)
    })
    
    ##Calculate Akaike's Information Criterion (AIC)
    # 兼容失败的 MLE
    if (inherits(m.mle, "mle_failed")) {
        aic.fit <- NA
        bic.fit <- NA
        m_mle_coef <- NA
        m_mle_val <- NA
    } else {
        aic.fit <- AIC(m.mle, k=2)
        bic.fit <- BIC(m.mle)
        m_mle_coef <- m.mle@coef['m']
        m_mle_val <- m.mle@details$value
    }
    
    ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
    freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
    Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
    
    pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Calculate AIC for binomial model
    bino.LL <- function(mu, sigma){
        R = freq - pbinom(d, N, p, lower.tail=FALSE)
        R = dnorm(R, mu, sigma)
        -sum(log(R))
    }
    
    # 【修改核心 3】：保护二项分布拟合
    bino.mle <- tryCatch({
        mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
    }, error = function(e) {
        fake_mle <- list(details = list(value = NA))
        class(fake_mle) <- "mle_failed"
        return(fake_mle)
    })
    
    if (inherits(bino.mle, "mle_failed")) {
        aic.bino <- NA
        bic.bino <- NA
        bino_val <- NA
    } else {
        aic.bino <- AIC(bino.mle, k=2)
        bic.bino <- BIC(bino.mle)
        bino_val <- bino.mle@details$value
    }
    
    ##Goodness of fit for binomial model
    bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
    Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
    
    bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Calculate AIC for Poisson model
    pois.LL <- function(mu, sigma){
        R = freq - ppois(d, N*p, lower.tail=FALSE)
        R = dnorm(R, mu, sigma)
        -sum(log(R))
    }
    
    # 【修改核心 4】：保护泊松分布拟合
    pois.mle <- tryCatch({
        mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
    }, error = function(e) {
        fake_mle <- list(details = list(value = NA))
        class(fake_mle) <- "mle_failed"
        return(fake_mle)
    })
    
    if (inherits(pois.mle, "mle_failed")) {
        aic.pois <- NA
        bic.pois <- NA
        pois_val <- NA
    } else {
        aic.pois <- AIC(pois.mle, k=2)
        bic.pois <- BIC(pois.mle)
        pois_val <- pois.mle@details$value
    }
    
    ##Goodness of fit for Poisson model
    pois.pred <- ppois(d, N*p, lower.tail=FALSE)
    Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
    
    pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Results
    if(stats==TRUE){
        fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
        
        # 动态处理可能的 NA 值
        m_ci_val <- ifelse(is.na(m.ci[1]), NA, coef(m.fit) - m.ci[1])
        
        fitstats[1,] <- c(coef(m.fit), m_ci_val, m_mle_coef, m_mle_val, bino_val, pois_val, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
        return(fitstats)
    } else {
        A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
        A <- as.data.frame(A)
        colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
        if(is.null(taxon)){
            B <- A[order(A[,1]),]
        } else {
            B <- merge(A, taxon, by=0, all=TRUE)
            row.names(B) <- B[,1]
            B <- B[,-1]
            B <- B[order(B[,1]),]
        }
        return(B)
    }
}