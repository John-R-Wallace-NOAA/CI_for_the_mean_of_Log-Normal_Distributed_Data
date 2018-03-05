This function is in the JRWToolBox R package, test code is below the function.

    Log.Norm.CI.Cox.t <- function(x = NULL, MeanLogNorm = NULL, SDLogNorm.Pop = NULL, N = NULL, BootReps = 10000, Prob = 0.95) {
    
    "  # Confidence intervals for the mean of a Log-Normal Distribution using the Cox method modified to use a t distribution is from:  "
    "  #         Journal of Statistics Education Volume 13, Number 1 (2005), www.amstat.org/publications/jse/v13n1/olsson.html  "
    "  #         In the paper a simulation study was performed to compare this method with others.  "
    
    "   # This function, including the parametric bootstrapping for the case when only the mean, SD, and N of the log-normal data are known, is by:  "
    "   #        John Wallace John.Wallace@noaa.gov  "
    
    "   # Note that if only the SD of the mean is available, then the argument 'SDLogNorm.Pop' could be set to SD * N .  "
    "   # Note also, that if the exact value of N is not known, a reasonably close estimate of N will still give a good estimate of the CI.  "
    
     if(!is.null(x)) {   # Apply the modified Cox function to log-normal data
    
        N <- length(x)
        MeanLogNorm <- mean(x)
        Var <- var(log(x))
     
        c(exp(mean(log(x)) + Var/2 - qt(1-(1-Prob)/2, N) * sqrt(Var/N + Var^2/(2*(N-1)))), MeanLogNorm,
            exp(mean(log(x)) + Var/2 + qt(1-(1-Prob)/2, N) *sqrt(Var/N + Var^2/(2*(N-1)))))
    
     } else {    # Use the mean, SD, and N from log-normal data and do a parametric bootstrap 
    
        SDNorm <- sqrt(log((SDLogNorm.Pop/MeanLogNorm)^2+1)) # SD in normal space, since CV = sqrt(exp(sigma^2) - 1) for a log-normal
        MeanNorm <- log(MeanLogNorm) - SDNorm^2/2  # Mean in normal space using correction
    
        x <- rlnorm(N * BootReps, MeanNorm, SDNorm) # Create random log-normals in a matrix of size N X BootReps
        x.mat <- matrix(x, ncol = BootReps)
    
        Boots <- apply(x.mat, 2, function(x) {  # Apply the modified Cox function to each column
    
          Var <- var(log(x))
          tmp <- qt(1-(1-Prob)/2, N) * sqrt(Var/N + Var^2/(2*(N-1)))
          c(exp(mean(log(x)) + Var/2 - tmp),
            exp(mean(log(x)) + Var/2 + tmp))
         })
    
     out <- apply(Boots, 1, mean)   # Take the average by row
     c(out[1], MeanLogNorm, out[2]) # Output answer
    
     }
    }
    
    # Data from the reference:
    
    DATA <- c(914.9, 1568.3, 50.5, 94.1, 199.5, 23.8, 70.5, 213.1,
             44.1, 331.7, 139.3, 115.6, 38.4, 357.1, 725.9, 253.2,
             905.6, 155.4, 138.1, 95.2, 75.2, 275.0, 401.1, 653.8,
             390.8, 483.5, 62.6, 128.5, 81.5, 218.5, 308.2, 41.2,
             60.3, 506.9, 221.8, 112.5, 93.7, 199.3, 210.6, 39.2)
    
    
    Log.Norm.CI.Cox.t(DATA) # Compares to 188.0 and 414.7 in the reference
    Log.Norm.CI.Cox.t(Mean = mean(DATA), SD = sqrt(var(DATA)), N = length(DATA), Boot = 10000)
    
    # -----
    
    Prob <- 0.95
    
    (X.LogNorm.Mean <- mean(DATA))
    (X.LogNorm.SD <-sqrt(var(DATA))) 
    
    out.lower <- exp(log(X.LogNorm.Mean) - qnorm(1-(1-Prob)/2) * sqrt(log((X.LogNorm.SD/X.LogNorm.Mean)^2+1)))
    out.upper <- exp(log(X.LogNorm.Mean) + qnorm(1-(1-Prob)/2) * sqrt(log((X.LogNorm.SD/X.LogNorm.Mean)^2+1)))
    
    c(out.lower, X.LogNorm.Mean, out.upper)
    
    
    (X.LogNorm.Mean <- mean(DATA))
    (X.LogNorm.SD <- sqrt(var(DATA))/length(DATA)) # SD of the mean; gets smaller with increasing sample size
    
    out.lower <- exp(log(X.LogNorm.Mean) - qnorm(1-(1-Prob)/2) * sqrt(log((X.LogNorm.SD/X.LogNorm.Mean)^2+1)))
    out.upper <- exp(log(X.LogNorm.Mean) + qnorm(1-(1-Prob)/2) * sqrt(log((X.LogNorm.SD/X.LogNorm.Mean)^2+1)))
    
    c(out.lower, X.LogNorm.Mean, out.upper)
    
    
    # --------------- Generated data ---------------------
    
    # DATA <- rlnorm(5, 5, 1.5); Boots <- 50000
    # DATA <- rlnorm(100, 5, 1.5); Boots <- 10000
    # DATA <- rlnorm(500, 5, 1.5); Boots <- 10000
    
    # DATA <- rlnorm(1000, 5, 1.5); Boots <- 10000
    # DATA <- rlnorm(1000, 5, 0.5); Boots <- 10000
    DATA <- rlnorm(1000, 5, 1); Boots <- 10000
    # DATA <- rlnorm(1000, 5, 3); Boots <- 10000
    
    # DATA <- rlnorm(10000, 5, 1.5); Boots <- 1000
    
    Log.Norm.CI.Cox.t(DATA) 
    Log.Norm.CI.Cox.t(Mean = mean(DATA), SD = sqrt(var(DATA)), N = length(DATA), Boot= Boots)
    
    # -----
    
    Prob <- 0.95
    
    (X.LogNorm.Mean <- mean(DATA))
    (X.LogNorm.SD <-sqrt(var(DATA))) 
    
    out.lower <- exp(log(X.LogNorm.Mean) - qnorm(1-(1-Prob)/2) * sqrt(log((X.LogNorm.SD/X.LogNorm.Mean)^2+1)))
    out.upper <- exp(log(X.LogNorm.Mean) + qnorm(1-(1-Prob)/2) * sqrt(log((X.LogNorm.SD/X.LogNorm.Mean)^2+1)))
    
    c(out.lower, X.LogNorm.Mean, out.upper)
    
    
    (X.LogNorm.Mean <- mean(DATA))
    (X.LogNorm.SD <- sqrt(var(DATA))/length(DATA)) # SD of the mean; gets smaller with increasing sample size
    
    out.lower <- exp(log(X.LogNorm.Mean) - qnorm(1-(1-Prob)/2) * sqrt(log((X.LogNorm.SD/X.LogNorm.Mean)^2+1)))
    out.upper <- exp(log(X.LogNorm.Mean) + qnorm(1-(1-Prob)/2) * sqrt(log((X.LogNorm.SD/X.LogNorm.Mean)^2+1)))
    
    c(out.lower, X.LogNorm.Mean, out.upper)
    


