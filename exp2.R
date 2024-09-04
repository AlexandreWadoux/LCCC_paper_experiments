
# Second limitation: The values can not be compared across studies since the $\rho_{c}$ 
# values are sensitive to the variance of the measured data

# function to calculate the CCC
eval <- function(x, y){
  
  # concordance correlation coefficient
  n <- length(x)
  sdx <- sd(x, na.rm = T)
  sdy <- sd(y, na.rm = T)
  r <- stats::cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')
  # scale shift
  v <- sdx / sdy
  sx2 <- var(x, na.rm = T) * (n - 1) / n
  sy2 <- var(y, na.rm = T) * (n - 1) / n
  # location shift relative to scale
  u <- (mean(x, na.rm = T) - mean(y, na.rm = T)) / ((sx2 * sy2)^0.25)
  Cb <- ((v + 1 / v + u^2)/2)^-1
  rCb <- r * Cb
  rhoC <- round(rCb, digits = 2)
  
  Cb <- round(Cb, digits = 2)
  r <- round(r, digits = 2)
  
  # return the results
  evalRes <- data.frame(rhoC = rhoC, Cb = Cb, r=r)
  
  return(evalRes)
}

# first dataset
measurements1 <- c(21, 19, 22, 21, 22, 19, 18, 20, 16, 19, 20, 21, 21, 20, 21, 21, 21, 22, 19, 17)
prediction1 <- c(21, 18, 19, 20, 22, 24, 20, 21, 19, 23, 22, 22, 21, 20, 23, 18, 24, 22, 20, 20)
residuals1 <- measurements1 - prediction1
mean(measurements1);mean(prediction1);mean(residuals1)
sd(measurements1);sd(prediction1);sd(residuals1)
eval(measurements1, prediction1)

# second dataset
measurements2 <- c(37, 19, 19, 21, 23, 10, 10, 24, 24, 24, 15, 17, 18, 12, 19, 
                  25, 11, 25, 20, 27)
prediction2 <- c(37, 18, 16, 20, 23, 15, 12, 25, 27, 28, 17, 18, 18, 12, 21, 
                22, 14, 25, 21, 30)
residuals2 <- measurements2 - prediction2
mean(measurements2);mean(prediction2);mean(residuals2)
sd(measurements2);sd(prediction2);sd(residuals2)
eval(measurements2, prediction2)
