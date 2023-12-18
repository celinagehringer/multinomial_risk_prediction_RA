#########################
# Section 3.3.1 of manuscript
# External validation sample size calculation
# code from: https://github.com/gscollins1973/Logistic-regression-external-validation-sample-size/blob/main/validation_sample_size_22_May_2021.R
#########################

### Sample size for O/E
N_OE <- function(phi, SE_OE = 0.051){
  # phi = prevalence
  # SE_OE = standard error of O/E
  round((1 - phi) / (phi * SE_OE^2), 0)
}

### Sample size for the calibration slope
N_slope <- function(LP, beta0 = 0, beta1 = 1, SE_slope = 0.051){
  # LP = linear predictor (vector)
  # beta = calibration intercept 
  # beta = calibration slope
  # SE_slope = standard error of the calibration slope
  I00 <- mean(exp(beta0 + (beta1 * LP)) / ((1 + exp(beta0 + (beta1 * LP)))^2))
  I01 <- mean(LP * exp(beta0 + (beta1 * LP)) / ((1 + exp(beta0 + (beta1 * LP)))^2))
  I11 <- mean(LP * LP * exp(beta0 + (beta1 * LP)) / ((1 + exp(beta0 + (beta1 * LP)))^2))
  round(I00 / (SE_slope^2 * (I00 * I11 - I01^2)), 0)
}

### Sample size for the c-statistic
N_C <- function(c.index = 0.7, phi = 0.1, target.se = 0.0255, min.opt = 100, max.opt = 10000){
  # needs c.index
  # phi = prevalence
  # target.se = se of the c.index
  se.c <- function(n, cc = c.index, phi2 = phi, target.se2 = target.se){
    zz <- sqrt((cc * (1 - cc)*(1 + (n/2 - 1)*((1 - cc)/(2 - cc)) + ((n/2 - 1) * cc / (1 + cc))))/(n^2 * phi2 * (1 - phi2)))
    abs(zz - target.se2)
  }
  round(optimize(se.c, c(min.opt, max.opt, tol = 0.0001))$minimum, 0)
}

### Sample size for net benefit
N_nb <- function(sens, spec, phi, pt, se.nb = 0.051){
  # sens = sensitivity
  # spec = specificity
  # phi = prevalence
  # pt = threshold
  # se.nb = standard error of net benefit
  w    <- (1 - phi) / phi * (pt / (1 - pt))
  round((sens * (1-sens)/phi + w^2*spec*(1-spec)/(1-phi) + w^2*(1-spec)^2/(phi*(1-phi))) / se.nb^2, 0)
}

########
# Minimum sample size for accurate estimate of Observed / Expected
######## 

#using expected CI width/SE as per paper
N_OE(phi = 0.45, SE_OE = 0.051) # prevalence = 0.44, CI width = 0.2 (i.e., se = 0.051) #470 (LDA outcome)
N_OE(phi = 0.10, SE_OE = 0.102) # prevalence = 0.1, CI width = 0.4 (i.e., se = 0.102) #865 (AEs outcome)

########
# Minimum sample size for accurate estimate of calibration slope
########

N <- 1000000
#for LP2
#using LP and sd from development data
LP  <- data.frame(id = 1:N, value = rnorm(N, 0.02, 0.9))
p   <- exp(LP$value)/(1 + exp(LP$value))
outcome <- rbinom(N, 1, prob=p)
N_slope(LP$value, SE_slope =  0.072) #1518

#for LP3
LP  <- data.frame(id = 1:N, value = rnorm(N, -0.33, 0.9))
p   <- exp(LP$value)/(1 + exp(LP$value))
outcome <- rbinom(N, 1, prob=p)
N_slope(LP$value, SE_slope = 0.072) #1542

########
# Minimum sample size for accurate estimate of concordance statistic
########

N_C(c.index = 0.72, phi = 0.45, target.se = 0.0255) #400
N_C(c.index = 0.53, phi = 0.10, target.se = 0.0300) #1025

### Summary
## sub-model 1
# N_OE = 470
# N_slope = 1518
# N_C = 400
# overall required n = 1518

## sub-model 2
# N_OE = 865
# N_slope = 1542
# N_C = 1025
# overall required n = 1542
