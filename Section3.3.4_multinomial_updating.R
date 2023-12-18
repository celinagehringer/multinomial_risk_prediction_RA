#########################
# Section 3.4 of manuscript
# Multinomial recalibration (line 6-378) and model refitting (line 378-777)
#########################

#### Model recalibration: updating the slope and intercept 

#results matrix
results_updt <- matrix(nrow = 36, ncol = 1)

#store LPs as matrix for equation below
LP2i <- as.matrix(LP_y2_0m$mean_LP2)
LP3i <- as.matrix(LP_y3_0m$mean_LP3) 

#regress outcome onto LPs
fitnp <- vgam(multi_outcome~s(LP2i, df=2) + s(LP3i, df=2), family=multinomial(refLevel=1), data = data)

fitnp <- vgam(multi_outcome~ LP2i + LP3i, family=multinomial(refLevel=1), data = data)

alpha0_1 <- coefficients(fitnp)[1]
alpha0_2 <- coefficients(fitnp)[2]
alpha1_1 <- coefficients(fitnp)[3]
alpha1_2 <- coefficients(fitnp)[4]
alpha2_1 <- coefficients(fitnp)[5]
alpha2_2 <- coefficients(fitnp)[6]

newLP2 <- alpha0_1 + (alpha1_1*LP2i) + (alpha2_1*LP3i)
newLP3 <- alpha0_2 + (alpha1_2*LP2i) + (alpha2_2*LP3i)

#new coefficients for predictors
#y2 vs y1
results_updt[1,] <- (alpha1_1*(-0.309527966)) + (alpha2_1*0.168876117) #sex
results_updt[2,] <- (alpha1_1*0.009278724) + (alpha2_1*0.002849258) #age
results_updt[3,] <- (alpha1_1*0.418708442) + (alpha2_1*0.207914561) #RF
results_updt[4,] <- (alpha1_1*(-0.580176307)) + (alpha2_1*0.187888068) #HAQ
results_updt[5,] <- (alpha1_1*(-0.371220631)) + (alpha2_1*(-0.227110461)) #DAS28
#y3 vs y1
results_updt[6,] <- (alpha1_2*(-0.309527966)) + (alpha2_2*0.168876117)
results_updt[7,] <- (alpha1_2*0.009278724) + (alpha2_2*0.002849258)
results_updt[8,] <- (alpha1_2*0.418708442) + (alpha2_2*0.207914561)
results_updt[9,] <- (alpha1_2*(-0.580176307)) + (alpha2_2*0.187888068)
results_updt[10,] <- (alpha1_2*(-0.371220631)) + (alpha2_2*(-0.227110461))
#new intercepts
results_updt[11,] <- alpha0_1 + (alpha1_1*1.636419978) + (alpha2_1*(-0.974627064))
results_updt[12,] <- alpha0_2 + (alpha1_2*1.636419978) + (alpha2_2*(-0.974627064))

#mean and sd of updated linear predictor
results_updt[13,] <- mean(newLP2)
results_updt[14,] <- sd(newLP2)
results_updt[15,] <- mean(newLP3)
results_updt[16,] <- sd(newLP3)

#defining predicted probabilities
prob_y1_0m_recal <- (1 / (1 + exp(newLP2) + exp(newLP3)))
prob_y2_0m_recal <- exp(newLP2) / (1 + exp(newLP2) + exp(newLP3))
prob_y3_0m_recal <- exp(newLP3) / (1 + exp(newLP2) + exp(newLP3))

p_data_0m_recal <- as.data.frame(cbind(prob_y1_0m_recal, prob_y2_0m_recal, prob_y3_0m_recal, data$multi_outcome))
colnames(p_data_0m_recal)[1] <- "prob(cat.1)"
colnames(p_data_0m_recal)[2] <- "prob(cat.2)"
colnames(p_data_0m_recal)[3] <- "prob(cat.3)"
colnames(p_data_0m_recal)[4] <- "outcome"

p_0m_recal <- p_data_0m_recal %>% dplyr::select("prob(cat.1)", "prob(cat.2)", "prob(cat.3)") %>% as.matrix(.) #needs to be a matrix to work in the function below

#creating the LP matrix for calibration function below
LP_data_0m_recal <- as.data.frame(cbind(newLP2, newLP3, data$multi_outcome))
colnames(LP_data_0m_recal)[1] <- "LP2vs1"
colnames(LP_data_0m_recal)[2] <- "LP3vs1"
colnames(LP_data_0m_recal)[3] <- "outcome"

LP_0m_recal <- LP_data_0m_recal %>% dplyr::select(LP2vs1, LP3vs1) %>% as.matrix(.) #needs to be a matrix to work in the function below

outcome <- as.numeric(data$multi_outcome) 

Polcal_recal <- function(outcome,k,p,LP,r=1,estimates=FALSE,dfr=2,plotoverall=TRUE,datapoints=TRUE,smoothing=TRUE,smoothpar=1,intercept=FALSE,slope=FALSE,test=FALSE){
  
  # checks
  if(k != length(table(outcome))){stop('The number of categories in outcome does not equal the specified number of categories.')}
  if(dim(p)[2]!=k){stop('The number of columns in p does not equal the specified number of categories.')}
  if(dim(LP)[2]!=k-1){stop('The number of columns in LP does not equal the specified number of categories minus 1.')}
  if(! r %in% 1:k){stop('The given reference category (r) is not possible.')}     
  if(!is.matrix(p)){stop('p is not a matrix.')}
  if(!is.matrix(LP)){stop('LP is not a matrix.')}
  if(isTRUE(plotoverall) && !isTRUE(datapoints) && !isTRUE(smoothing)){stop('For the overall (non-)parametric calibration plots either datapoints or smoothed lines should be requested.')}
  
  # if tests for perfect calibration are requested, automatically calibration intercepts and calibration slopes 
  # are given
  if(isTRUE(test)){intercept<-slope<-TRUE}
  
  # probabilities
  probs <- split(p,col(p))    
  
  # linear predictors necessary for non-parametric calibration plot - give a name to each linear predictor 
  # seperately
  lps <- split(LP,col(LP))
  for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps[[i]]))}
  
  ###############################################
  # parametric logistic recalibration framework 
  # cf. section 2.2.1.                          
  ###############################################
  
  # reference category r
  # LP = matrix with linear predictors
  
  fitp<-vglm(outcome~LP,family=multinomial(refLevel=r))
  if(isTRUE(estimates)){est<-coefficients(fitp)
  names(est) <- paste('EST',names(est),sep='.')}
  
  #######################################################
  # non-parametric recalibration framework (using df=2) 
  # cf. section 2.2.1.                                  
  #######################################################
  
  # reference category r
  # lpi = ith linear predictor
  # for k outcome categories, there are k-1 linear predictors and the code should be adapted to:
  # fitnp <- vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+...+s(lpk-1,df=dfr),family=multinomial(refLevel=r))
  
  fitnp<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr),family=multinomial(refLevel=r))
  
  ###############################################                  
  # Separate (non-)parametric calibration plots
  ###############################################
  
  windows()
  # when running code on Mac: use quartz() instead of windows()
  #quartz()
  
  par(mfrow=c(ceiling(k/2),2))
  
  # non-parametric calibration plot 
  # cf. section 2.2.2.              
  ###################################
  
  par(mfrow=c(ceiling(k/2),2))
  for(i in 1:k){p <- unlist(probs[[i]])
  if(isTRUE(smoothing)){color<-'grey'}else{color<-1+i}
  matplot(p,fitted(fitnp)[,i],type="p",pch=i,col=color,lwd=1,ylab="",xlab="",xlim=0:1,ylim=0:1)
  par(new=T)
  ref <- rbind(c(0,0),c(1,1))
  matplot(ref,ref,type="l",col=1,lwd=2,ylab="Observed proportions",xlab="Predicted probabilities",xlim=0:1,ylim=0:1)
  # smoother for calibration plots 
  ##################################
  # a = smoothing parameter
  if(isTRUE(smoothing)){
    a = smoothpar
    points(smooth.spline(p, fitted(fitnp)[,i],spar=a), type="l", col=(1+i), lwd = 4)}
  # legend
  legende <- c(paste("cat ", i, sep = ""))
  legend(x=0.6, y=(0.2),col=(1+i),lty =1,legend=legende)
  title(main = "Recalibration: Calibration plot")
  par(new=F)}
  
  #############################################            
  # Overall (non-)parametric calibration plot 
  #############################################
  
  if(isTRUE(plotoverall)){
    
    # non-parametric calibration plot 
    # cf. section 2.2.2.              
    ###################################
    
    windows()
    # when running code on Mac: use quartz() instead of windows()
    #quartz() 
    
    if(isTRUE(datapoints)){for(i in 1:k){p <- unlist(probs[[i]])
    matplot(p,fitted(fitnp)[,i],type="p",pch=i,col=(1+i),lwd=1,ylab="",xlab="",xlim=0:1,ylim=0:1)
    par(new=T)}}
    ref <- rbind(c(0,0),c(1,1))
    matplot(ref,ref,type="l",col=1,lwd=2,ylab="Observed proportions",xlab="Predicted  probabilities",xlim=0:1,ylim=0:1)
    # smoother for calibration plots 
    ##################################
    # a = smoothing parameter
    if(isTRUE(smoothing)){a = smoothpar
    for(i in 1:k){p <- unlist(probs[[i]])
    points(smooth.spline(p, fitted(fitnp)[,i],spar=a), type="l", col=(1+i), lwd = 4)}}
    # legend
    for(i in 1:k){if(i <= 3){legende <- c("no LDA","LDA", "AEs")}
      if(i > 3){legende <- c(legende,paste("cat ", i, sep = ""))}}
    legend(x=0.7, y=(0.20+(k-3)*0.05),col=2:(k+1),lty =1,legend=legende)
    title(main = "Recalibration: Calibration plot")
    par(new=F)}
  
  ########################################
  # estimation of calibration intercepts 
  # cf. section 2.2.3. and 2.2.4.        
  ########################################
  
  if(isTRUE(intercept)){int<-vgam(outcome~1,offset=LP,family=multinomial(refLevel=r))
  coeffint<-coefficients(int)
  se<-sqrt(diag(vcov(int)))
  ci1i <- cbind(LL1 = coeffint[1] - qnorm(0.975) * se[1], UL1 = coeffint[1] + qnorm(0.975) * se[1])
  ci2i <- cbind(LL2 = coeffint[2] - qnorm(0.975) * se[2], UL2 = coeffint[2] + qnorm(0.975) * se[2])
  estint <- c(coeffint[1],ci1i,coeffint[2],ci2i)
  names(estint) <- paste('CALINT',c('int1','LLint1','ULint1','int2','LLint2','ULint2'),sep='.')}
  
  ####################################
  # estimation of calibration slopes 
  # cf. section 2.2.3. and 2.2.4.    
  ####################################
  
  # we used constraints to fix some coefficients to zero as appropriate
  # for k outcome categories this code should be changed to:
  # i <- diag(k-1)
  # i2 <- cbind(c(1,rep(0,k-2)))
  # i3 <- cbind(c(0,1,rep(0,k-1)))
  # i4 <- cbind(c(0,0,1,rep(0,k-2)))
  # ... (ij <- cbind(c(rep(0,j-2),1,rep(0,k-j)))
  # ik <- cbind(c(rep(0,k-2),1))
  # clist<-list("(Intercept)"=i,"lp1"=i2,"lp2"=i3,...,"lpk-1"=ik)
  # slopes<-vgam(outcome~lp1+lp2+...+lpk-1,family=multinomial(refLevel=r),constraints=clist)
  
  if(isTRUE(slope)){i<-diag(k-1)
  i2<-rbind(1,0)
  i3<-rbind(0,1)
  clist<-list("(Intercept)"=i,"lp1"=i2,"lp2"=i3)
  slopes<-vgam(outcome~lp1+lp2,family=multinomial(refLevel=r),constraints=clist)
  coeffslopes<-coefficients(slopes)[k:length(coefficients(slopes))]
  se<-sqrt(diag(vcov(slopes)))
  ci1s <- cbind(LL1 = coeffslopes[1] - qnorm(0.975) * se[3], UL1 = coeffslopes[1] + qnorm(0.975) * se[3])
  ci2s <- cbind(LL2 = coeffslopes[2] - qnorm(0.975) * se[4], UL2 = coeffslopes[2] + qnorm(0.975) * se[4])
  estslopes <- c(coeffslopes[1],ci1s,coeffslopes[2],ci2s)
  names(estslopes) <- paste('CALSLOPES',c('lp1','LLlp1','ULlp1','lp2','LLlp2','ULlp2'),sep='.')}
  
  #################################
  # calibration testing          
  # cf. section 2.2.3. and 2.2.4. 
  #################################
  
  # this code requires the bayesm library developed by Peter Rossi
  
  if(isTRUE(test)){
    
    # -2 log-likelihood of model without adaptations
    # for k outcome categories this code should be changed to:
    # alphas <- rep(0,k-1) #(i.e. all intercepts zero)
    # beta1 <- c(1,rep(0,k-2)) #(i.e. first linear predictor for first equation)
    # beta2 <- c(0,1,rep(0,k-3)) #(i.e. second linear predictor for second equation)      
    # betaj <- c(rep(0,j-1),1,rep(0,k-1-j)) #(i.e. jth linear predictor for jth equation)
    # betak <- c(rep(0,k-2),1) #(i.e. kth linear predictor for kth equation)
    # parametersk <- c(alphas, beta1, beta2, ..., betak)
    
    parametersk <- c(0,0,1,0,0,1) #c(alpha1,alpha2,b22,b23,b32,b33)
    Xdk=LP
    x <- createX(p=k,na=0,nd=k-1,Xa=NULL,Xd=Xdk,INT=TRUE,DIFF=FALSE,base=1)
    deviancewithout <- -2*llmnl(parametersk,outcome,x)
    names(deviancewithout)<-c('original deviance')
    
    devint <- deviance(int)
    names(devint)<-c('intercept deviance')
    devslopes <- deviance(slopes)
    names(devslopes)<-c('slopes deviance')
    
    # overall calibration (i.e. calibration intercepts and slopes) 
    ################################################################
    
    poverall<- pchisq(deviancewithout - devslopes, df = 2*(k-1), lower.tail = FALSE)
    
    # calibration intercepts 
    ##########################
    
    pint<- pchisq(deviancewithout - devint, df = k-1, lower.tail = FALSE)
    
    # calibration slopes 
    ######################
    
    pslopes<- pchisq(devint - devslopes, df = k-1, lower.tail = FALSE)
    names(poverall)<-c('p overall')
    names(pint)<-c('p int')
    names(pslopes)<-c('p slopes')}
  
  # Printing of results
  # The probabilities of calibration intercepts and slopes are only shown when the hypothesis of perfect 
  # calibration is rejected.
  
  results<-list(if(isTRUE(estimates)){est}else{'Not requested'},if(isTRUE(intercept)){estint}else{'Not requested'},if(isTRUE(slope)){estslopes}else{'Not requested'},if(isTRUE(test)){c(deviancewithout,devint,devslopes)}else{'Not requested'},if(isTRUE(test)){c(poverall,if(poverall<0.05){c(pint,pslopes)})}else{'Not requested'})
  names(results)<-c("Coefficients of parametric recalibration framework","Calibration Intercepts with 95% CI","Calibration Slopes with 95% CI","Deviances","P-values")
  n <- 1:5
  selection <- c(isTRUE(estimates),isTRUE(intercept),isTRUE(slope),isTRUE(test),isTRUE(test))
  results[n[selection]]
}

Polcal_results_0m_recal <- Polcal_recal(outcome=outcome, k=3, p=p_0m_recal, LP=LP_0m_recal, r=1, estimates=FALSE, dfr=2, plotoverall=TRUE, datapoints=TRUE, smoothing=TRUE, smoothpar=1, intercept=TRUE, slope=TRUE, test=FALSE)
Polcal_results_0m_recal

#Please save or screenshot these calibration plots. There should be 2 non-parametric plots. 

#print results into dataframe
#calibration slope
results_updt[17,] <- Polcal_results_0m_recal$`Calibration Slopes with 95% CI`[[1]] #c-slope for y=2 vs 1
results_updt[18,] <- Polcal_results_0m_recal$`Calibration Slopes with 95% CI`[[2]] #lower CI for y=2 vs 1
results_updt[19,] <- Polcal_results_0m_recal$`Calibration Slopes with 95% CI`[[3]] #upper CI for y=2 vs 1
results_updt[20,] <- Polcal_results_0m_recal$`Calibration Slopes with 95% CI`[[4]] #c-slope for for y=3 vs 1
results_updt[21,] <- Polcal_results_0m_recal$`Calibration Slopes with 95% CI`[[5]] #lower CI for y=3 vs 1
results_updt[22,] <- Polcal_results_0m_recal$`Calibration Slopes with 95% CI`[[6]] #upper CI for y=3 vs 1
#same for calibration intercepts
results_updt[23,] <- Polcal_results_0m_recal$`Calibration Intercepts with 95% CI`[[1]]
results_updt[24,] <- Polcal_results_0m_recal$`Calibration Intercepts with 95% CI`[[2]]
results_updt[25,] <- Polcal_results_0m_recal$`Calibration Intercepts with 95% CI`[[3]]
results_updt[26,] <- Polcal_results_0m_recal$`Calibration Intercepts with 95% CI`[[4]]
results_updt[27,] <- Polcal_results_0m_recal$`Calibration Intercepts with 95% CI`[[5]]
results_updt[28,] <- Polcal_results_0m_recal$`Calibration Intercepts with 95% CI`[[6]]

#### Pairwise c-statistics
#outcome LDA (2) vs no LDA (1)
#filter out outcome 3
cstat_y2_0m_recal <- p_data_0m_recal %>% filter(outcome != 3) %>%
  rename(prob1 = "prob(cat.1)",
         prob2 = "prob(cat.2)",
         prob3 = "prob(cat.3)")
cstat_y2_0m_recal$outcome <- as.factor(cstat_y2_0m_recal$outcome)
#drop the unused outcome category
#cstat_y2_0m$outcome <- droplevels(cstat_y2_0m)
cstat_mod_y2_0m_recal <- roc(outcome~prob2,data=cstat_y2_0m_recal)
results_updt[29,] <- cstat_mod_y2_0m_recal$auc #AUC

se_recal <- sqrt(var(cstat_mod_y2_0m_recal))
cstat_y2_0m_LCI_recal <- (cstat_mod_y2_0m_recal$auc) - 1.96*se_recal 
results_updt[30,] <- cstat_y2_0m_LCI_recal #lower CI
cstat_y2_0m_UCI_recal <- (cstat_mod_y2_0m_recal$auc) + 1.96*se_recal
results_updt[31,] <- cstat_y2_0m_UCI_recal #upper CI

#outcome AEs (3) vs no LDA (1)
#filter out outcome 2
cstat_y3_0m_recal <- p_data_0m_recal %>% filter(outcome != 2) %>%
  rename(prob1 = "prob(cat.1)",
         prob2 = "prob(cat.2)",
         prob3 = "prob(cat.3)")
cstat_y3_0m_recal$outcome <- as.factor(cstat_y3_0m_recal$outcome)
#drop the unused outcome category
cstat_y3_0m_recal$outcome <- droplevels(cstat_y3_0m_recal$outcome)
cstat_mod_y3_0m_recal <- roc(outcome~prob3,data=cstat_y3_0m_recal)
results_updt[32,] <- cstat_mod_y3_0m_recal$auc #AUC

se_recal <- sqrt(var(cstat_mod_y3_0m_recal))
cstat_y3_0m_LCI_recal <- (cstat_mod_y3_0m_recal$auc) - 1.96*se_recal 
results_updt[33,] <- cstat_y3_0m_LCI_recal #lower CI
cstat_y3_0m_UCI_recal <- (cstat_mod_y3_0m_recal$auc) + 1.96*se_recal
results_updt[34,] <- cstat_y3_0m_UCI_recal #upper CI

#PDI
outcome <- as.data.frame(outcome)
pdi_recal <- pdi3cat0(outcome, prob_y1_0m_recal)
pdi_recal 

pdi3_recal <- pdi3v0(outcome, p_0m_recal)
results_updt[35,] <- pdi3_recal$pdi

#Variance explained - R2
#regress LPs onto outcome
citl_model_0m_recal <- multinom(outcome ~ (newLP2 + newLP3), data=LP_data_0m_recal) 

LR_0m_recal <- -2 * (as.numeric(logLik(multinom(outcome ~ 1, data=LP_data_0m_recal))) - 
                       as.numeric(logLik(citl_model_0m_recal)))

R2_coxsnell_0m_recal <- 1 - exp(-LR_0m_recal/length(LP_data_0m_recal$outcome))
R2_coxsnell_0m_recal 

#Nagelkerke calculations - adjusting the cox snell r2.
Nagel_0m_recal <- prev_y1_0m * log(prev_y1_0m/n_0m) + (n_0m-prev_y1_0m) * (log(1-(prev_y1_0m/n_0m)))
R2_0m_recal <- 1-exp((2*Nagel_0m_recal)/n_0m)
final_R2_recal <- R2_coxsnell_0m_recal/R2_0m_recal
results_updt[36,] <- final_R2_recal

#label results dataframe
rownames(results_updt) <- c("beta_sex_y2", "beta_age_y2", "beta_RF_y2", "beta_HAQ_y2", "beta_DAS28_y2", "beta_sex_y3", "beta_age_y3", "beta_RF_y3", "beta_HAQ_y3", 
                            "beta_DAS28_y3", "new_int_y2", "new_int_y3","mean_LP2_0m", "sd_LP2_0m", "mean_LP3_0m", "sd_LP3_0m", "cslope_y2_0m", "cslope_LCI_y2_0m", 
                            "cslope_UCI_y2_0m", "cslope_y3_0m", "cslope_LCI_y3_0m", "cslope_UCI_y3_0m", "cint_y2_0m", "cint_LCI_y2_0m", "cint_UCI_y2_0m", "cint_y3_0m",
                            "cint_LCI_y3_0m", "cint_UCI_y3_0m", "cstat_y2_0m", "cstat_LCI_y2_0m", "cstat_UCI_y2_0m", "cstat_y3_0m", "cstat_LCI_y3_0m", "cstat_UCI_y3_0m",
                            "PDI", "R2_0m")
#limit to 5 decimal places
results_updt <- as.data.frame(results_updt) %>% mutate(across(where(is.numeric), ~ round(., 5)))

#---------------------------------------------------------------------------------------------------

### Model refitting 

#results matrix
results_refit <- matrix(nrow = 25, ncol = 1)

# fitting the multinomial model to the imputed dataset
multi_mo <- with(multi_imp_0m, 
                 multinom(multi_outcome ~ sex + age + HAQ + RF + DAS28_CRP_0m, model = T))

### summarising the model by pooling estimates across the 10 imputed datasets
pooled_sum <- summary(pool(multi_mo)) %>% rename(outcome = y.level)
#OR's and CI's
pooled_sum$odds.ratio <- exp(pooled_sum$estimate)
pooled_sum$confint_lower <- pooled_sum$odds.ratio - 1.96 * pooled_sum$std.error
pooled_sum$confint_upper <- pooled_sum$odds.ratio + 1.96 * pooled_sum$std.error

#check AIC of model
results_refit[1,] <- mean(sapply(multi_mo$analyses, AIC))

### Calculating probabilities and Linear predictors (LP)
#filtering out estimates and covariates for outcome 1 (vs 2)
ests2 <- pooled_sum %>% data.frame("term" = pooled_sum$term,  "estimate" = pooled_sum$estimate) %>% filter(outcome==2)
form <- as.formula(multi_outcome ~ sex + age + HAQ + RF + DAS28_CRP_0m)
Design_matrix <- model.matrix(form, data= stacked_data_0m)

LP2 <- as.numeric(cbind(1,
                        Design_matrix[,ests2$term[-1]]) %*% ests2$estimate)

#filtering out estimates and covariates for outcome 3 (vs 1)
ests3 <- pooled_sum %>% data.frame("term" = pooled_sum$term,  "estimate" = pooled_sum$estimate) %>% filter(outcome==3)

LP3 <- as.numeric(cbind(1,
                        Design_matrix[,ests3$term[-1]]) %*% ests3$estimate)

#### Calculating probabilities for each outcome using the LP's defined above
#probabilities for outcome 1
OutSample_Predictions1 <- stacked_data_0m %>%
  filter(.imp != 0) %>%
  mutate("Pi1" = 1 / (1 + exp(LP2) + exp(LP3)))

#probabilities for y=2
OutSample_Predictions2 <- stacked_data_0m %>%
  filter(.imp != 0) %>%
  mutate("LP2" = LP2,
         "Pi2" = exp(LP2) / (1+ exp(LP2) + exp(LP3)))

#probabilities for y=3
OutSample_Predictions3 <- stacked_data_0m %>%
  filter(.imp != 0) %>%
  mutate("LP3" = LP3,
         "Pi3" = exp(LP3) / (1+ exp(LP2) + exp(LP3)))

### summarising probabilities and LPs across the 10 imputed datasets
# for outcome 1
prob.cat.1 <- OutSample_Predictions1 %>% group_by(.id) %>% summarise(pr_test_mean = mean(Pi1), unique(multi_outcome))
# for outcome 2
prob.cat.2 <- OutSample_Predictions2 %>% group_by(.id) %>% summarise(pr_test_mean = mean(Pi2), unique(multi_outcome)) 
LP2vs1 <- OutSample_Predictions2 %>% group_by(.id) %>% summarise(LP2vs1 = mean(LP2), unique(multi_outcome)) 
# for outcome 3
prob.cat.3 <- OutSample_Predictions3 %>% group_by(.id) %>% summarise(pr_test_mean = mean(Pi3), unique(multi_outcome)) 
LP3vs1 <- OutSample_Predictions3 %>% group_by(.id) %>% summarise(LP3vs1 = mean(LP3), unique(multi_outcome)) 

p_data_0m_refit <- as.data.frame(cbind(prob.cat.1$pr_test_mean, prob.cat.2$pr_test_mean, prob.cat.3$pr_test_mean))
colnames(p_data_0m_refit)[1] <- "prob(cat.1)"
colnames(p_data_0m_refit)[2] <- "prob(cat.2)"
colnames(p_data_0m_refit)[3] <- "prob(cat.3)"

p_0m_refit <- p_data_0m_refit %>% dplyr::select("prob(cat.1)", "prob(cat.2)", "prob(cat.3)") %>% as.matrix(.) #needs to be a matrix to work in the function below

#creating the LP matrix for calibration function below
LP_data_0m_refit <- as.data.frame(cbind(LP2vs1$LP2vs1, LP3vs1$LP3vs1))
colnames(LP_data_0m_refit)[1] <- "LP2vs1"
colnames(LP_data_0m_refit)[2] <- "LP3vs1"

#mean and sd of linear predictor
results_refit[2,] <- mean(LP_data_0m_refit$LP2vs1)
results_refit[3,] <- sd(LP_data_0m_refit$LP2vs1)
results_refit[4,] <- mean(LP_data_0m_refit$LP3vs1)
results_refit[5,] <- sd(LP_data_0m_refit$LP3vs1)

LP_0m_refit <- LP_data_0m_refit %>% dplyr::select(LP2vs1, LP3vs1) %>% as.matrix(.) #needs to be a matrix to work in the function below

outcome <- as.numeric(data$multi_outcome) 

Polcal_refit <- function(outcome,k,p,LP,r=1,estimates=FALSE,dfr=2,plotoverall=TRUE,datapoints=TRUE,smoothing=TRUE,smoothpar=1,intercept=FALSE,slope=FALSE,test=FALSE){
  
  # checks
  if(k != length(table(outcome))){stop('The number of categories in outcome does not equal the specified number of categories.')}
  if(dim(p)[2]!=k){stop('The number of columns in p does not equal the specified number of categories.')}
  if(dim(LP)[2]!=k-1){stop('The number of columns in LP does not equal the specified number of categories minus 1.')}
  if(! r %in% 1:k){stop('The given reference category (r) is not possible.')}     
  if(!is.matrix(p)){stop('p is not a matrix.')}
  if(!is.matrix(LP)){stop('LP is not a matrix.')}
  if(isTRUE(plotoverall) && !isTRUE(datapoints) && !isTRUE(smoothing)){stop('For the overall (non-)parametric calibration plots either datapoints or smoothed lines should be requested.')}
  
  # if tests for perfect calibration are requested, automatically calibration intercepts and calibration slopes 
  # are given
  if(isTRUE(test)){intercept<-slope<-TRUE}
  
  # probabilities
  probs <- split(p,col(p))    
  
  # linear predictors necessary for non-parametric calibration plot - give a name to each linear predictor 
  # seperately
  lps <- split(LP,col(LP))
  for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps[[i]]))}
  
  
  ###############################################
  # parametric logistic recalibration framework 
  # cf. section 2.2.1.                          
  ###############################################
  
  # reference category r
  # LP = matrix with linear predictors
  
  fitp<-vglm(outcome~LP,family=multinomial(refLevel=r))
  if(isTRUE(estimates)){est<-coefficients(fitp)
  names(est) <- paste('EST',names(est),sep='.')}
  
  
  #######################################################
  # non-parametric recalibration framework (using df=2) 
  # cf. section 2.2.1.                                  
  #######################################################
  
  # reference category r
  # lpi = ith linear predictor
  # for k outcome categories, there are k-1 linear predictors and the code should be adapted to:
  # fitnp <- vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr)+...+s(lpk-1,df=dfr),family=multinomial(refLevel=r))
  
  fitnp<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr),family=multinomial(refLevel=r))
  
  
  ###############################################                  
  # Separate (non-)parametric calibration plots
  ###############################################
  
  #windows()
  # when running code on Mac: use quartz() instead of windows()
  quartz()
  
  par(mfrow=c(ceiling(k/2),2))
  
  # non-parametric calibration plot 
  # cf. section 2.2.2.              
  ###################################
  #windows()
  # when running code on Mac: use quartz() instead of windows()
  
  par(mfrow=c(ceiling(k/2),2))
  for(i in 1:k){p <- unlist(probs[[i]])
  if(isTRUE(smoothing)){color<-'grey'}else{color<-1+i}
  matplot(p,fitted(fitnp)[,i],type="p",pch=i,col=color,lwd=1,ylab="",xlab="",xlim=0:1,ylim=0:1)
  par(new=T)
  ref <- rbind(c(0,0),c(1,1))
  matplot(ref,ref,type="l",col=1,lwd=2,ylab="Observed proportions",xlab="Predicted probabilities",xlim=0:1,ylim=0:1)
  # smoother for calibration plots 
  ##################################
  # a = smoothing parameter
  if(isTRUE(smoothing)){
    a = smoothpar
    points(smooth.spline(p, fitted(fitnp)[,i],spar=a), type="l", col=(1+i), lwd = 4)}
  # legend
  legende <- c(paste("cat ", i, sep = ""))
  legend(x=0.6, y=(0.2),col=(1+i),lty =1,legend=legende)
  title(main = "Refitting: Non-parametric calibration plot")
  par(new=F)}
  
  
  #############################################            
  # Overall (non-)parametric calibration plot 
  #############################################
  
  if(isTRUE(plotoverall)){
    #windows()
    # when running code on Mac: use quartz() instead of windows()
    
    # non-parametric calibration plot 
    # cf. section 2.2.2.              
    ###################################
    
    #windows()
    # when running code on Mac: use quartz() instead of windows()
    quartz() 
    
    if(isTRUE(datapoints)){for(i in 1:k){p <- unlist(probs[[i]])
    matplot(p,fitted(fitnp)[,i],type="p",pch=i,col=(1+i),lwd=1,ylab="",xlab="",xlim=0:1,ylim=0:1)
    par(new=T)}}
    ref <- rbind(c(0,0),c(1,1))
    matplot(ref,ref,type="l",col=1,lwd=2,ylab="Observed proportions",xlab="Predicted  probabilities",xlim=0:1,ylim=0:1)
    # smoother for calibration plots 
    ##################################
    # a = smoothing parameter
    if(isTRUE(smoothing)){a = smoothpar
    for(i in 1:k){p <- unlist(probs[[i]])
    points(smooth.spline(p, fitted(fitnp)[,i],spar=a), type="l", col=(1+i), lwd = 4)}}
    # legend
    for(i in 1:k){if(i <= 3){legende <- c("no LDA","LDA", "AEs")}
      if(i > 3){legende <- c(legende,paste("cat ", i, sep = ""))}}
    legend(x=0.7, y=(0.20+(k-3)*0.05),col=2:(k+1),lty =1,legend=legende)
    title(main = "Refitting: Non-parametric calibration plot")
    par(new=F)}
  
  
  ########################################
  # estimation of calibration intercepts 
  # cf. section 2.2.3. and 2.2.4.        
  ########################################
  
  if(isTRUE(intercept)){int<-vgam(outcome~1,offset=LP,family=multinomial(refLevel=r))
  coeffint<-coefficients(int)
  se<-sqrt(diag(vcov(int)))
  ci1i <- cbind(LL1 = coeffint[1] - qnorm(0.975) * se[1], UL1 = coeffint[1] + qnorm(0.975) * se[1])
  ci2i <- cbind(LL2 = coeffint[2] - qnorm(0.975) * se[2], UL2 = coeffint[2] + qnorm(0.975) * se[2])
  estint <- c(coeffint[1],ci1i,coeffint[2],ci2i)
  names(estint) <- paste('CALINT',c('int1','LLint1','ULint1','int2','LLint2','ULint2'),sep='.')}
  
  
  ####################################
  # estimation of calibration slopes 
  # cf. section 2.2.3. and 2.2.4.    
  ####################################
  
  # we used constraints to fix some coefficients to zero as appropriate
  # for k outcome categories this code should be changed to:
  # i <- diag(k-1)
  # i2 <- cbind(c(1,rep(0,k-2)))
  # i3 <- cbind(c(0,1,rep(0,k-1)))
  # i4 <- cbind(c(0,0,1,rep(0,k-2)))
  # ... (ij <- cbind(c(rep(0,j-2),1,rep(0,k-j)))
  # ik <- cbind(c(rep(0,k-2),1))
  # clist<-list("(Intercept)"=i,"lp1"=i2,"lp2"=i3,...,"lpk-1"=ik)
  # slopes<-vgam(outcome~lp1+lp2+...+lpk-1,family=multinomial(refLevel=r),constraints=clist)
  
  if(isTRUE(slope)){i<-diag(k-1)
  i2<-rbind(1,0)
  i3<-rbind(0,1)
  clist<-list("(Intercept)"=i,"lp1"=i2,"lp2"=i3)
  slopes<-vgam(outcome~lp1+lp2,family=multinomial(refLevel=r),constraints=clist)
  coeffslopes<-coefficients(slopes)[k:length(coefficients(slopes))]
  se<-sqrt(diag(vcov(slopes)))
  ci1s <- cbind(LL1 = coeffslopes[1] - qnorm(0.975) * se[3], UL1 = coeffslopes[1] + qnorm(0.975) * se[3])
  ci2s <- cbind(LL2 = coeffslopes[2] - qnorm(0.975) * se[4], UL2 = coeffslopes[2] + qnorm(0.975) * se[4])
  estslopes <- c(coeffslopes[1],ci1s,coeffslopes[2],ci2s)
  names(estslopes) <- paste('CALSLOPES',c('lp1','LLlp1','ULlp1','lp2','LLlp2','ULlp2'),sep='.')}
  
  
  #################################
  # calibration testing          
  # cf. section 2.2.3. and 2.2.4. 
  #################################
  
  # this code requires the bayesm library developed by Peter Rossi
  
  if(isTRUE(test)){
    
    # -2 log-likelihood of model without adaptations
    # for k outcome categories this code should be changed to:
    # alphas <- rep(0,k-1) #(i.e. all intercepts zero)
    # beta1 <- c(1,rep(0,k-2)) #(i.e. first linear predictor for first equation)
    # beta2 <- c(0,1,rep(0,k-3)) #(i.e. second linear predictor for second equation)      
    # betaj <- c(rep(0,j-1),1,rep(0,k-1-j)) #(i.e. jth linear predictor for jth equation)
    # betak <- c(rep(0,k-2),1) #(i.e. kth linear predictor for kth equation)
    # parametersk <- c(alphas, beta1, beta2, ..., betak)
    
    parametersk <- c(0,0,1,0,0,1) #c(alpha1,alpha2,b22,b23,b32,b33)
    Xdk=LP
    x <- createX(p=k,na=0,nd=k-1,Xa=NULL,Xd=Xdk,INT=TRUE,DIFF=FALSE,base=1)
    deviancewithout <- -2*llmnl(parametersk,outcome,x)
    names(deviancewithout)<-c('original deviance')
    
    devint <- deviance(int)
    names(devint)<-c('intercept deviance')
    devslopes <- deviance(slopes)
    names(devslopes)<-c('slopes deviance')
    
    # overall calibration (i.e. calibration intercepts and slopes) 
    ################################################################
    
    poverall<- pchisq(deviancewithout - devslopes, df = 2*(k-1), lower.tail = FALSE)
    
    # calibration intercepts 
    ##########################
    
    pint<- pchisq(deviancewithout - devint, df = k-1, lower.tail = FALSE)
    
    
    # calibration slopes 
    ######################
    
    pslopes<- pchisq(devint - devslopes, df = k-1, lower.tail = FALSE)
    names(poverall)<-c('p overall')
    names(pint)<-c('p int')
    names(pslopes)<-c('p slopes')}
  
  # Printing of results
  # The probabilities of calibration intercepts and slopes are only shown when the hypothesis of perfect 
  # calibration is rejected.
  
  results<-list(if(isTRUE(estimates)){est}else{'Not requested'},if(isTRUE(intercept)){estint}else{'Not requested'},if(isTRUE(slope)){estslopes}else{'Not requested'},if(isTRUE(test)){c(deviancewithout,devint,devslopes)}else{'Not requested'},if(isTRUE(test)){c(poverall,if(poverall<0.05){c(pint,pslopes)})}else{'Not requested'})
  names(results)<-c("Coefficients of parametric recalibration framework","Calibration Intercepts with 95% CI","Calibration Slopes with 95% CI","Deviances","P-values")
  n <- 1:5
  selection <- c(isTRUE(estimates),isTRUE(intercept),isTRUE(slope),isTRUE(test),isTRUE(test))
  results[n[selection]]
  
}

Polcal_results_0m_refit <- Polcal_refit(outcome=outcome, k=3, p=p_0m_refit, LP=LP_0m_refit, r=1, estimates=FALSE, dfr=2, plotoverall=TRUE, datapoints=TRUE, smoothing=TRUE, smoothpar=1, intercept=TRUE, slope=TRUE, test=FALSE)
Polcal_results_0m_refit

#once the calibration plots are produced, please save or screenshot these. There should be 2 non-parametric plots. 

#print results into dataframe
#calibration slope
results_refit[6,] <- Polcal_results_0m_refit$`Calibration Slopes with 95% CI`[[1]] #c-slope for y=2 vs 1
results_refit[7,] <- Polcal_results_0m_refit$`Calibration Slopes with 95% CI`[[2]] #lower CI for y=2 vs 1
results_refit[8,] <- Polcal_results_0m_refit$`Calibration Slopes with 95% CI`[[3]] #upper CI for y=2 vs 1
results_refit[9,] <- Polcal_results_0m_refit$`Calibration Slopes with 95% CI`[[4]] #c-slope for for y=3 vs 1
results_refit[10,] <- Polcal_results_0m_refit$`Calibration Slopes with 95% CI`[[5]] #lower CI for y=3 vs 1
results_refit[11,] <- Polcal_results_0m_refit$`Calibration Slopes with 95% CI`[[6]] #upper CI for y=3 vs 1
#same for calibration intercepts
results_refit[12,] <- Polcal_results_0m_refit$`Calibration Intercepts with 95% CI`[[1]]
results_refit[13,] <- Polcal_results_0m_refit$`Calibration Intercepts with 95% CI`[[2]]
results_refit[14,] <- Polcal_results_0m_refit$`Calibration Intercepts with 95% CI`[[3]]
results_refit[15,] <- Polcal_results_0m_refit$`Calibration Intercepts with 95% CI`[[4]]
results_refit[16,] <- Polcal_results_0m_refit$`Calibration Intercepts with 95% CI`[[5]]
results_refit[17,] <- Polcal_results_0m_refit$`Calibration Intercepts with 95% CI`[[6]]

#### Pairwise c-statistics
#outcome LDA (2) vs no LDA (1)
#filter out outcome 3
cstat_y2_0m_refit <- prob.cat.2 %>% filter(outcome != 3)
colnames(cstat_y2_0m_refit)[3] <- "outcome"
cstat_y2_0m_refit$outcome <- as.factor(cstat_y2_0m_refit$outcome)
#drop the unused outcome category
#cstat_y2_0m$outcome <- droplevels(cstat_y2_0m)
cstat_mod_y2_0m_refit <- roc(outcome~pr_test_mean,data=cstat_y2_0m_refit)
results_refit[18,] <- cstat_mod_y2_0m_refit$auc #AUC

se_refit <- sqrt(var(cstat_mod_y2_0m_refit))
cstat_y2_0m_LCI_refit <- (cstat_mod_y2_0m_refit$auc) - 1.96*se_refit 
results_refit[19,] <- cstat_y2_0m_LCI_refit #lower CI
cstat_y2_0m_UCI_refit <- (cstat_mod_y2_0m_refit$auc) + 1.96*se_refit 
results_refit[20,] <- cstat_y2_0m_UCI_refit #upper CI

#outcome AEs (3) vs no LDA (1)
#filter out outcome 2
cstat_y3_0m_refit <- prob.cat.3 %>% filter(outcome != 2)
colnames(cstat_y3_0m_refit)[3] <- "outcome"
cstat_y3_0m_refit$outcome <- as.factor(cstat_y3_0m_refit$outcome)
#drop the unused outcome category
cstat_y3_0m_refit$outcome <- droplevels(cstat_y3_0m_refit$outcome)
cstat_mod_y3_0m_refit <- roc(outcome~pr_test_mean,data=cstat_y3_0m_refit)
results_refit[21,] <- cstat_mod_y3_0m_refit$auc #AUC

se_refit <- sqrt(var(cstat_mod_y3_0m_refit))
cstat_y3_0m_LCI_refit <- (cstat_mod_y3_0m_refit$auc) - 1.96*se_refit 
results_refit[22,] <- cstat_y3_0m_LCI_refit #lower CI
cstat_y3_0m_UCI_refit <- (cstat_mod_y3_0m_refit$auc) + 1.96*se_refit
results_refit[23,] <- cstat_y3_0m_UCI_refit #upper CI

#PDI
prob1_refit <- prob.cat.1$pr_test_mean
outcome <- as.data.frame(outcome)
pdi_refit <- pdi3cat0(outcome, prob1_refit)
pdi_refit 

pdi3_refit <- pdi3v0(outcome, p_0m_refit)
results_refit[24,] <- pdi3_refit$pdi

#Variance explained - R2
#regress LPs onto outcome
R2_data <- as.data.frame(cbind(LP_0m_refit, data$multi_outcome))
colnames(R2_data)[3] <- "outcome"
citl_model_0m_refit <- multinom(outcome ~ (LP2vs1 + LP3vs1), data=R2_data) 

LR_0m_refit <- -2 * (as.numeric(logLik(multinom(outcome ~ 1, data=R2_data))) - 
                       as.numeric(logLik(citl_model_0m_refit)))

R2_coxsnell_0m_refit <- 1 - exp(-LR_0m_refit/length(R2_data$outcome))
R2_coxsnell_0m_refit 

#Nagelkerke calculations - adjusting the cox snell r2.
Nagel_0m_refit <- prev_y1_0m * log(prev_y1_0m/n_0m) + (n_0m-prev_y1_0m) * (log(1-(prev_y1_0m/n_0m)))
R2_0m_refit <- 1-exp((2*Nagel_0m_refit)/n_0m)
final_R2_refit <- R2_coxsnell_0m_refit/R2_0m_refit
results_refit[25,] <- final_R2_refit

#label results dataframe
rownames(results_refit) <- c("AIC_refit", "mean_LP2_0m", "sd_LP2_0m", "mean_LP3_0m", "sd_LP3_0m", "cslope_y2_0m", "cslope_LCI_y2_0m", "cslope_UCI_y2_0m", 
                             "cslope_y3_0m", "cslope_LCI_y3_0m", "cslope_UCI_y3_0m", "cint_y2_0m", "cint_LCI_y2_0m", "cint_UCI_y2_0m", "cint_y3_0m", "cint_LCI_y3_0m", 
                             "cint_UCI_y3_0m", "cstat_y2_0m", "cstat_LCI_y2_0m", "cstat_UCI_y2_0m", "cstat_y3_0m", "cstat_LCI_y3_0m", "cstat_UCI_y3_0m", "PDI", "R2_0m")
#limit to 5 decimal places
results_refit <- as.data.frame(results_refit) %>% mutate(across(where(is.numeric), ~ round(., 5)))

write.csv(results_updt, file="recalibration_results.csv")
write.csv(results_refit, file="refitting_results.csv")