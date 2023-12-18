#########################
# Section 3.3.4 of manuscript
# Externally validating and updating multinomial prediction model for MTX outcomes in RA
#########################

#load required libraries
library(tidyverse)
library(VGAM)
library(bayesm)
library(pROC)
library(mice)
library(nnet)
library(stats)

#record sample size (for R2)
n_0m <- nrow(data)

#record outcome prevalence (for R2)
prev <- as.data.frame(table(data$multi_outcome))
prev_y1_0m <- prev[1,"Freq"]

#filter variables we need for multiple imputation and model fitting 
data <- data %>% dplyr::select(pat_ID, multi_outcome, sex, age, TJC28, SJC28, CRP, RF, HAQ, DAS28_CRP_0m, PGA)

##################################
# Original model developed in RAMS
##################################

#multiple imputation 
multi_imp_0m <- mice(data, m=10, seed = 28382)
#stack datasets
stacked_data_0m <- complete(multi_imp_0m, action = "long") 

#fit linear predictor (LP) across all imputed datasets
fit_LP_y2_0m <- stacked_data_0m %>% filter(.imp != 0) %>% mutate(stacked_LP_y2_0m = 1.636419978 - (0.309527966*sex) + (0.009278724*age) 
                                                                 + (0.418708442*RF) - (0.580176307*HAQ) - (0.371220631*DAS28_CRP_0m))
#take average LP for each patient across the 10 imputed datasets
LP_y2_0m <- fit_LP_y2_0m %>% group_by(pat_ID) %>% summarise(mean_LP2 = mean(stacked_LP_y2_0m))
LP_y2_0m <- LP_y2_0m[ ,2]

#fit LP across all imputed datasets
fit_LP_y3_0m <- stacked_data_0m %>% filter(.imp != 0) %>% mutate(stacked_LP_y3_0m = -0.974627064 + (0.168876117*sex) + (0.002849258*age) 
                                                                 + (0.207914561*RF) + (0.187888068*HAQ) - (0.227110461*DAS28_CRP_0m))

#take average LP for each patient across the 10 imputed datasets
LP_y3_0m <- fit_LP_y3_0m %>% group_by(pat_ID) %>% summarise(mean_LP3 = mean(stacked_LP_y3_0m))
LP_y3_0m <- LP_y3_0m[ ,2]

#results matrix
results <- matrix(nrow = 34, ncol = 1)

#mean and sd of linear predictor
results[1,] <- mean(LP_y2_0m$mean_LP2)
results[2,] <- sd(LP_y2_0m$mean_LP2)
results[3,] <- mean(LP_y3_0m$mean_LP3)
results[4,] <- sd(LP_y3_0m$mean_LP3)

#####################
# Plotting predictor distributions to compare case-mix
#####################

#Age
ggplot(data, aes(x=age)) + 
  geom_histogram(color="black", fill="white", bins = 15) +
  geom_vline(aes(xintercept=median(age)),
             color="blue", linetype="dashed", size=1) +
  xlab("Age (years)") + 
  ggtitle("Distribution of age in NOR-DMARD")
#HAQ
haq_plot <- filter(data, !is.na(HAQ)) %>% select(HAQ)
ggplot(haq_plot, aes(x=HAQ)) + 
  geom_histogram(color="black", fill="white", bins=14) +
  geom_vline(aes(xintercept=median(HAQ)),
             color="blue", linetype="dashed", size=1) +
  xlab("HAQ score") + 
  ggtitle("Distribution of HAQ in NOR-DMARD") 
#DAS28
DAS28_plot <- filter(data, !is.na(DAS28_CRP_0m)) %>% select(DAS28_CRP_0m)
ggplot(DAS28_plot, aes(x=DAS28_CRP_0m)) + 
  geom_histogram(color="black", fill="white", bins=15) +
  geom_vline(aes(xintercept=median(DAS28_CRP_0m)),
             color="blue", linetype="dashed", size=1) +
  xlab("DAS28 score") + 
  ggtitle("Distribution of DAS28 in NOR-DMARD")

########
# Distribution of Linear Predictors
ggplot(LP_y2_0m, aes(mean_LP2)) + 
  geom_histogram(color="black", fill="white", bins = 20) +
  geom_vline(aes(xintercept=median(mean_LP2)),
             color="blue", linetype="dashed", size=1) +
  xlab("Linear Predictor y2 vs y1") + 
  ggtitle("Distribution of linear Predictor y2 vs y1 in RAMS")

ggplot(LP_y3_0m, aes(mean_LP3)) + 
  geom_histogram(color="black", fill="white", bins = 20) +
  geom_vline(aes(xintercept=median(mean_LP3)),
             color="blue", linetype="dashed", size=1) +
  xlab("Linear Predictor y3 vs y1") + 
  ggtitle("Distribution of linear Predictor y3 vs y1 in NOR-DMARD")

#####################
# Univariable associations between predictors and outcome
#####################

univ_age <- summary(multinom(multi_outcome ~ age, data = data))
results[5,] <- univ_age$coefficients[3]
results[6,] <- univ_age$coefficients[4]

univ_sex <- summary(multinom(multi_outcome ~ sex, data = data))
results[7,] <- univ_sex$coefficients[3]
results[8,] <- univ_sex$coefficients[4]

univ_haq <- summary(multinom(multi_outcome ~ HAQ, data = data))
results[9,] <- univ_haq$coefficients[3]
results[10,] <- univ_haq$coefficients[4]

univ_das28 <- summary(multinom(multi_outcome ~ DAS28_CRP_0m, data = data))
results[11,] <- univ_das28$coefficients[3]
results[12,] <- univ_das28$coefficients[4]

univ_rf <- summary(multinom(multi_outcome ~ RF, data = data))
results[13,] <- univ_rf$coefficients[3]
results[14,] <- univ_rf$coefficients[4]

#####-----------------------------------------------------------------

#defining predicted probabilities
prob_y1_0m <- (1 / (1 + exp(LP_y2_0m) + exp(LP_y3_0m)))
prob_y2_0m <- exp(LP_y2_0m) / (1 + exp(LP_y2_0m) + exp(LP_y3_0m))
prob_y3_0m <- exp(LP_y3_0m) / (1 + exp(LP_y2_0m) + exp(LP_y3_0m))

p_data_0m <- as.data.frame(cbind(prob_y1_0m, prob_y2_0m, prob_y3_0m, data$multi_outcome))
colnames(p_data_0m)[1] <- "prob(cat.1)"
colnames(p_data_0m)[2] <- "prob(cat.2)"
colnames(p_data_0m)[3] <- "prob(cat.3)"
colnames(p_data_0m)[4] <- "outcome"

p_0m <- p_data_0m %>% dplyr::select("prob(cat.1)", "prob(cat.2)", "prob(cat.3)") %>% as.matrix(.) #needs to be a matrix to work in the function below

#creating the LP matrix for calibration function below
LP_data_0m <- as.data.frame(cbind(LP_y2_0m, LP_y3_0m, multi_outcome))
colnames(LP_data_0m)[1] <- "LP2vs1"
colnames(LP_data_0m)[2] <- "LP3vs1"
colnames(LP_data_0m)[3] <- "outcome"

LP_0m <- LP_data_0m %>% dplyr::select(LP2vs1, LP3vs1) %>% as.matrix(.) #needs to be a matrix to work in the function below

outcome <- as.numeric(data$multi_outcome) 

#Assessing calibration of a multinomial prediction model - code from van Hoorde et al. https://doi.org/10.1002/sim.6114
Polcal <- function(outcome,k,p,LP,r=1,estimates=FALSE,dfr=2,plotoverall=TRUE,datapoints=TRUE,smoothing=TRUE,smoothpar=1,intercept=FALSE,slope=FALSE,test=FALSE){
  
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
  title(main = "External: Calibration plot")
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
    title(main = "External: Calibration plot")
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

Polcal_results_0m <- Polcal(outcome=outcome, k=3, p=p_0m, LP=LP_0m, r=1, estimates=FALSE, dfr=2, plotoverall=TRUE, datapoints=TRUE, smoothing=TRUE, smoothpar=1, intercept=TRUE, slope=TRUE, test=FALSE)
Polcal_results_0m

#Please save or screenshot these calibration plots. There should be 2 non-parametric plots. 

#print results into dataframe
#calibration slope
results[15,] <- Polcal_results_0m$`Calibration Slopes with 95% CI`[[1]] #c-slope for y=2 vs 1
results[16,] <- Polcal_results_0m$`Calibration Slopes with 95% CI`[[2]] #lower CI for y=2 vs 1
results[17,] <- Polcal_results_0m$`Calibration Slopes with 95% CI`[[3]] #upper CI for y=2 vs 1
results[18,] <- Polcal_results_0m$`Calibration Slopes with 95% CI`[[4]] #c-slope for for y=3 vs 1
results[19,] <- Polcal_results_0m$`Calibration Slopes with 95% CI`[[5]] #lower CI for y=3 vs 1
results[20,] <- Polcal_results_0m$`Calibration Slopes with 95% CI`[[6]] #upper CI for y=3 vs 1
#same for calibration intercepts
results[21,] <- Polcal_results_0m$`Calibration Intercepts with 95% CI`[[1]]
results[22,] <- Polcal_results_0m$`Calibration Intercepts with 95% CI`[[2]]
results[23,] <- Polcal_results_0m$`Calibration Intercepts with 95% CI`[[3]]
results[24,] <- Polcal_results_0m$`Calibration Intercepts with 95% CI`[[4]]
results[25,] <- Polcal_results_0m$`Calibration Intercepts with 95% CI`[[5]]
results[26,] <- Polcal_results_0m$`Calibration Intercepts with 95% CI`[[6]]

#### Pairwise c-statistics
#outcome LDA (2) vs no LDA (1)
#filter out outcome 3
cstat_y2_0m <- p_data_0m %>% filter(outcome != 3) %>%
  rename(prob1 = "prob(cat.1)",
         prob2 = "prob(cat.2)",
         prob3 = "prob(cat.3)")
cstat_y2_0m$outcome <- as.factor(cstat_y2_0m$outcome)
#drop the unused outcome category
#cstat_y2_0m$outcome <- droplevels(cstat_y2_0m)
cstat_mod_y2_0m <- roc(outcome~prob2,data=cstat_y2_0m)
results[27,] <- cstat_mod_y2_0m$auc #AUC

se <- sqrt(var(cstat_mod_y2_0m))
cstat_y2_0m_LCI <- (cstat_mod_y2_0m$auc) - 1.96*se 
results[28,] <- cstat_y2_0m_LCI #lower CI
cstat_y2_0m_UCI <- (cstat_mod_y2_0m$auc) + 1.96*se 
results[29,] <- cstat_y2_0m_UCI #upper CI

#outcome AEs (3) vs no LDA (1)
#filter out outcome 2
cstat_y3_0m <- p_data_0m %>% filter(outcome != 2) %>%
  rename(prob1 = "prob(cat.1)",
         prob2 = "prob(cat.2)",
         prob3 = "prob(cat.3)")
cstat_y3_0m$outcome <- as.factor(cstat_y3_0m$outcome)
#drop the unused outcome category
cstat_y3_0m$outcome <- droplevels(cstat_y3_0m$outcome)
cstat_mod_y3_0m <- roc(outcome~prob3,data=cstat_y3_0m)
results[30,] <- cstat_mod_y3_0m$auc #AUC

se <- sqrt(var(cstat_mod_y3_0m))
cstat_y3_0m_LCI <- (cstat_mod_y3_0m$auc) - 1.96*se 
results[31,] <- cstat_y3_0m_LCI #lower CI
cstat_y3_0m_UCI <- (cstat_mod_y3_0m$auc) + 1.96*se
results[32,] <- cstat_y3_0m_UCI #upper CI

#PDI
prob1 <- prob_y1_0m
pdi3cat0 <- function(outcome, prob1){
  
  #----------------------
  #-- outcome: outcome vector (3 category, 1 is of interest)
  #-- probs: probs of the category 1
  #----------------------
  
  p1=prob1[outcome==1] ; n1=length(p1) ; #w1=weight[outcome==1] ;
  p2=prob1[outcome==2] ; n2=length(p2) ; #w2=weight[outcome==2] ;
  p3=prob1[outcome==3] ; n3=length(p3) ; #w3=weight[outcome==3] ;
  
  wk=ties=rep(0,length(p1))
  for (i in 1:n1){
    count2=round(sum(p2<p1[i]))
    count3=round(sum(p3<p1[i]))
    
    #--- ties ---
    ties2=sum(p2==p1[i])
    ties3=sum(p3==p1[i])
    
    tmp=0
    tmp=tmp + ( ties2*count3)/2
    tmp=tmp + (count2* ties3)/2
    
    tmp=tmp + ( ties2* ties3)/3
    
    wk[i]=round(count2)*round(count3) + tmp
    
    
  }
  pdii=sum(wk)/(round(n1)*round(n2)*round(n3))
}
outcome <- as.data.frame(outcome)
pdi <- pdi3cat0(outcome, prob1)
pdi 

# This function calculates a PDI for each outcome category
pdi3v0=function(outcome, probs){
  
  aa=function(xx){
    if(xx==1){
      y=outcome ;
      prob1=probs[,1]
    }
    if(xx==2){
      y=outcome-1 ;
      y[outcome==1]=3 ;
      prob1=probs[,2]
    }
    if(xx==3){
      y=outcome-2 ;
      y[outcome==1]=2 ;
      y[outcome==2]=3 ;
      prob1=probs[,3]
    }
    Z=pdi3cat0(y, prob1)
    Z  
  }
  
  a=sapply(1:3, aa)
  
  w=table(outcome)
  
  Z=list()
  Z$pdii=a
  Z$pdi=mean(a)
  Z$pdiw=t(as.matrix(w))%*%as.matrix(a)/sum(as.matrix(w))
  return(Z)
}

pdi3 <- pdi3v0(outcome, p_0m)
results[33,] <- pdi3$pdi

#Variance explained - R2
#regress LPs onto outcome
citl_model_0m <- multinom(outcome ~ (LP2vs1 + LP3vs1), data=LP_data_0m) 

LR_0m <- -2 * (as.numeric(logLik(multinom(outcome ~ 1, data=LP_data_0m))) - 
                 as.numeric(logLik(citl_model_0m)))

R2_coxsnell_0m <- 1 - exp(-LR_0m/length(LP_data_0m$outcome))
R2_coxsnell_0m 

#Nagelkerke calculations - adjusting the cox snell r2.
Nagel_0m <- prev_y1_0m * log(prev_y1_0m/n_0m) + (n_0m-prev_y1_0m) * (log(1-(prev_y1_0m/n_0m)))
R2_0m <- 1-exp((2*Nagel_0m)/n_0m)
final_R2 <- R2_coxsnell_0m/R2_0m
results[34,] <- final_R2


#label results dataframe
rownames(results) <- c("mean_LP2_0m", "sd_LP2_0m", "mean_LP3_0m", "sd_LP3_0m", "univar_age_y2", "univar_age_y3", "univar_sex_y2", "univar_sex_y3", "univar_haq_y2", 
                       "univar_haq_y3", "univar_das_y2", "univar_das_y3", "univar_rf_y2", "univar_rf_y3", "cslope_y2_0m", "cslope_LCI_y2_0m", "cslope_UCI_y2_0m", 
                       "cslope_y3_0m", "cslope_LCI_y3_0m", "cslope_UCI_y3_0m", "cint_y2_0m", "cint_LCI_y2_0m", "cint_UCI_y2_0m", "cint_y3_0m", "cint_LCI_y3_0m", 
                       "cint_UCI_y3_0m", "cstat_y2_0m", "cstat_LCI_y2_0m", "cstat_UCI_y2_0m", "cstat_y3_0m", "cstat_LCI_y3_0m", "cstat_UCI_y3_0m", "PDI", "R2_0m")
#limit to 5 decimal places
results <- as.data.frame(results) %>% mutate(across(where(is.numeric), ~ round(., 5)))

#save each dataframe with results as excel document
write.csv(results, "main_results.xlsx")

