#########################
# Section 3.2.2 to 3.2.5 of manuscript
# Developing and internally validating (bootstrap) a multinomial prediction model
#########################

library(mice)
library(tidyverse)
library(pROC)
library(purrr)
library(labelled)
library(MASS)
library(nnet)
library(bayesm)
library(VGAM)
library(mcca)

#create dataframe of variables we will use in our model
predictors <- analysis_data %>% dplyr::select(outcome_multi_comb_LDA, sex, age_consent, tjc28_value, sjc28_value, DAS28_crp_value_cleaned, cfelab_rf_cat_dot, haq_score_sdi_ac, DAS28_CRP_0m, vaspat)
  
colnames(predictors) <- c("study_id", "outcome", "sex", "age", "TJC28", "SJC28", "CRP",
                          "RF", "haq", "DAS28", "PGA")

table(predictors$outcome)
predictors$outcome <- as.factor(predictors$outcome)
predictors$outcome <- relevel(predictors$outcome, ref = "1")

# Univariable associations (these will be compared between development and validation data)

univar_age <- summary(multinom(outcome ~ age, data = predictors))
univar_age$coefficients[3]#0.007342036
univar_age$coefficients[4]#0.003186375

univar_sex <- summary(multinom(outcome ~ sex, data = predictors))
univar_sex$coefficients[3]#-0.5724117
univar_sex$coefficients[4]#0.1788559

univar_haq <- summary(multinom(outcome ~ haq, data = predictors))
univar_haq$coefficients[3]#-0.8952755
univar_haq$coefficients[4]#0.0110848

univar_DAS28 <- summary(multinom(outcome ~ DAS28, data = predictors))
univar_DAS28$coefficients[3]#-0.5306533
univar_DAS28$coefficients[4]#-0.2152469

univar_RF <- summary(multinom(outcome ~ RF, data = predictors))
univar_RF$coefficients[3]#0.3113714
univar_RF$coefficients[4]#0.2258929

###########################
# MI & MLR REGRESSION 
########################
#impute data set
multi.imp <- mice(predictors, m=10, seed = 28382)
#stack datasets
stacked_data_og <- complete(multi.imp, action = "long") 

# fitting the multinomial model to the imputed dataset
multi_mo <- with(multi.imp, 
                   multinom(outcome ~  sex + age + RF + haq + DAS28, model = T))

### If interested in using data driven variable selection techniques, see code below

#FinalModel <- with(multi.imp, 
##                   expression(FullModel <- multinom(outcome ~ sex + age + RF + haq + ESR + 
#                                                 DAS28 + TJC28 + SJC28 + CRP),
##                              NullModel <- multinom(outcome ~ sex + age + RF + haq + DAS28),
#                               ModelSelection <- step(FullModel, k=2, scope = list("upper" = FullModel,
#                                                                                   "lower" = NullModel),
#                                                     direction = "backward",
#                                                     trace = TRUE)))
#Select the variables that are selected a majority of time across the imputed datasets
#Variables_Slected_NuM <- FinalModel$analyses %>%
#  map(formula) %>%
#  map(terms) %>%
#  map(labels) %>%
#  unlist() %>%
#  table()
#Variables_Slected <- names(Variables_Slected_NuM[which(Variables_Slected_NuM >= 10)]) #select majority

#FinalModel <- with(multi.imp,
#                   expression(form <- as.formula(paste("outcome ~", paste(Variables_Slected, collapse = "+"))),
#                              Mod <- multinom(formula = form)))


### summarising the model by pooling estimates across the 10 imputed datasets
pooled_sum <- summary(pool(multi_mo)) %>% 
              rename(outcome = y.level) %>%
              mutate(odds.ratio = exp(pooled_sum$estimate),
                     confint_lower = pooled_sum$odds.ratio * exp(-1.96 * pooled_sum$std.error),
                     confint_ipper = pooled_sum$odds.ratio * exp(1.96 * pooled_sum$std.error))

### Calculating probabilities and Linear predictors (LP)
#filtering out estimates and covariates for outcome 1 (vs 2)
ests2 <- pooled_sum %>% 
         data.frame("term" = pooled_sum$term,  "estimate" = pooled_sum$estimate) %>% 
         filter(outcome==2)
form <- as.formula(outcome ~ sex + age + RF + haq + DAS28)
Design_matrix <- model.matrix(form, data= stacked_data_og)

LP2 <- as.numeric(cbind(1,
                       Design_matrix[,ests2$term[-1]]) %*% ests2$estimate)

#filtering out estimates and covariates for outcome 3 (vs 1)
ests3 <- pooled_sum %>% 
         data.frame("term" = pooled_sum$term,  "estimate" = pooled_sum$estimate) %>% 
         filter(outcome==3)

LP3 <- as.numeric(cbind(1,
                       Design_matrix[,ests3$term[-1]]) %*% ests3$estimate)

#### Calculating probabilities for each outcome using the LP's defined above
#probabilities for outcome 1
OutSample_Predictions1 <- stacked_data_og %>%
  filter(.imp != 0) %>%
  mutate("Pi1" = 1 / (1 + exp(LP2) + exp(LP3)))

#probabilities for y=2
OutSample_Predictions2 <- stacked_data_og %>%
  filter(.imp != 0) %>%
  mutate("LP2" = LP2,
         "Pi2" = exp(LP2) / (1+ exp(LP2) + exp(LP3)))

#probabilities for y=3
OutSample_Predictions3 <- stacked_data_og %>%
  filter(.imp != 0) %>%
  mutate("LP3" = LP3,
         "Pi3" = exp(LP3) / (1+ exp(LP2) + exp(LP3)))

### summarising probabilities and LPs across the 10 imputed datasets
# for outcome 1
prob.cat.1 <- OutSample_Predictions1 %>% 
              group_by(.id) %>% 
              summarise(pr_test_mean = mean(Pi1), unique(outcome))

# for outcome 2
prob.cat.2 <- OutSample_Predictions2 %>% 
              group_by(.id) %>% 
              summarise(pr_test_mean = mean(Pi2), unique(outcome)) 

LP2vs1 <- OutSample_Predictions2 %>% 
          group_by(.id) %>% 
          summarise(LP2vs1 = mean(LP2), unique(outcome)) 

# for outcome 3
prob.cat.3 <- OutSample_Predictions3 %>% 
              group_by(.id) %>% 
              summarise(pr_test_mean = mean(Pi3), unique(outcome)) 

LP3vs1 <- OutSample_Predictions3 %>% 
          group_by(.id) %>% 
          summarise(LP3vs1 = mean(LP3), unique(outcome)) 

#creating the probability matrix for the calibration function
p <- as.data.frame(cbind(prob.cat.1$pr_test_mean, prob.cat.2$pr_test_mean, prob.cat.3$pr_test_mean)) %>%
      rename("prob(cat.1)" = V1,
             "prob(cat.2)" = V2,
             "prob(cat.3)" = V3)
          
# check that all probabilities add up to 1
table(p$`prob(cat.1)` + p$`prob(cat.2)` + p$`prob(cat.3)`)

p <- as.matrix((p)) #needs to be a matrix to work in the function below

#creating the LP matrix for calibration function below
LP <- as.data.frame(cbind(LP2vs1$LP2vs1, LP3vs1$LP3vs1)) %>%
      rename(LP2vs1 = V1,
             LP3vs1 = V2)

LP <- as.matrix(LP) #needs to be a matrix to work in the function below

#Assessing calibration of a multinomial prediction model - code from van Hoorde et al. https://doi.org/10.1002/sim.6114
##################################################################################
# outcome = column vector containing the outcome for every case, with values 1 to k (i.c. k=3)
outcome <- as.numeric(predictors$outcome)
# k = number of outcome categories (i.c. 3)
# p = matrix with the probabilities of the prediction model, ordered from prob(cat. 1) to prob(cat. k)
# LP = matrix with all the linear predictors with respect to the chosen reference category, ordered (e.g. LP2vs1 and LP3vs1)
# r = reference category (default: category 1)
# estimates = indicates whether the coefficients of the parametric recalibration framework are desired (default=FALSE)
# dfr =	degrees of freedom for the non-parametric calibration (default=2)
# plotoverall = indicates whether overall (non-)parametric calibration plots are constructed (default=TRUE)
# datapoints = indicates whether the individual datapoints are shown in the overall (non-)parametric calibration plots (default = TRUE)
# smoothing = indicates whether a smoothed line (using cubic splines) is added to the calibration plots (default=TRUE)
# smoothpar = smoothing parameter for the smoothed line (default=1)
# intercept = indicates whether calibration intercepts are desired (default=FALSE)
# slope = indicates whether calibration slopes are desired (default=FALSE)
# test = indicates whether statistical tests for calibration are desired (default=FALSE)
##################################################################################
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
  par(mfrow=c(ceiling(k/2),2))
  
  # parametric calibration plots 
  # cf. section 2.2.2.           
  ################################
  
  for(i in 1:k){p <- unlist(probs[[i]])
  if(isTRUE(smoothing)){color<-'grey'}else{color<-1+i}
  matplot(p,fitted(fitp)[,i],type="p",pch=i,col=color,lwd=1,ylab="",xlab="",xlim=0:1,ylim=0:1)
  par(new=T)
  ref <- rbind(c(0,0),c(1,1))
  matplot(ref,ref,type="l",col=1,lwd=2,ylab="Observed proportions",xlab="Predicted probabilities",xlim=0:1,ylim=0:1)
  # smoother for calibration plots 
  ##################################
  # a = smoothing parameter
  if(isTRUE(smoothing)){
    a = smoothpar
    points(smooth.spline(p, fitted(fitp)[,i],spar=a), type="l", col=(1+i), lwd = 4)}
  # legend
  legende <- c(paste("cat ", i, sep = ""))
  legend(x=0.6, y=(0.2),col=(1+i),lty =1,legend=legende)
  title(main = "Parametric calibration plot")
  par(new=F)}
  
  # non-parametric calibration plot 
  # cf. section 2.2.2.              
  ###################################
  windows()
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
  title(main = "Internal: Non-parametric calibration plot")
  par(new=F)}
  
  
  #############################################            
  # Overall (non-)parametric calibration plot 
  #############################################
  
  if(isTRUE(plotoverall)){
    windows()
    
    # parametric calibration plot 
    # cf. section 2.2.2.          
    ###############################
    
    if(isTRUE(datapoints)){for(i in 1:k){p <- unlist(probs[[i]])
    matplot(p,fitted(fitp)[,i],type="p",pch=i,col=(1+i),lwd=1,ylab="",xlab="",xlim=0:1,ylim=0:1)
    par(new=T)}}
    ref <- rbind(c(0,0),c(1,1))
    matplot(ref,ref,type="l",col=1,lwd=2,ylab="Observed proportions",xlab="Predicted probabilities",xlim=0:1,ylim=0:1)
    # smoother for calibration plots 
    ##################################
    # a = smoothing parameter
    if(isTRUE(smoothing)){
      a = smoothpar
      for(i in 1:k){p <- unlist(probs[[i]])
      points(smooth.spline(p, fitted(fitp)[,i],spar=a), type="l", col=(1+i), lwd = 4)}}
    # legend
    for(i in 1:k){if(i <= 3){legende <- c("no LDA", "LDA", "AEs")}
      if(i > 3){legende <- c(legende,paste("cat ", i, sep = ""))}}
    legend(x=0.7, y=(0.20+(k-3)*0.05),col=2:(k+1),lty =1,legend=legende)
    title(main = "Parametric calibration plot")
    par(new=F)
    
    # non-parametric calibration plot 
    # cf. section 2.2.2.              
    ###################################
    
    windows()
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
    for(i in 1:k){if(i <= 3){legende <- c("no LDA", "LDA", "AEs")}
      if(i > 3){legende <- c(legende,paste("cat ", i, sep = ""))}}
    legend(x=0.6, y=(0.20+(k-3)*0.05),col=2:(k+1),lty =1,legend=legende)
    title(main = "Internal: Non-parametric calibration plot")
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

Polcal.results <- Polcal(outcome = outcome,k = 3,p = p,LP = LP,r=1,estimates=TRUE,dfr=2,plotoverall=TRUE,datapoints=TRUE,smoothing=TRUE,smoothpar=1,intercept=TRUE,slope=TRUE,test=FALSE)
Polcal.results
recal <- Polcal.results$`Coefficients of parametric recalibration framework`

LP2new <- recal[1] + (recal[2]*LP$LP2vs1) + (recal[3]*LP$LP3vs1)
LP3new <- recal[3] + (recal[5]*LP$LP2vs1) + (recal[6]*LP$LP3vs1)

### PDI
prob1 <- OutSample_Predictions1$Pi1
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

pdi3 <- pdi3v0(outcome, p)
pdi3
pdi3$pdi
#outcome 1: 0.490
#outcome 2: 0.590
#outcome 3: 0.472
#overall PDI: 0.515

#Pairwise c-statistics
#Submodel 1
#filter out outcome 3
cstat.2 <- prob.cat.2 %>% filter(outcome != 3)
colnames(cstat.2)[3] <- "outcome"
#drop the unused outcome category
cstat.2$outcome <- droplevels(cstat.2$outcome)
cstat.model.2 <- roc(outcome~pr_test_mean,data=cstat.2)
cstat.model.2$auc #0.72
se <- sqrt(var(cstat.model.2))
cstat.2.LCI <- (cstat.model.2$auc) - 1.96*se 
cstat.2.LCI #0.70
cstat.2.UCI <- (cstat.model.2$auc) + 1.96*se 
cstat.2.UCI #0.75

#Submodel 21
#filter out outcome 2
cstat.3 <- prob.cat.3 %>% filter(outcome != 2)
colnames(cstat.3)[3] <- "outcome"
#drop the unused outcome category
cstat.3$outcome <- droplevels(cstat.3$outcome)
cstat.model.3 <- roc(outcome~pr_test_mean,data=cstat.3)
cstat.model.3$auc #0.56
se <- sqrt(var(cstat.model.3))
cstat.3.LCI <- (cstat.model.3$auc) - 1.96*se #0.51
cstat.3.LCI
cstat.3.UCI <- (cstat.model.3$auc) + 1.96*se #0.61
cstat.3.UCI

#Nagelkerke R2
LP_for_R2 <- as.data.frame(cbind(LP2vs1$LP2vs1, LP3vs1$LP3vs1, LP3vs1$`unique(outcome)`)) %>% 
              rename(LP2vs1 = V1,
                     LP3vs1 = V2,
                     outcome = V3)

#regress LPs onto outcome
citl_model <- multinom(outcome ~ (LP2vs1 + LP3vs1), data=LP_for_R2) 

LR <- -2 * (as.numeric(logLik(multinom(outcome ~ 1, data=LP_for_R2))) - 
              as.numeric(logLik(citl_model)))

R2_coxsnell <- 1 - exp(-LR/length(LP_for_R2$outcome))
R2_coxsnell #0.15

##########################
##########################
# Internal Validation
## BOOTSTRAP PROGRAMME
##########################
##########################

# Start with the original dataset, in which some data is missing
# boostrap from the original dataset
# do everything in the bootstrap dataset that we did to develop the main model in the original dataset
# a. impute 10 times
# b. calculate performance measures averaged across the 10 imputed datasets
# apply the model developed in this bootstrap dataset to the 10 multiple imputations of the og dataset that we used to develop the main model,
# averaging performance measures across the 10 imputed versions of the og dataset
# calculate the difference (optimism) in the performance measures produced in the bootstrap (development) in the original (test) dataset

manual_boot_bw <- function(data, outsample, samples){
  results <- matrix(nrow = samples,ncol = 32)
  set.seed(231398)
  for (i in 1:samples) {
    samp_index <- sample(1:nrow(data), nrow(data), rep=TRUE) # create a sampling index vector
    
    bs_samp <- data[samp_index,] # index the orignal dataset using the sampling vector to give the bs sample
    imp_bs <- mice(bs_samp, m=10) #impute bootstrap sample
    stacked_data <- complete(imp_bs, action = "long") #stack bootstrap datasets
    
    multi_mo <- with(imp_bs,  
                     multinom(outcome ~ sex + age + RF + haq + DAS28, model = T)) #fit multinomial model
    
    ### summarising the model by pooling estimates across the 10 imputed datasets
    pooled_sum <- summary(pool(multi_mo))
    colnames(pooled_sum)[1] <- "outcome" 
    
    ###---------- APPARENT MEASURES IN BOOTSTRAP SAMPLE
    
    # Calculating probabilities and Linear predictors (LP)
    #filtering out estimates and covariates for outcome 1 (vs 2)
    ests2 <- pooled_sum %>% data.frame("term" = pooled_sum$term,  "estimate" = pooled_sum$estimate) %>% filter(outcome==2)
    
    form <- as.formula(outcome ~ sex + age + RF + haq + DAS28)
    Boot_Design_matrix <- model.matrix(form, data= stacked_data)
    
    #Linear predictor for outcome 2 vs 1
    LP2 <- as.numeric(cbind(1,
                            Boot_Design_matrix[,ests2$term[-1]]) %*% ests2$estimate)
    
    #filtering out estimates and covariates for outcome 3 (vs 1)
    ests3 <- pooled_sum %>% data.frame("term" = pooled_sum$term,  "estimate" = pooled_sum$estimate) %>% filter(outcome==3)
    
    #Linear predictor for outcome 3 vs 1
    LP3 <- as.numeric(cbind(1,
                            Boot_Design_matrix[,ests3$term[-1]]) %*% ests3$estimate)
    
    #probabilities for y=1
    Boot_Predictions1 <- stacked_data %>%
      filter(.imp != 0) %>%
      mutate("Pi1" = 1 / (1 + exp(LP2) + exp(LP3)))
    
    #probabilities for y=2
    Boot_Predictions2 <- stacked_data %>%
      filter(.imp != 0) %>%
      mutate("LP2" = LP2,
             "Pi2" = exp(LP2) / (1+ exp(LP2) + exp(LP3)))
    
    #probabilities for y=3
    Boot_Predictions3 <- stacked_data %>%
      filter(.imp != 0) %>%
      mutate("LP3" = LP3,
             "Pi3" = exp(LP3) / (1+ exp(LP2) + exp(LP3)))
    
    # summarising across all 10 imputed datasets
    # for outcome 1
    boot.prob.cat.1 <- Boot_Predictions1 %>% 
                       group_by(.id) %>% 
                       summarise(pr_bs_mean = mean(Pi1), unique(outcome))
    # for outcome 2
    boot.prob.cat.2 <- Boot_Predictions2 %>% 
                       group_by(.id) %>% 
                       summarise(pr_bs_mean = mean(Pi2), unique(outcome)) 
    
    boot.LP2vs1 <- Boot_Predictions2 %>% 
                   group_by(.id) %>% 
                   summarise(LP2vs1 = mean(LP2), unique(outcome)) 
    
    # for outcome 3
    boot.prob.cat.3 <- Boot_Predictions3 %>% 
                       group_by(.id) %>% 
                       summarise(pr_bs_mean = mean(Pi3), unique(outcome)) 
    
    boot.LP3vs1 <- Boot_Predictions3 %>% 
                   group_by(.id) %>% 
                  summarise(LP3vs1 = mean(LP3), unique(outcome)) 
    
    #creating the probability matrix for the calibration function
    p.bs <- as.data.frame(cbind(boot.prob.cat.1$pr_bs_mean, boot.prob.cat.2$pr_bs_mean, boot.prob.cat.3$pr_bs_mean)) %>%
      rename("prob(cat.1)" = V1,
             "prob(cat.2)" = V2,
             "prob(cat.3)" = V3)
    
    p.bs <- as.matrix((p.bs))
    
    #creating the LP matrix for calibration function below
    LP.bs <- as.data.frame(cbind(boot.LP2vs1$LP2vs1, boot.LP3vs1$LP3vs1)) %>% 
            rename("LP2vs1" = V1,
                   "LP3vs1" = V2)

    LP.bs <- as.matrix(LP.bs)
    
    outcome <- as.numeric(bs_samp$outcome)
  
    #PDI measure
    pdi3v0 <- function(outcome, probs){
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
    pdi3.bs <- pdi3v0(outcome, p.bs)
    results[i,1] <- pdi3.bs$pdi
    
    ### CALIBRATION
    #this function is slightly adjusted from the Polcal function on line 148-407 as we are not interested in calibration plots
    Polcal2 <- function(outcome,k,p,LP,r=1,estimates=FALSE,dfr=2,intercept=FALSE,slope=FALSE,test=FALSE,...){
      if(isTRUE(test)){intercept<-slope<-TRUE}
      # probabilities
      probs <- split(p,col(p))    
      # linear predictors necessary for non-parametric calibration plot - give a name to each linear predictor 
      # seperately
      lps <- split(LP,col(LP))
      for(i in 1:(k-1)){assign(paste("lp", i, sep = ""),unlist(lps[[i]]))}
      fitp<-vglm(outcome~LP,family=multinomial(refLevel=r))
      if(isTRUE(estimates)){est<-coefficients(fitp)
      names(est) <- paste('EST',names(est),sep='.')}
      fitnp<-vgam(outcome~s(lp1,df=dfr)+s(lp2,df=dfr),family=multinomial(refLevel=r))
      
      if(isTRUE(intercept)){int<-vgam(outcome~1,offset=LP,family=multinomial(refLevel=r))
      coeffint<-coefficients(int)
      se<-sqrt(diag(vcov(int)))
      ci1i <- cbind(LL1 = coeffint[1] - qnorm(0.975) * se[1], UL1 = coeffint[1] + qnorm(0.975) * se[1])
      ci2i <- cbind(LL2 = coeffint[2] - qnorm(0.975) * se[2], UL2 = coeffint[2] + qnorm(0.975) * se[2])
      estint <- c(coeffint[1],ci1i,coeffint[2],ci2i)
      names(estint) <- paste('CALINT',c('int1','LLint1','ULint1','int2','LLint2','ULint2'),sep='.')}
      
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
      
      if(isTRUE(test)){
        
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
        pint<- pchisq(deviancewithout - devint, df = k-1, lower.tail = FALSE)
        # calibration slopes 
        pslopes<- pchisq(devint - devslopes, df = k-1, lower.tail = FALSE)
        names(poverall)<-c('p overall')
        names(pint)<-c('p int')
        names(pslopes)<-c('p slopes')}
      
      # Printing of results
      results<-list(if(isTRUE(estimates)){est}else{'Not requested'},if(isTRUE(intercept)){estint}else{'Not requested'},if(isTRUE(slope)){estslopes}else{'Not requested'},if(isTRUE(test)){c(deviancewithout,devint,devslopes)}else{'Not requested'},if(isTRUE(test)){c(poverall,if(poverall<0.05){c(pint,pslopes)})}else{'Not requested'})
      names(results)<-c("Coefficients of parametric recalibration framework","Calibration Intercepts with 95% CI","Calibration Slopes with 95% CI","Deviances","P-values")
      n <- 1:5
      selection <- c(isTRUE(estimates),isTRUE(intercept),isTRUE(slope),isTRUE(test),isTRUE(test))
      results[n[selection]]
      
    }
    Polcal2.results <- Polcal2(outcome=outcome, k=3, p=p.bs, LP=LP.bs, r=1, estimates=FALSE, dfr=2, intercept=TRUE, slope=TRUE, test=FALSE)
    
    #Nagelkerke R2
    LP_for_R2_bs <- as.data.frame(cbind(boot.LP2vs1$LP2vs1, boot.LP3vs1$LP3vs1, boot.LP3vs1$`unique(outcome)`)) %>% 
      rename(LP2vs1 = V1,
             LP3vs1 = V2,
             outcome = V3)
    #regress LPs onto outcome
    citl_model_bs <- multinom(outcome ~ (LP2vs1 + LP3vs1), data=LP_for_R2_bs) 
    LR_bs <- -2 * (as.numeric(logLik(multinom(outcome ~ 1, data=LP_for_R2_bs))) - 
                  as.numeric(logLik(citl_model_bs)))
    R2_coxsnell_bs <- 1 - exp(-LR_bs/length(LP_for_R2_bs$outcome))
    
    #pairwise c-statistics
    #outcome 2 vs 1.
    #filter out outcome 3
    cstat.2 <- boot.prob.cat.2 %>% filter(outcome != 3)
    colnames(cstat.2)[3] <- "outcome"
    #drop the unused outcome category
    cstat.2$outcome <- droplevels(cstat.2$outcome)
    cstat.model.2.bs <- roc(outcome~pr_bs_mean,data=cstat.2)
    cstat.model.2.bs$auc[1] 
    
    #outcome 3 vs 1.
    #filter out outcome 2
    cstat.3 <- boot.prob.cat.3 %>% filter(outcome != 2)
    colnames(cstat.3)[3] <- "outcome"
    #drop the unused outcome category
    cstat.3$outcome <- droplevels(cstat.3$outcome)
    cstat.model.3.bs <- roc(outcome~pr_bs_mean,data=cstat.3)
    cstat.model.3.bs$auc 
    
    #entering calibration output into results matrix
    #calibration slopes and confidence interval
    results[i,2] <- Polcal2.results$`Calibration Slopes with 95% CI`[[1]] #c-slope for y=2 vs 1
    results[i,3] <- Polcal2.results$`Calibration Slopes with 95% CI`[[2]] #lower CI for y=2 vs 1
    results[i,4] <- Polcal2.results$`Calibration Slopes with 95% CI`[[3]] #upper CI for y=2 vs 1
    results[i,5] <- Polcal2.results$`Calibration Slopes with 95% CI`[[4]] #c-slope for for y=3 vs 1
    results[i,6] <- Polcal2.results$`Calibration Slopes with 95% CI`[[5]] #lower CI for y=3 vs 1
    results[i,7] <- Polcal2.results$`Calibration Slopes with 95% CI`[[6]] #upper CI for y=3 vs 1
    #same for calibration intercepts
    results[i,8] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[1]]
    results[i,9] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[2]]
    results[i,10] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[3]]
    results[i,11] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[4]]
    results[i,12] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[5]]
    results[i,13] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[6]]
    
    results[i,14] <- R2_coxsnell_bs
    
    results[i,15] <- cstat.model.2.bs$auc[1] 
    results[i,16] <- cstat.model.3.bs$auc[1]
    
    ####---------------- OUTSAMPLE MEASURES
    
    # this uses data on which the original model was developed
    # Calculating probabilities and Linear predictors (LP)
    Test_Design_matrix <- model.matrix(form, data = outsample)
    #Linear predictor for outcome 2 vs 1
    LP2 <- as.numeric(cbind(1,
                            Test_Design_matrix[,ests2$term[-1]]) %*% ests2$estimate)
    #Linear predictor for outcome 3 vs 1
    LP3 <- as.numeric(cbind(1,
                            Test_Design_matrix[,ests3$term[-1]]) %*% ests3$estimate)
    #probabilities for y=1
    Test_Predictions1 <- outsample %>%
      filter(.imp != 0) %>%
      mutate("Pi1" = 1 / (1 + exp(LP2) + exp(LP3)))
    
    #probabilities for y=2
    Test_Predictions2 <- outsample %>%
      filter(.imp != 0) %>%
      mutate("LP2" = LP2,
             "Pi2" = exp(LP2) / (1+ exp(LP2) + exp(LP3)))
    
    #probabilities for y=3
    Test_Predictions3 <- outsample %>%
      filter(.imp != 0) %>%
      mutate("LP3" = LP3,
             "Pi3" = exp(LP3) / (1+ exp(LP2) + exp(LP3)))
    
    ### summarising probabilities and LPs across the 10 imputed datasets
    #average patients indivual risks across the 10 imputed datasets
    #--- for test sample 
    # for outcome 1
    test.prob.cat.1 <- Test_Predictions1 %>% group_by(.id) %>% summarise(pr_test_mean = mean(Pi1), unique(outcome))
    # for outcome 2
    test.prob.cat.2 <- Test_Predictions2 %>% group_by(.id) %>% summarise(pr_test_mean = mean(Pi2), unique(outcome)) 
    test.LP2vs1 <- Test_Predictions2 %>% group_by(.id) %>% summarise(LP2vs1 = mean(LP2), unique(outcome)) 
    # for outcome 3
    test.prob.cat.3 <- Test_Predictions3 %>% group_by(.id) %>% summarise(pr_test_mean = mean(Pi3), unique(outcome)) 
    test.LP3vs1 <- Test_Predictions3 %>% group_by(.id) %>% summarise(LP3vs1 = mean(LP3), unique(outcome)) 
    
    
    #creating the probability matrix for the calibration function
    p.test <- as.data.frame(cbind(test.prob.cat.1$pr_test_mean, test.prob.cat.2$pr_test_mean, test.prob.cat.3$pr_test_mean)) %>%
      rename("prob(cat.1)" = V1,
             "prob(cat.2)" = V2,
             "prob(cat.3)" = V3)
    
    p.test <- as.matrix((p.test))
    
    #creating the LP matrix for calibration function below
    LP.test <- as.data.frame(cbind(test.LP2vs1$LP2vs1, test.LP3vs1$LP3vs1)) %>%
      rename("LP2vs1" = V1,
             "LP3vs1" = V2)

    LP.test <- as.matrix(LP.test)
  
    ### PDI and Calibration measures using function defined on lines 592-620 (PDI) and 625-690 (Calibration)
    #PDI measure
    outcome <- as.numeric(data$outcome)
    pdi3.test <- pdi3v0(outcome, p.test)
    results[i,17] <- pdi3.test$pdi
    #calibration
    Polcal2.results <- Polcal2(outcome=outcome, k=3, p=p.test, LP=LP.test, r=1, estimates=FALSE, dfr=2, intercept=TRUE, slope=TRUE, test=FALSE)
    
    #Nagelkerke R2
    LP_for_R2_test <- as.data.frame(cbind(test.LP2vs1$LP2vs1, test.LP3vs1$LP3vs1, test.LP3vs1$`unique(outcome)`)) %>% 
      rename(LP2vs1 = V1,
             LP3vs1 = V2,
             outcome = V3)
    #regress LPs onto outcome
    citl_model_test <- multinom(outcome ~ (LP2vs1 + LP3vs1), data=LP_for_R2_test) 
    LR_test <- -2 * (as.numeric(logLik(multinom(outcome ~ 1, data=LP_for_R2_test))) - 
                     as.numeric(logLik(citl_model_test)))
    R2_coxsnell_test <- 1 - exp(-LR_test/length(LP_for_R2_test$outcome))
    
    #pairwise c-statistics
    #outcome 2 vs 1.
    #filter out outcome 3
    cstat.2 <- test.prob.cat.2 %>% filter(outcome != 3)
    colnames(cstat.2)[3] <- "outcome"
    #drop the unused outcome category
    cstat.2$outcome <- droplevels(cstat.2$outcome)
    cstat.model.2.test <- roc(outcome~pr_test_mean,data=cstat.2)
    cstat.model.2.test$auc 
    
    #outcome 3 vs 1.
    #filter out outcome 2
    cstat.3 <- test.prob.cat.3 %>% filter(outcome != 2)
    colnames(cstat.3)[3] <- "outcome"
    #drop the unused outcome category
    cstat.3$outcome <- droplevels(cstat.3$outcome)
    cstat.model.3.test <- roc(outcome~pr_test_mean,data=cstat.3)
    cstat.model.3.test$auc 
    
    results[i,18] <- Polcal2.results$`Calibration Slopes with 95% CI`[[1]]
    results[i,19] <- Polcal2.results$`Calibration Slopes with 95% CI`[[2]]
    results[i,20] <- Polcal2.results$`Calibration Slopes with 95% CI`[[3]]
    results[i,21] <- Polcal2.results$`Calibration Slopes with 95% CI`[[4]]
    results[i,22] <- Polcal2.results$`Calibration Slopes with 95% CI`[[5]]
    results[i,23] <- Polcal2.results$`Calibration Slopes with 95% CI`[[6]]
    
    results[i,24] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[1]]
    results[i,25] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[2]]
    results[i,26] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[3]]
    results[i,27] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[4]]
    results[i,28] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[5]]
    results[i,29] <- Polcal2.results$`Calibration Intercepts with 95% CI`[[6]]
    
    results[i,30] <- R2_coxsnell_test
    
    results[i,31] <- cstat.model.2.test$auc[1] 
    results[i,32] <- cstat.model.3.test$auc[1]
    print(i)
  }
  results2 <- as.data.frame(results)
  colnames(results2) <- c("boot_PDI", "boot_cslope1", "boot_cslope1_LCI", "boot_cslope1_UCI", "boot_cslope3", 
                          "boot_cslope3_LCI", "boot_cslope3_UCI", "boot_intercept1", "boot_intercept1_LCI",
                          "boot_intercept1_UCI", "boot_intercept3", "boot_intercept3_LCI", "boot_intercept3_UCI", "boot_R2", "cstat2_bs", "cstat3_bs",
                          "test_PDI", "test_cslope1", "test_cslope1_LCI", "test_cslope1_UCI", "test_cslope3", 
                          "test_cslope3_LCI", "test_cslope3_UCI", "test_intercept1", "test_intercept1_LCI",
                          "test_intercept1_UCI", "test_intercept3", "test_intercept3_LCI", "test_intercept3_UCI", "test_R2", "cstat2_test", "cstat3_test")
  return(results2)

}
boot_results_bw <- manual_boot_bw(predictors, stacked_data_og, 1000) #using 50 samples just to get an indication of results, will use at least 500 in the real run

# data = predictors ---> original dataset that is not imputed - which will then get bootstrapped and imputed and used for model building
# outsample = stacked_data_og ---> outsample dataset or original dataset on which we want to test performance of the bootstrap model
# samples = number of bootstraps performed

#### optimism adjusted statistics
# Original measure - mean(boot - original)
pdi.adj <- pdi3$pdi - (mean(boot_results_bw$boot_PDI) - mean(boot_results_bw$test_PDI))
pdi.adj #0.49

#c-slope for y=2 vs 1
Calslope.adj.1 <- Polcal.results$`Calibration Slopes with 95% CI`[[1]] - (mean(boot_results_bw$boot_cslope1) - mean(boot_results_bw$test_cslope1))
Calslope.adj.1 #1.01 (0.87, 1.14)

#manual CI
boot_LCI_1 <- mean(boot_results_bw$boot_cslope1) - 1.96*(sd(boot_results_bw$boot_cslope1) / sqrt(1000))
test_LCI_1 <- mean(boot_results_bw$test_cslope1) - 1.96*(sd(boot_results_bw$test_cslope1) / sqrt(1000))
LCI_1 <- Polcal.results$`Calibration Slopes with 95% CI`[[2]] - (boot_LCI_1 - test_LCI_1)
LCI_1 #0.87

boot_UCI_1 <- mean(boot_results_bw$boot_cslope1) + 1.96*(sd(boot_results_bw$boot_cslope1) / sqrt(1000))
test_UCI_1 <- mean(boot_results_bw$test_cslope1) + 1.96*(sd(boot_results_bw$test_cslope1) / sqrt(1000))
UCI_1 <- Polcal.results$`Calibration Slopes with 95% CI`[[3]] - (boot_LCI_1 - test_LCI_1)
UCI_1 #1.14

#c-slope for y=3 vs 1
Calslope.adj.3 <- Polcal.results$`Calibration Slopes with 95% CI`[[4]] - (mean(boot_results_bw$boot_cslope3) - mean(boot_results_bw$test_cslope3))
Calslope.adj.3 #0.83 (0.30, 1.34)

#manual CI
boot_LCI_3 <- mean(boot_results_bw$boot_cslope3) - 1.96*(sd(boot_results_bw$boot_cslope3) / sqrt(1000))
test_LCI_3 <- mean(boot_results_bw$test_cslope3) - 1.96*(sd(boot_results_bw$test_cslope3) / sqrt(1000))
LCI_3 <- Polcal.results$`Calibration Slopes with 95% CI`[[5]] - (boot_LCI_3 - test_LCI_3)
LCI_3 #0.30

boot_UCI_3 <- mean(boot_results_bw$boot_cslope3) + 1.96*(sd(boot_results_bw$boot_cslope3) / sqrt(1000))
test_UCI_3 <- mean(boot_results_bw$test_cslope3) + 1.96*(sd(boot_results_bw$test_cslope3) / sqrt(1000))
UCI_3 <- Polcal.results$`Calibration Slopes with 95% CI`[[6]] - (boot_LCI_3 - test_LCI_3)
UCI_3 #1.34

#Cox-Snell R2
R2_adj <- R2_coxsnell - (mean(boot_results_bw$boot_R2) - mean(boot_results_bw$test_R2))
R2_adj #0.14

#Nagelkerke calculations - adjusting the cox snell r2.
L0.1 <- 730 * log(730/1632) + (1632-730) * (log(1-(730/1632)))
L0.1

1-exp((2*L0.1)/1632)

R2_adj/0.75 # 19%

# calibration intercept
cal_int_2 <- Polcal.results$`Calibration Intercepts with 95% CI`[[1]] - (mean(boot_results_bw$boot_intercept1) - mean(boot_results_bw$test_intercept1))
cal_int_2 #0.001 (-0.111, 0.114)

boot_cint2_LCI <- mean(boot_results_bw$boot_intercept1) - 1.96*(sd(boot_results_bw$boot_intercept1) / sqrt(1000))
test_cint2_LCI <- mean(boot_results_bw$test_intercept1) - 1.96*(sd(boot_results_bw$test_intercept1) / sqrt(1000))
LCI_1 <- Polcal.results$`Calibration Intercepts with 95% CI`[[2]] - (boot_cint2_LCI - test_cint2_LCI)
LCI_1 #-0.111
boot_cint2_UCI <- mean(boot_results_bw$boot_intercept1) + 1.96*(sd(boot_results_bw$boot_intercept1) / sqrt(1000))
test_cint2_UCI <- mean(boot_results_bw$test_intercept1) + 1.96*(sd(boot_results_bw$test_intercept1) / sqrt(1000))
UCI_1 <- Polcal.results$`Calibration Intercepts with 95% CI`[[3]] - (boot_cint2_UCI - test_cint2_UCI)
UCI_1 #0.114

cal_int_3 <- Polcal.results$`Calibration Intercepts with 95% CI`[[4]] - (mean(boot_results_bw$boot_intercept3) - mean(boot_results_bw$test_intercept3))
cal_int_3 #0.004 (-0.180, 0.188)

boot_cint3_LCI <- mean(boot_results_bw$boot_intercept3) - 1.96*(sd(boot_results_bw$boot_intercept3) / sqrt(1000))
test_cint3_LCI <- mean(boot_results_bw$test_intercept3) - 1.96*(sd(boot_results_bw$test_intercept3) / sqrt(1000))
LCI_1 <- Polcal.results$`Calibration Intercepts with 95% CI`[[5]] - (boot_cint3_LCI - test_cint3_LCI)
LCI_1 #-0.180
boot_cint3_UCI <- mean(boot_results_bw$boot_intercept3) + 1.96*(sd(boot_results_bw$boot_intercept3) / sqrt(1000))
test_cint3_UCI <- mean(boot_results_bw$test_intercept3) + 1.96*(sd(boot_results_bw$test_intercept3) / sqrt(1000))
UCI_1 <- Polcal.results$`Calibration Intercepts with 95% CI`[[6]] - (boot_cint3_UCI - test_cint3_UCI)
UCI_1 #0.188

# pairwise c-statistics
cstat_2vs1_adj <- cstat.model.2$auc[1] - (mean(boot_results_bw$cstat2_bs) - mean(boot_results_bw$cstat2_test))
cstat_2vs1_adj #0.72 (0.69, 0.74)

boot_cstat2_LCI_1 <- mean(boot_results_bw$cstat2_bs) - 1.96*(sd(boot_results_bw$cstat2_bs) / sqrt(1000))
test_cstat2_LCI_1 <- mean(boot_results_bw$cstat2_test) - 1.96*(sd(boot_results_bw$cstat2_test) / sqrt(1000))
LCI_1 <- cstat.2.LCI - (boot_cstat2_LCI_1 - test_cstat2_LCI_1)
LCI_1 #0.69
boot_cstat2_UCI_1 <- mean(boot_results_bw$cstat2_bs) + 1.96*(sd(boot_results_bw$cstat2_bs) / sqrt(1000))
test_cstat2_UCI_1 <- mean(boot_results_bw$cstat2_test) + 1.96*(sd(boot_results_bw$cstat2_test) / sqrt(1000))
UCI_1 <- cstat.2.UCI - (boot_cstat2_UCI_1 - test_cstat2_UCI_1)
UCI_1 #0.74

cstat_3vs1_adj <- cstat.model.3$auc[1] - (mean(boot_results_bw$cstat3_bs) - mean(boot_results_bw$cstat3_test))
cstat_3vs1_adj #0.53 (0.48, 0.59)

#manual CI
boot_cstat3_LCI_1 <- mean(boot_results_bw$cstat3_bs) - 1.96*(sd(boot_results_bw$cstat3_bs) / sqrt(1000))
test_cstat3_LCI_1 <- mean(boot_results_bw$cstat3_test) - 1.96*(sd(boot_results_bw$cstat3_test) / sqrt(1000))
LCI_1 <- cstat.3.LCI - (boot_cstat3_LCI_1 - test_cstat3_LCI_1)
LCI_1 #0.48

boot_cstat3_UCI_1 <- mean(boot_results_bw$cstat3_bs) + 1.96*(sd(boot_results_bw$cstat3_bs) / sqrt(1000))
test_cstat3_UCI_1 <- mean(boot_results_bw$cstat3_test) + 1.96*(sd(boot_results_bw$cstat3_test) / sqrt(1000))
UCI_1 <- cstat.3.UCI - (boot_cstat3_UCI_1 - test_cstat3_UCI_1)
UCI_1 #0.59

#### optimism on its own
#PDI
mean(boot_results_bw$boot_PDI) - mean(boot_results_bw$test_PDI)
#0.012

#c-slope y=2(vs1)
mean(boot_results_bw$boot_cslope1) - mean(boot_results_bw$test_cslope1)
#0.028

#c-slope y=3(vs1)
mean(boot_results_bw$boot_cslope3) - mean(boot_results_bw$test_cslope3)
#0.229

#c-intercept y=2(vs1)
mean(boot_results_bw$boot_intercept1) - mean(boot_results_bw$test_intercept1)
#-0.000

#c-intercept y=3(vs1)
mean(boot_results_bw$boot_intercept3) - mean(boot_results_bw$test_intercept3)
#-0.003

#c-statistic y=2(vs1)
mean(boot_results_bw$cstat2_bs) - mean(boot_results_bw$cstat2_test)
#0.004

#c-statistic y=3(vs1)
mean(boot_results_bw$cstat3_bs) - mean(boot_results_bw$cstat3_test)
#0.023

#c-statistic y=3(vs1)
mean(boot_results_bw$boot_R2) - mean(boot_results_bw$test_R2)
#0.009

# adjusting model coefficients for optimism
### Submodel 1
# multiply beta coeff by adjusted c-slope
pooled_sum_adj_2 <- pooled_sum %>% filter(outcome==2)
# multiplying beta coefficients by adjusted calibration slope (shrinkage)
pooled_sum_adj_2$adj_coeff <- pooled_sum_adj_2$estimate * Calslope.adj.1
# calculate OR by exponentiating the adjusted beta coefficients
pooled_sum_adj_2$adj_OR <- exp(pooled_sum_adj_2$adj_coeff)

### Submodel 2
# multiply beta coeff by adjusted c-slope
pooled_sum_adj_3 <- pooled_sum %>% filter(outcome==3)
pooled_sum_adj_3$adj_coeff <- pooled_sum_adj_3$estimate * Calslope.adj.3
# calculate OR by exponentiating the adjusted beta coefficients
pooled_sum_adj_3$adj_OR <- exp(pooled_sum_adj_3$adj_coeff)

### re-joining both dataframes
pooled_sum_adj <- rbind(pooled_sum_adj_2, pooled_sum_adj_3)