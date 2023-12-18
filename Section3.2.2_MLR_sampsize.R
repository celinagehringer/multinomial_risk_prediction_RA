#########################
# Section 3.2.1 of manuscript
# Sample size calculation for the development of a multinomial prediction model
# Code from Pate et al.: https://doi.org/10.1177/09622802231151220)
#########################

set.seed(101)

### Load relevant packages
library(pROC)

###################################################################################
##### STEP'S 1 and 2: Identify values for Q, p_k, p_k_r, max(R2_CS), R2_CS_adj and 
##### R2_CS_adj_k_r
###################################################################################

#First, define the number of events in each category
EV1 <- 756 #non-response
EV2 <- 730 #response
EV3 <- 146 #adverse events

EV1 <- 620 #non-response
EV2 <- 580 #response
EV3 <- 111 #adverse events

############
### Define Q - number of predictor parameters
############
Q <- 5

########################
### Define p_k and p_k_r
########################

## Define total number of events
n.total <- EV1 + EV2 + EV3 

## Calculate p_k
p.1 <- EV1/n.total
p.2 <- EV2/n.total
p.3 <- EV3/n.total

p.1
p.2
p.3


## Calculate p_k_r
p.1.2 <- (EV1 + EV2)/n.total
p.1.3 <- (EV1 + EV3)/n.total
p.2.3 <- (EV2 + EV3)/n.total


########################
### Calculate max(R2_CS)
########################
max_R2_CS <- 1 - ((p.1^p.1)*(p.2^p.2)*(p.3^p.3))^2
max_R2_CS


########################
### Calculate R2_CS_adj
########################

### Calculate an estimte of R2_CS_app, based off R2_NAGEL = 0.15
R2_CS_adj <- 0.15*max_R2_CS
R2_CS_adj

###########################
### Calculate R2_CS_adj_k_r
###########################

###Since no prior information on c-statistics is available, we will use the conservative approach of R2=0.15

## Calculate pairwise outcome proportions (phi), of category k relative to category i
phi.1.2 <- EV2/(EV1 + EV2)
phi.1.3 <- EV3/(EV1 + EV3)
phi.2.3 <- EV3/(EV2 + EV3)

phi.1.2
phi.1.3
phi.2.3


R2_app_kr.1.2 <- 0.15 * (1 - ((phi.1.2^(phi.1.2))*((1-phi.1.2)^(1-phi.1.2)))^2)
R2_app_kr.1.2

R2_app_kr.1.3 <- 0.15 * (1 - ((phi.1.3^(phi.1.3))*((1-phi.1.3)^(1-phi.1.3)))^2)
R2_app_kr.1.3

R2_app_kr.2.3 <- 0.15 * (1 - ((phi.2.3^(phi.2.3))*((1-phi.2.3)^(1-phi.2.3)))^2)
R2_app_kr.2.3


###########################
###########################
##### STEP 3: Criterion (i)
###########################
###########################

## Let S be the level of shrinkage we are targeting
S <- 0.9

## Calculate m_k_r
m.1.2 <- Q/((S - 1)*log(1 - R2_app_kr.1.2/S))
m.1.3 <- Q/((S - 1)*log(1 - R2_app_kr.1.3/S))
m.2.3 <- Q/((S - 1)*log(1 - R2_app_kr.2.3/S))


### Calculate n_k_r for criterion (i) for each submodel
N_C1.1.2 <- m.1.2/p.1.2
N_C1.1.3 <- m.1.3/p.1.3
N_C1.2.3 <- m.2.3/p.2.3


N_C1.1.2
N_C1.1.3
N_C1.2.3



### Take the ceiling of the maximum of these as the sample size for criteiron (i)
N_C1 <- ceiling(max(N_C1.1.2, N_C1.1.3, N_C1.2.3))
N_C1


### Now calculate number of each event we expect to see in a datast of this size
N_C1*p.1
N_C1*p.2
N_C1*p.3



############################
############################
##### STEP 4: Criterion (ii)
############################
############################

N_C2 <- 4*Q/((R2_CS_adj/(R2_CS_adj + 0.05*max_R2_CS) - 1)*log(1 - R2_CS_adj - 0.05*max_R2_CS))
N_C2 <- ceiling(N_C2)
N_C2

### Now calculate number of each event we expect to see in a datast of this size
N_C2*p.1
N_C2*p.2
N_C2*p.3



#############################
#############################
##### STEP 5: Criterion (iii)
#############################
#############################

N_C3.1 <- qchisq(1-0.05/5, 1)*p.1*(1-p.1)/0.05^2 
N_C3.2 <- qchisq(1-0.05/5, 1)*p.2*(1-p.2)/0.05^2 
N_C3.3 <- qchisq(1-0.05/5, 1)*p.3*(1-p.3)/0.05^2 

N_C3.1
N_C3.2
N_C3.3


N_C3 <- ceiling(max(N_C3.1, N_C3.2, N_C3.3))
N_C3

### Now calculate number of each event we expect to see in a dataset of this size
N_C3*p.1
N_C3*p.2
N_C3*p.3



#####################################################################
#####################################################################
##### STEP 6: Take the maximum sample size across all three criteria
#####################################################################
#####################################################################

N_C1
N_C2
N_C3

N <- max(N_C1, N_C2, N_C3)
N

### Minimum sample size: n=1430
