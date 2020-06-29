pkgname <- "riskRegression"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('riskRegression')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CSC")
### * CSC

flush(stderr()); flush(stdout())

### Name: CSC
### Title: Cause-specific Cox proportional hazard regression
### Aliases: CSC
### Keywords: survival

### ** Examples


library(prodlim)
library(survival)
data(Melanoma)
## fit two cause-specific Cox models
## different formula for the two causes
fit1 <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
            data=Melanoma)
print(fit1)
## Not run: 
##D library(Publish)
##D publish(fit1)
## End(Not run)

## model hazard of all cause mortality instead of hazard of type 2 
fit1a <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
             data=Melanoma,
             surv.type="surv")

## the predicted probabilities are similar 
plot(predictRisk(fit1,times=500,cause=1,newdata=Melanoma),
     predictRisk(fit1a,times=500,cause=1,newdata=Melanoma))

## special case where cause 2 has no covariates
fit1b <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~1),
             data=Melanoma)
print(fit1b)
predict(fit1b,cause=1,times=100,newdata=Melanoma)


## same formula for both causes
fit2 <- CSC(Hist(time,status)~invasion+epicel+age,
            data=Melanoma)
print(fit2)

## combine a cause-specific Cox regression model for cause 2 
## and a Cox regression model for the event-free survival:
## different formula for cause 2 and event-free survival
fit3 <- CSC(list(Hist(time,status)~sex+invasion+epicel+age,
                 Hist(time,status)~invasion+epicel+age),
            surv.type="surv",
            data=Melanoma)
print(fit3)

## same formula for both causes
fit4 <- CSC(Hist(time,status)~invasion+epicel+age,
            data=Melanoma,
            surv.type="surv")
print(fit4)

## strata
fit5 <- CSC(Hist(time,status)~invasion+epicel+age+strata(sex),
            data=Melanoma,
            surv.type="surv")
print(fit5)

## sanity checks

cox1 <- coxph(Surv(time,status==1)~invasion+epicel+age+strata(sex),data=Melanoma)
cox2 <- coxph(Surv(time,status!=0)~invasion+epicel+age+strata(sex),data=Melanoma)
all.equal(coef(cox1),coef(fit5$models[[1]]))
all.equal(coef(cox2),coef(fit5$models[[2]]))

## predictions
##
## surv.type = "hazard": predictions for both causes can be extracted
## from the same fit
fit2 <- CSC(Hist(time,status)~invasion+epicel+age, data=Melanoma)
predict(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
predictRisk(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
predictRisk(fit2,cause=2,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
predict(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
predict(fit2,cause=2,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))

## surv.type = "surv" we need to change the cause of interest 
library(survival)
fit5.2 <- CSC(Hist(time,status)~invasion+epicel+age+strata(sex),
            data=Melanoma,
            surv.type="surv",cause=2)
## now this does not work
try(predictRisk(fit5.2,cause=1,newdata=Melanoma,times=4))

## but this does
predictRisk(fit5.2,cause=2,newdata=Melanoma,times=100)
predict(fit5.2,cause=2,newdata=Melanoma,times=100)
predict(fit5.2,cause=2,newdata=Melanoma[4,],times=100)




cleanEx()
nameEx("Ctree")
### * Ctree

flush(stderr()); flush(stdout())

### Name: Ctree
### Title: S3-Wrapper for ctree.
### Aliases: Ctree

### ** Examples

library(prodlim)
library(party)
library(survival)
set.seed(50)
d <- SimSurv(50)
nd <- data.frame(X1=c(0,1,0),X2=c(-1,0,1))
f <- Ctree(Surv(time,status)~X1+X2,data=d)
predictRisk(f,newdata=nd,times=c(3,8))




cleanEx()
nameEx("FGR")
### * FGR

flush(stderr()); flush(stdout())

### Name: FGR
### Title: Formula wrapper for crr from cmprsk
### Aliases: FGR
### Keywords: survival

### ** Examples


library(prodlim)
library(survival)
library(cmprsk)
library(lava)
d <- prodlim::SimCompRisk(100)
f1 <- FGR(Hist(time,cause)~X1+X2,data=d)
print(f1)

## crr allows that some covariates are multiplied by
## a function of time (see argument tf of crr)
## by FGR uses the identity matrix
f2 <- FGR(Hist(time,cause)~cov2(X1)+X2,data=d)
print(f2)

## same thing, but more explicit:
f3 <- FGR(Hist(time,cause)~cov2(X1)+cov1(X2),data=d)
print(f3)

## both variables can enter cov2:
f4 <- FGR(Hist(time,cause)~cov2(X1)+cov2(X2),data=d)
print(f4)

## change the function of time
qFun <- function(x){x^2}
noFun <- function(x){x}
sqFun <- function(x){x^0.5}

## multiply X1 by time^2 and X2 by time:
f5 <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2),data=d)
print(f5)
print(f5$crrFit)
## same results as crr
with(d,crr(ftime=time,
           fstatus=cause,
           cov2=d[,c("X1","X2")],
           tf=function(time){cbind(qFun(time),time)}))

## still same result, but more explicit
f5a <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2,tf=noFun),data=d)
f5a$crrFit

## multiply X1 by time^2 and X2 by sqrt(time)
f5b <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2,tf=sqFun),data=d,cause=1)

## additional arguments for crr
f6<- FGR(Hist(time,cause)~X1+X2,data=d, cause=1,gtol=1e-5)
f6
f6a<- FGR(Hist(time,cause)~X1+X2,data=d, cause=1,gtol=0.1)
f6a



cleanEx()
nameEx("IPA")
### * IPA

flush(stderr()); flush(stdout())

### Name: IPA
### Title: Explained variation for settings with binary, survival and
###   competing risk outcome
### Aliases: IPA rsquared rsquared.default rsquared.glm rsquared.coxph
###   rsquared.CauseSpecificCox IPA.default IPA.glm IPA.coxph
###   IPA.CauseSpecificCox

### ** Examples

library(prodlim)
library(data.table)
# binary outcome
library(lava)
set.seed(18)
learndat <- sampleData(48,outcome="binary")
lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
IPA(lr1)

## validation data
valdat=sampleData(94,outcome="binary")
IPA(lr1,newdata=valdat)

## predicted risks externally given
p1=predictRisk(lr1,newdata=valdat)
IPA(p1,formula=Y~1,valdat)

# survival
library(survival)
data(pbc)
pbc=na.omit(pbc)
pbctest=(1:NROW(pbc)) %in% sample(1:NROW(pbc),size=.632*NROW(pbc))
pbclearn=pbc[pbctest,]
cox1= coxph(Surv(time,status!=0)~age+sex+log(bili)+log(albumin)+log(protime),
      data=pbclearn,x=TRUE)

## same data
IPA(cox1,formula=Surv(time,status!=0)~1,times=1000)

## validation data
pbcval=pbc[!pbctest,]
IPA(cox1,formula=Surv(time,status!=0)~1,newdata=pbcval,times=1000)

## predicted risks externally given
p2=predictRisk(cox1,newdata=pbcval,times=1000)
IPA(cox1,formula=Surv(time,status!=0)~1,newdata=pbcval,times=1000)
 
# competing risks
data(Melanoma)
Melanomatest=(1:NROW(Melanoma)) %in% sample(1:NROW(Melanoma),size=.632*NROW(Melanoma))
Melanomalearn=Melanoma[Melanomatest,]
fit1 <- CSC(list(Hist(time,status)~sex,
                 Hist(time,status)~invasion+epicel+age),
                 data=Melanoma)
IPA(fit1,times=1000,cause=2)

## validation data
Melanomaval=Melanoma[!Melanomatest,]
IPA(fit1,formula=Hist(time,status)~1,newdata=Melanomaval,times=1000)

## predicted risks externally given
p3= predictRisk(fit1,cause=1,newdata=Melanomaval,times=1000)
IPA(p3,formula=Hist(time,status)~1,cause=1,newdata=Melanomaval,times=1000)
 



cleanEx()
nameEx("Melanoma")
### * Melanoma

flush(stderr()); flush(stdout())

### Name: Melanoma
### Title: Malignant melanoma data
### Aliases: Melanoma
### Keywords: datasets

### ** Examples


data(Melanoma)



cleanEx()
nameEx("Paquid")
### * Paquid

flush(stderr()); flush(stdout())

### Name: Paquid
### Title: Paquid sample
### Aliases: Paquid
### Keywords: datasets

### ** Examples

data(Paquid)



cleanEx()
nameEx("Score.list")
### * Score.list

flush(stderr()); flush(stdout())

### Name: Score.list
### Title: Score risk predictions
### Aliases: Score.list Score

### ** Examples

# binary outcome
library(lava)
set.seed(18)
learndat <- sampleData(48,outcome="binary")
testdat <- sampleData(40,outcome="binary")

## score logistic regression models
lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
lr2 = glm(Y~X3+X5,data=learndat,family=binomial)
Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5)"=lr2),formula=Y~1,data=testdat)

## ROC curve and calibration plot
xb=Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5+X6)"=lr2),formula=Y~1,
         data=testdat,plots=c("calibration","ROC"))
## Not run: 
##D plotROC(xb)
##D plotCalibration(xb)
## End(Not run)

## compute AUC for a list of continuous markers
markers = as.list(testdat[,.(X6,X7,X8,X9,X10)])
Score(markers,formula=Y~1,data=testdat,metrics=c("auc"))

# cross-validation
## Not run: 
##D     learndat=sampleData(400,outcome="binary")
##D     lr1a = glm(Y~X6,data=learndat,family=binomial)
##D     lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
##D     ## bootstrap cross-validation
##D     x1=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=100)
##D     x1
##D     ## leave-one-out and leave-pair-out bootstrap
##D     x2=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,
##D              split.method="loob",
##D              B=100,plots="calibration")
##D     x2
## End(Not run)
# survival outcome

# Score Cox regression models
## Not run: 
##D library(survival)
##D library(rms)
##D library(prodlim)
##D set.seed(18)
##D trainSurv <- sampleData(100,outcome="survival")
##D testSurv <- sampleData(40,outcome="survival")
##D cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
##D cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
##D xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##D          formula=Surv(time,event)~1,data=testSurv,conf.int=FALSE,times=c(5,8))
##D xs
## End(Not run)

# Integrated Brier score
## Not run: 
##D xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##D          formula=Surv(time,event)~1,data=testSurv,conf.int=FALSE,
##D          summary="ibs",
##D          times=sort(unique(testSurv$time)))
## End(Not run)

# time-dependent AUC for list of markers
## Not run: 
##D survmarkers = as.list(testSurv[,.(X6,X7,X8,X9,X10)])
##D Score(survmarkers,
##D       formula=Surv(time,event)~1,metrics="auc",data=testSurv,
##D       conf.int=TRUE,times=c(5,8))
##D 
##D # compare models on test data
##D Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##D       formula=Surv(time,event)~1,data=testSurv,conf.int=TRUE,times=c(5,8))
## End(Not run)
# crossvalidation models in traindata
## Not run: 
##D     library(survival)
##D     set.seed(18)
##D     trainSurv <- sampleData(400,outcome="survival")
##D     cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
##D     cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
##D     x1 = Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##D                formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=c(5,8),
##D                split.method="loob",B=100,plots="calibration")
##D 
##D     x2= Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##D               formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=c(5,8),
##D               split.method="bootcv",B=100)
## End(Not run)

# restrict number of comparisons
## Not run: 
##D     Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##D           formula=Surv(time,event)~1,data=trainSurv,contrasts=TRUE,
##D           null.model=FALSE,conf.int=TRUE,times=c(5,8),split.method="bootcv",B=3)
##D 
##D     # competing risks outcome
##D     set.seed(18)
##D     trainCR <- sampleData(40,outcome="competing.risks")
##D     testCR <- sampleData(40,outcome="competing.risks")
##D     library(riskRegression)
##D     library(cmprsk)
##D     # Cause-specific Cox regression
##D     csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR)
##D     csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR)
##D     # Fine-Gray regression
##D     fgr1 = FGR(Hist(time,event)~X1+X2+X7+X9,data=trainCR,cause=1)
##D     fgr2 = FGR(Hist(time,event)~X3+X5+X6,data=trainCR,cause=1)
##D     Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2,
##D                "FGR(X1+X2+X7+X9)"=fgr1,"FGR(X3+X5+X6)"=fgr2),
##D           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(5,8))
## End(Not run)



## Not run: 
##D     # reproduce some results of Table IV of Blanche et al. Stat Med 2013
##D     data(Paquid)
##D     ResPaquid <- Score(list("DSST"=-Paquid$DSST,"MMSE"=-Paquid$MMSE),
##D                        formula=Hist(time,status)~1,
##D                        data=Paquid,
##D                        null.model = FALSE,
##D                        conf.int=TRUE,
##D                        metrics=c("auc"),
##D                        times=c(3,5,10),
##D                        plots="ROC")
##D     ResPaquid
##D     plotROC(ResPaquid,time=5)
## End(Not run)



cleanEx()
nameEx("SuperPredictor")
### * SuperPredictor

flush(stderr()); flush(stdout())

### Name: SuperPredictor
### Title: Formula interface for SuperLearner::SuperLearner
### Aliases: SuperPredictor

### ** Examples

## Not run: 
##D library(SuperLearner)
##D library(data.table)
##D d = sampleData(338, outcome="binary")
##D spfit = SuperPredictor(Y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=d)
##D predictRisk(spfit)
##D x <- Score(list(spfit),data=d,formula=Y~1)
## End(Not run)



cleanEx()
nameEx("SurvResponseVar")
### * SurvResponseVar

flush(stderr()); flush(stdout())

### Name: SurvResponseVar
### Title: Extract the time and event variable from a Cox model
### Aliases: SurvResponseVar

### ** Examples

## Not run: 
##D SurvResponseVar(Surv(time,event)~X1+X2)
##D SurvResponseVar(Hist(time,event==0)~X1+X2)
##D SurvResponseVar(Surv(start,time, status,type="counting") ~ X3+X5)
##D SurvResponseVar(Surv(start,event=status, time2=time,type="counting") ~ X3+X5)
##D 
##D SurvResponseVar(survival::Surv(start,event=status, time2=time,type="counting") ~ X3+X5)
##D SurvResponseVar(status ~ X3+X5)
##D SurvResponseVar(I(status == 1) ~ X3+X5)
##D SurvResponseVar(list(Hist(time, event) ~ X1+X6,Hist(time, event) ~ X6))
## End(Not run)



cleanEx()
nameEx("ate")
### * ate

flush(stderr()); flush(stdout())

### Name: ate
### Title: Compute the Average Treatment Effects Via
### Aliases: ate

### ** Examples

library(survival)
library(rms)
library(prodlim)
set.seed(10)

#### Survival settings  ####
#### ATE with Cox model ####

## generate data
n <- 100
dtS <- sampleData(n, outcome="survival")
dtS$time <- round(dtS$time,1)
dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))

## estimate the Cox model
fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)

## compute the ATE at times 5, 6, 7, and 8 using X1 as the treatment variable
## Not run: 
##D ## only point estimate (argument se = FALSE)
##D ateFit1a <- ate(fit, data = dtS, treatment = "X1", times = 5:8,
##D                se = FALSE)
##D 
##D ## standard error / confidence intervals computed using the influence function
##D ## (argument se = TRUE and B = 0)
##D ateFit1b <- ate(fit, data = dtS, treatment = "X1", times = 5:8,
##D                se = TRUE, B = 0)
##D 
##D ## same as before with in addition the confidence bands for the ATE
##D ## (argument band = TRUE)
##D ateFit1c <- ate(fit, data = dtS, treatment = "X1", times = 5:8,
##D                se = TRUE, band = TRUE, B = 0)
##D 
##D ## standard error / confidence intervals computed using 100 boostrap samples
##D ## (argument se = TRUE and B = 100) 
##D ateFit1d <- ate(fit, data = dtS, treatment = "X1",
##D                 times = 5:8, se = TRUE, B = 100)
##D ## NOTE: for real applications 100 bootstrap samples is not enougth 
##D 
##D ## same but using 2 cpus for generating and analyzing the boostrap samples
##D ## (parallel computation, argument mc.cores = 2) 
##D ateFit1e <- ate(fit, data = dtS, treatment = "X1",
##D                 times = 5:8, se = TRUE, B = 100, mc.cores = 2)
## End(Not run)

#### Survival settings without censoring ####
#### ATE with glm                        ####

## generate data
n <- 100
dtB <- sampleData(n, outcome="binary")
dtB[, X2 := as.numeric(X2)]

## estimate a logistic regression model
fit <- glm(formula = Y ~ X1+X2, data=dtB, family = "binomial")

## compute the ATE using X1 as the treatment variable
## only point estimate (argument se = FALSE)
ateFit1a <- ate(fit, data = dtB, treatment = "X1", se = FALSE)

## Not run: 
##D ## standard error / confidence intervals computed using the influence function
##D ateFit1b <- ate(fit, data = dtB, treatment = "X1",
##D                times = 5, ## just for having a nice output not used in computations
##D                se = TRUE, B = 0)
##D 
##D ## standard error / confidence intervals computed using 100 boostrap samples
##D ateFit1d <- ate(fit, data = dtB, treatment = "X1",
##D                 times = 5, se = TRUE, B = 100)
##D 
##D ## using the lava package
##D ateLava <- estimate(fit, function(p, data){
##D a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X2"] ;
##D R.X11 <- expit(a + b + c * data[["X2"]])
##D R.X10 <- expit(a + c * data[["X2"]])
##D list(risk0=R.X10,risk1=R.X11,riskdiff=R.X11-R.X10)},
##D average=TRUE)
##D ateLava
##D 
##D ateFit1b$meanRisk
## End(Not run)

#### Competing risks settings               ####
#### ATE with cause specific Cox regression ####

## Not run: 
##D ## generate data
##D n <- 500
##D set.seed(10)
##D dt <- sampleData(n, outcome="competing.risks")
##D dt$time <- round(dt$time,1)
##D dt$X1 <- factor(rbinom(n, prob = c(0.2,0.3) , size = 2), labels = paste0("T",0:2))
##D 
##D ## estimate cause specific Cox model
##D fitCR <-  CSC(Hist(time,event)~ X1+X8,data=dt,cause=1)
##D 
##D ## compute the ATE at times 10, 15, 20 using X1 as the treatment variable
##D ateFit2a <- ate(fitCR, data = dt, treatment = "X1", times = c(10,15,20),
##D                 cause = 1, se = FALSE)
##D 
##D ## standard error / confidence intervals computed using the influence function
##D ## (argument se = TRUE and B = 0)
##D ateFit2b <- ate(fitCR, data = dt, treatment = "X1", times = c(10,15,20),
##D                 cause = 1, se = TRUE, B = 0)
##D 
##D ## same as before with in addition the confidence bands for the ATE
##D ## (argument band = TRUE)
##D ateFit2c <- ate(fitCR, data = dt, treatment = "X1", times = c(10,15,20), 
##D                cause = 1, se = TRUE, band = TRUE, B = 0)
##D 
##D ## standard error / confidence intervals computed using 100 boostrap samples
##D ## (argument se = TRUE and B = 100) 
##D ateFit2d <- ate(fitCR, data = dt, treatment = "X1", times = c(10,15,20), 
##D                 cause = 1, se = TRUE, B = 100)
##D ## NOTE: for real applications 100 bootstrap samples is not enougth 
##D 
##D ## same but using 2 cpus for generating and analyzing the boostrap samples
##D ## (parallel computation, argument mc.cores = 2) 
##D ateFit2e <- ate(fitCR, data = dt, treatment = "X1", times = c(10,15,20), 
##D                 cause = 1, se = TRUE, B = 100, mc.cores = 2)
## End(Not run)

#### time-dependent covariates ###
## Not run: 
##D library(survival)
##D fit <- coxph(Surv(time, status) ~ celltype+karno + age + trt, veteran)
##D vet2 <- survSplit(Surv(time, status) ~., veteran,
##D                        cut=c(60, 120), episode ="timegroup")
##D fitTD <- coxph(Surv(tstart, time, status) ~ celltype+karno + age + trt,
##D                data= vet2,x=1)
##D set.seed(16)
##D resVet <- ate(fitTD,formula=Hist(entry=tstart,time=time,event=status)~1,
##D           data = vet2, treatment = "celltype", contrasts = NULL,
##D         times=5,verbose=1,
##D         landmark = c(0,30,60,90), cause = 1, B = 10, se = 1,
##D         band = FALSE, mc.cores=1)
##D resVet
## End(Not run)

## Not run: 
##D set.seed(137)
##D d=sampleDataTD(127)
##D library(survival)
##D d[,status:=1*(event==1)]
##D d[,X3:=as.factor(X3)]
##D ## ignore competing risks
##D cox1TD <- coxph(Surv(start,time, status,type="counting") ~ X3+X5+X6+X8,
##D                 data=d, x = TRUE)
##D resTD1 <- ate(cox1TD,formula=Hist(entry=start,time=time,event=status)~1,
##D         data = d, treatment = "X3", contrasts = NULL,
##D         times=.5,verbose=1,
##D         landmark = c(0,0.5,1), B = 20, se = 1,
##D         band = FALSE, mc.cores=1)
##D resTD1
##D ## account for competing risks
##D cscTD <- CSC(Hist(time=time, event=event,entry=start) ~ X3+X5+X6+X8, data=d)
##D set.seed(16)
##D resTD <- ate(cscTD,formula=Hist(entry=start,time=time,event=event)~1,
##D         data = d, treatment = "X3", contrasts = NULL,
##D         times=.5,verbose=1,
##D         landmark = c(0,0.5,1), cause = 1, B = 20, se = 1,
##D         band = FALSE, mc.cores=1)
##D resTD
## End(Not run)



cleanEx()
nameEx("autoplot.Score")
### * autoplot.Score

flush(stderr()); flush(stdout())

### Name: autoplot.Score
### Title: ggplot AUC curve
### Aliases: autoplot.Score

### ** Examples

library(survival)
library(ggplot2)
d=sampleData(100,outcome="survival")
nd=sampleData(100,outcome="survival")
f1=coxph(Surv(time,event)~X1+X6+X8,data=d,x=TRUE,y=TRUE)
f2=coxph(Surv(time,event)~X2+X5+X9,data=d,x=TRUE,y=TRUE)
xx=Score(list(f1,f2), formula=Surv(time,event)~1,
data=nd, metrics="auc", null.model=FALSE, times=seq(3:10))
g <- autoplot(xx)
print(g)
aucgraph <- plotAUC(xx)
plotAUC(xx,conf.int=TRUE)
plotAUC(xx,which="contrasts")
plotAUC(xx,which="contrasts",conf.int=TRUE)





cleanEx()
nameEx("autoplot.ate")
### * autoplot.ate

flush(stderr()); flush(stdout())

### Name: autoplot.ate
### Title: Plot Average Risks
### Aliases: autoplot.ate

### ** Examples

## Not run: 
##D library(survival)
##D library(rms)
##D library(ggplot2)
##D #### simulate data ####
##D n <- 1e2
##D set.seed(10)
##D dtS <- sampleData(n,outcome="survival")
##D 
##D 
##D #### Cox model ####
##D fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
##D 
##D #### Average treatment effect ####
##D seqTimes <- sort(unique(fit$y[,1]))
##D seqTimes5 <- seqTimes[seqTimes>5 & seqTimes<10]
##D ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
##D               times = seqTimes, B = 0, band = TRUE, nsim.band = 500, y = TRUE,
##D               mc.cores=1)
##D 
##D #### display #### 
##D ggplot2::autoplot(ateFit)
##D 
##D outGG <- autoplot(ateFit, band = TRUE, ci = TRUE, alpha = 0.1)
##D dd <- as.data.frame(outGG$data[treatment == 0])
##D outGG$plot + facet_wrap(~treatment, labeller = label_both)
## End(Not run)



cleanEx()
nameEx("autoplot.predictCSC")
### * autoplot.predictCSC

flush(stderr()); flush(stdout())

### Name: autoplot.predictCSC
### Title: Plot Predictions From a Cause-specific Cox Proportional Hazard
###   Regression
### Aliases: autoplot.predictCSC

### ** Examples

library(survival)
library(rms)
library(ggplot2)
library(prodlim)
#### simulate data ####
set.seed(10)
d <- sampleData(1e2, outcome = "competing.risks")

#### CSC model ####
m.CSC <- CSC(Hist(time,event)~ X1 + X2 + X6, data = d)

pred.CSC <- predict(m.CSC, newdata = d[1:2,], time = 1:5, cause = 1)#'
autoplot(pred.CSC)


#### stratified CSC model ####
m.SCSC <- CSC(Hist(time,event)~ strata(X1) + strata(X2) + X6,
              data = d)
pred.SCSC <- predict(m.SCSC, time = 1:3, newdata = d[1:4,],
                     cause = 1, keep.newdata = TRUE, keep.strata = TRUE)
autoplot(pred.SCSC, group.by = "strata")



cleanEx()
nameEx("autoplot.predictCox")
### * autoplot.predictCox

flush(stderr()); flush(stdout())

### Name: autoplot.predictCox
### Title: Plot Predictions From a Cox Model
### Aliases: autoplot.predictCox

### ** Examples

library(survival)
library(ggplot2)

#### simulate data ####
set.seed(10)
d <- sampleData(1e2, outcome = "survival")

#### Cox model ####
m.cox <- coxph(Surv(time,event)~ X1 + X2 + X3,
                data = d, x = TRUE, y = TRUE)

## display baseline hazard
e.basehaz <- predictCox(m.cox)

autoplot(e.basehaz, type = "cumhazard")

## display predicted survival
pred.cox <- predictCox(m.cox, newdata = d[1:4,],
  times = 1:5, type = "survival", keep.newdata = TRUE)
autoplot(pred.cox)
autoplot(pred.cox, group.by = "covariates")
autoplot(pred.cox, group.by = "covariates", reduce.data = TRUE)

## predictions with confidence interval/bands
pred.cox <- predictCox(m.cox, newdata = d[1,,drop=FALSE],
  times = 1:5, type = "survival", band = TRUE, se = TRUE, keep.newdata = TRUE)
autoplot(pred.cox, ci = TRUE, band = TRUE)
autoplot(pred.cox, ci = TRUE, band = TRUE, alpha = 0.1)

#### Stratified Cox model ####
m.cox.strata <- coxph(Surv(time,event)~ strata(X1) + strata(X2) + X3 + X6,
                      data = d, x = TRUE, y = TRUE)

pred.cox.strata <- predictCox(m.cox.strata, newdata = d[1:5,,drop=FALSE],
                              time = 1:5, keep.newdata = TRUE)

## display
res <- autoplot(pred.cox.strata, type = "survival", group.by = "strata")

## customize display
res$plot + facet_wrap(~strata, labeller = label_both)
res$plot %+% res$data[strata == "0, 1"]



cleanEx()
nameEx("boot2pvalue")
### * boot2pvalue

flush(stderr()); flush(stdout())

### Name: boot2pvalue
### Title: Compute the p.value from the distribution under H1
### Aliases: boot2pvalue

### ** Examples

set.seed(10)

#### no effect ####
x <- rnorm(1e3) 
boot2pvalue(x, null = 0, estimate = mean(x), alternative = "two.sided")
## expected value of 1
boot2pvalue(x, null = 0, estimate = mean(x), alternative = "greater")
## expected value of 0.5
boot2pvalue(x, null = 0, estimate = mean(x), alternative = "less")
## expected value of 0.5

#### positive effect ####
x <- rnorm(1e3, mean = 1) 
boot2pvalue(x, null = 0, estimate = 1, alternative = "two.sided")
## expected value of 0.32 = 2*pnorm(q = 0, mean = -1) = 2*mean(x<=0)
boot2pvalue(x, null = 0, estimate = 1, alternative = "greater")  
## expected value of 0.16 = pnorm(q = 0, mean = 1) = mean(x<=0)
boot2pvalue(x, null = 0, estimate = 1, alternative = "less")
## expected value of 0.84 = 1-pnorm(q = 0, mean = 1) = mean(x>=0)

#### negative effect ####
x <- rnorm(1e3, mean = -1) 
boot2pvalue(x, null = 0, estimate = -1, alternative = "two.sided") 
## expected value of 0.32 = 2*(1-pnorm(q = 0, mean = -1)) = 2*mean(x>=0)
boot2pvalue(x, null = 0, estimate = -1, alternative = "greater")
## expected value of 0.84 = pnorm(q = 0, mean = -1) = mean(x<=0)
boot2pvalue(x, null = 0, estimate = -1, alternative = "less") # pnorm(q = 0, mean = -1)
## expected value of 0.16 = 1-pnorm(q = 0, mean = -1) = mean(x>=0)



cleanEx()
nameEx("boxplot.Score")
### * boxplot.Score

flush(stderr()); flush(stdout())

### Name: boxplot.Score
### Title: Boxplot risk quantiles
### Aliases: boxplot.Score

### ** Examples

# binary outcome
library(data.table)
library(prodlim)
db=sampleData(40,outcome="binary")
fitconv=glm(Y~X3+X5,data=db,family=binomial)
fitnew=glm(Y~X1+X3+X5+X6+X7,data=db,family=binomial)
x=Score(list(new=fitnew,conv=fitconv),
        formula=Y~1,contrasts=list(c(2,1)),
               data=db,plots="box",null.model=FALSE)
boxplot(x)

# survival outcome
library(survival)
ds=sampleData(40,outcome="survival")
fit=coxph(Surv(time,event)~X6+X9,data=ds,x=TRUE,y=TRUE)
## Not run: 
##D  
##D scoreobj=Score(list("Cox"=fit),
##D                 formula=Hist(time,event)~1, data=ds,
##D                 metrics=NULL, plots="box",
##D                 times=c(1,5),null.model=FALSE)
##D boxplot(scoreobj,timepoint=5)
##D boxplot(scoreobj,timepoint=1)
##D 
## End(Not run)

# competing risks outcome
library(survival)
data(Melanoma, package = "riskRegression")
fit = CSC(Hist(time,event,cens.code="censored")~invasion+age+sex,data=Melanoma)
scoreobj=Score(list("CSC"=fit),
               formula=Hist(time,event,cens.code="censored")~1,
               data=Melanoma,plots="box",times=5*365.25,null.model=FALSE)
par(mar=c(4,12,4,4))
boxplot(scoreobj,timepoint=5*365.25)

# more than 2 competing risks
m=lava::lvm(~X1+X2+X3)
lava::distribution(m, "eventtime1") <- lava::coxWeibull.lvm(scale = 1/100)
lava::distribution(m, "eventtime2") <- lava::coxWeibull.lvm(scale = 1/100)
lava::distribution(m, "eventtime3") <- lava::coxWeibull.lvm(scale = 1/100)
lava::distribution(m, "censtime") <- lava::coxWeibull.lvm(scale = 1/100)
lava::regression(m,eventtime2~X3)=1.3
m <- lava::eventTime(m,
time ~ min(eventtime1 = 1, eventtime2 = 2, eventtime3 = 3, censtime = 0), "event")
set.seed(101)
dcr=as.data.table(lava::sim(m,101))
fit = CSC(Hist(time,event)~X1+X2+X3,data=dcr)
scoreobj=Score(list("my model"=fit),
               formula=Hist(time,event)~1,
               data=dcr,plots="box",times=5,null.model=FALSE)
boxplot(scoreobj)





graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("colCenter_cpp")
### * colCenter_cpp

flush(stderr()); flush(stdout())

### Name: colCenter_cpp
### Title: Apply - by column
### Aliases: colCenter_cpp

### ** Examples

x <- matrix(1,6,5)
sweep(x, MARGIN = 1, FUN = "-", STATS = 1:6)
colCenter_cpp(x, 1:6 )



cleanEx()
nameEx("colCumProd")
### * colCumProd

flush(stderr()); flush(stdout())

### Name: colCumProd
### Title: Apply cumprod in each column
### Aliases: colCumProd

### ** Examples

x <- matrix(1:8,ncol=2)
colCumProd(x)



cleanEx()
nameEx("colCumSum")
### * colCumSum

flush(stderr()); flush(stdout())

### Name: colCumSum
### Title: Apply cumsum in each column
### Aliases: colCumSum

### ** Examples

x <- matrix(1:8,ncol=2)
colCumSum(x)



cleanEx()
nameEx("colMultiply_cpp")
### * colMultiply_cpp

flush(stderr()); flush(stdout())

### Name: colMultiply_cpp
### Title: Apply * by column
### Aliases: colMultiply_cpp

### ** Examples

x <- matrix(1,6,5)
sweep(x, MARGIN = 1, FUN = "*", STATS = 1:6)
colMultiply_cpp(x, 1:6 )



cleanEx()
nameEx("colScale_cpp")
### * colScale_cpp

flush(stderr()); flush(stdout())

### Name: colScale_cpp
### Title: Apply / by column
### Aliases: colScale_cpp

### ** Examples

x <- matrix(1,6,5)
sweep(x, MARGIN = 1, FUN = "/", STATS = 1:6)
colScale_cpp(x, 1:6 )



cleanEx()
nameEx("colSumsCrossprod")
### * colSumsCrossprod

flush(stderr()); flush(stdout())

### Name: colSumsCrossprod
### Title: Apply crossprod and colSums
### Aliases: colSumsCrossprod

### ** Examples

x <- matrix(1:8,ncol=2)
y <- matrix(1:16,ncol=8)
colSumsCrossprod(x,y,0)

x <- matrix(1:8,ncol=2)
y <- matrix(1:16,ncol=2)
colSumsCrossprod(x,y,1)



cleanEx()
nameEx("confint.ate")
### * confint.ate

flush(stderr()); flush(stdout())

### Name: confint.ate
### Title: Confidence Intervals and Confidence Bands for the Predicted
###   Absolute Risk (Cumulative Incidence Function)
### Aliases: confint.ate

### ** Examples

library(survival)
library(data.table)

## ## generate data ####
set.seed(10)
d <- sampleData(70,outcome="survival")
d[, X1 := paste0("T",rbinom(.N, size = 2, prob = c(0.51)))]
## table(d$X1)

#### stratified Cox model ####
fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
             data=d, ties="breslow", x = TRUE, y = TRUE)

#### average treatment effect ####
fit.ate <- ate(fit, treatment = "X1", times = 1:3, data = d,
               se = TRUE, iid = TRUE, band = TRUE)
print(fit.ate, type = "meanRisk")
dt.ate <- as.data.table(fit.ate)

## manual calculation of se
dd <- copy(d)
dd$X1 <- rep(factor("T0", levels = paste0("T",0:2)), NROW(dd))
out <- predictCox(fit, newdata = dd, se = TRUE, times = 1:3, average.iid = TRUE)
term1 <- -out$survival.average.iid
term2 <- sweep(1-out$survival, MARGIN = 2, FUN = "-", STATS = colMeans(1-out$survival))
sqrt(colSums((term1 + term2/NROW(d))^2)) 
## fit.ate$meanRisk[treatment=="T0",meanRisk.se]

## note
out2 <- predictCox(fit, newdata = dd, se = TRUE, times = 1:3, iid = TRUE)
mean(out2$survival.iid[,1,1])
out$survival.average.iid[1,1]

## check confidence intervals (no transformation)
dt.ate[,.(lower = pmax(0,value + qnorm(0.025) * se),
          lower2 = lower,
          upper = value + qnorm(0.975) * se,
          upper2 = upper)]

## add confidence intervals computed on the log-log scale
## and backtransformed
outCI <- confint(fit.ate,
                 meanRisk.transform = "loglog", diffRisk.transform = "atanh",
                 ratioRisk.transform = "log")
print(outCI, type = "meanRisk")

dt.ate[type == "ate", newse := se/(value*log(value))]
dt.ate[type == "ate", .(lower = exp(-exp(log(-log(value)) - 1.96 * newse)),
                        upper = exp(-exp(log(-log(value)) + 1.96 * newse)))]



cleanEx()
nameEx("confint.predictCSC")
### * confint.predictCSC

flush(stderr()); flush(stdout())

### Name: confint.predictCSC
### Title: Confidence Intervals and Confidence Bands for the Predicted
###   Absolute Risk (Cumulative Incidence Function)
### Aliases: confint.predictCSC

### ** Examples

library(survival)
library(prodlim)
#### generate data ####
set.seed(10)
d <- sampleData(100) 

#### estimate a stratified CSC model ###
fit <- CSC(Hist(time,event)~ X1 + strata(X2) + X6, data=d)

#### compute individual specific risks
fit.pred <- predict(fit, newdata=d[1:3], times=c(3,8), cause = 1,
                    se = TRUE, iid = TRUE, band = TRUE)
fit.pred

## check confidence intervals
newse <- fit.pred$absRisk.se/(-fit.pred$absRisk*log(fit.pred$absRisk))
cbind(lower = as.double(exp(-exp(log(-log(fit.pred$absRisk)) + 1.96 * newse))),
      upper = as.double(exp(-exp(log(-log(fit.pred$absRisk)) - 1.96 * newse)))
)

#### compute confidence intervals without transformation
confint(fit.pred, absRisk.transform = "none")
cbind(lower = as.double(fit.pred$absRisk - 1.96 * fit.pred$absRisk.se),
      upper = as.double(fit.pred$absRisk + 1.96 * fit.pred$absRisk.se)
)





cleanEx()
nameEx("confint.predictCox")
### * confint.predictCox

flush(stderr()); flush(stdout())

### Name: confint.predictCox
### Title: Confidence Intervals and Confidence Bands for the predicted
###   Survival/Cumulative Hazard
### Aliases: confint.predictCox

### ** Examples

library(survival)

#### generate data ####
set.seed(10)
d <- sampleData(40,outcome="survival") 

#### estimate a stratified Cox model ####
fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
             data=d, ties="breslow", x = TRUE, y = TRUE)

#### compute individual specific survival probabilities  
fit.pred <- predictCox(fit, newdata=d[1:3], times=c(3,8), type = "survival",
                       se = TRUE, iid = TRUE, band = TRUE)
fit.pred

## check standard error
sqrt(rowSums(fit.pred$survival.iid[1,,]^2)) ## se for individual 1

## check confidence interval
newse <- fit.pred$survival.se/(-fit.pred$survival*log(fit.pred$survival))
cbind(lower = as.double(exp(-exp(log(-log(fit.pred$survival)) + 1.96 * newse))),
      upper = as.double(exp(-exp(log(-log(fit.pred$survival)) - 1.96 * newse)))
)

#### compute confidence intervals without transformation
confint(fit.pred, survival.transform = "none")
cbind(lower = as.double(fit.pred$survival - 1.96 * fit.pred$survival.se),
      upper = as.double(fit.pred$survival + 1.96 * fit.pred$survival.se)
)




cleanEx()
nameEx("getSplitMethod")
### * getSplitMethod

flush(stderr()); flush(stdout())

### Name: getSplitMethod
### Title: Input for data splitting algorithms
### Aliases: getSplitMethod

### ** Examples

# 3-fold crossvalidation
getSplitMethod("cv3",B=4,N=37)

# bootstrap with replacement
getSplitMethod("loob",B=4,N=37)

# bootstrap without replacement
getSplitMethod("loob",B=4,N=37,M=20)




cleanEx()
nameEx("iidCox")
### * iidCox

flush(stderr()); flush(stdout())

### Name: iidCox
### Title: Extract iid decomposition from a Cox model
### Aliases: iidCox iidCox.coxph iidCox.cph iidCox.phreg
###   iidCox.CauseSpecificCox

### ** Examples

library(survival)
library(data.table)
library(prodlim)
set.seed(10)
d <- sampleData(100, outcome = "survival")[,.(eventtime,event,X1,X6)]
setkey(d, eventtime)

m.cox <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
system.time(IF.cox <- iidCox(m.cox))
system.time(IF.cox_approx <- iidCox(m.cox, store.iid = "approx"))


IF.cox.all <- iidCox(m.cox, tau.hazard = sort(unique(c(7,d$eventtime))))
IF.cox.beta <- iidCox(m.cox, baseline.iid = FALSE)




cleanEx()
nameEx("influenceTest")
### * influenceTest

flush(stderr()); flush(stdout())

### Name: influenceTest
### Title: Influence test [Experimental!!]
### Aliases: influenceTest influenceTest.list influenceTest.default

### ** Examples

library(lava)
library(survival)
library(prodlim)
library(data.table)
n <- 100

#### Under H1
set.seed(1)
newdata <- data.frame(X1=0:1)

## simulate non proportional hazard using lava
m <- lvm()
regression(m) <- y ~ 1
regression(m) <- s ~ exp(-2*X1)
distribution(m,~X1) <- binomial.lvm()
distribution(m,~cens) <- coxWeibull.lvm(scale=1)
distribution(m,~y) <- coxWeibull.lvm(scale=1,shape=~s)
eventTime(m) <- eventtime ~ min(y=1,cens=0)
d <- as.data.table(sim(m,n))
setkey(d, eventtime)

## fit cox models
m.cox <- coxph(Surv(eventtime, status) ~ X1, 
               data = d, y = TRUE, x = TRUE)

mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1), 
                     data = d, y = TRUE, x = TRUE)

## compare models
# one time point
outIF <- influenceTest(list(m.cox, mStrata.cox), 
              type = "survival", newdata = newdata, times = 0.5)
confint(outIF)
                                 
# several timepoints
outIF <- influenceTest(list(m.cox, mStrata.cox), 
              type = "survival", newdata = newdata, times = c(0.5,1,1.5))
confint(outIF)

#### Under H0 (Cox) ####
set.seed(1)
## simulate proportional hazard using lava
m <- lvm()
regression(m) <- y ~ 1
distribution(m,~X1) <- binomial.lvm()
distribution(m,~cens) <- coxWeibull.lvm()
distribution(m,~y) <- coxWeibull.lvm()
eventTime(m) <- eventtime ~ min(y=1,cens=0)
d <- as.data.table(sim(m,n))
setkey(d, eventtime)

## fit cox models
Utime <- sort(unique(d$eventtime))
m.cox <- coxph(Surv(eventtime, status) ~ X1, 
               data = d, y = TRUE, x = TRUE)

mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1), 
                     data = d, y = TRUE, x = TRUE)

p.cox <- predictCox(m.cox, newdata = newdata, time = Utime, type = "survival")
p.coxStrata <- predictCox(mStrata.cox, newdata = newdata, time = Utime, type = "survival")

## display
library(ggplot2)
autoplot(p.cox)
autoplot(p.coxStrata)
 
## compare models
outIF <- influenceTest(list(m.cox, mStrata.cox), 
                       type = "survival", newdata = newdata, times = Utime[1:6])
confint(outIF)

#### Under H0 (CSC) ####
set.seed(1)
ff <- ~ f(X1,2) + f(X2,-0.033)
ff <- update(ff, ~ .+ f(X3,0) + f(X4,0) + f(X5,0))
ff <- update(ff, ~ .+ f(X6,0) + f(X7,0) + f(X8,0) + f(X9,0))
d <- sampleData(n, outcome = "competing.risk", formula = ff)
d[,X1:=as.numeric(as.character(X1))]
d[,X2:=as.numeric(as.character(X2))]
d[,X3:=as.numeric(as.character(X3))]
d[,X4:=as.numeric(as.character(X4))]
d[,X5:=as.numeric(as.character(X5))]
setkey(d, time)

Utime <- sort(unique(d$time))

## fit cox models
m.CSC <- CSC(Hist(time, event) ~ X1 + X2, data = d)
mStrata.CSC <- CSC(Hist(time, event) ~ strata(X1) + X2 + X3, data = d)

## compare models
outIF <- influenceTest(list(m.CSC, mStrata.CSC), 
             cause = 1, newdata = unique(d[,.(X1,X2,X3)]), times = Utime[1:5])
confint(outIF)



cleanEx()
nameEx("ipcw")
### * ipcw

flush(stderr()); flush(stdout())

### Name: ipcw
### Title: Estimation of censoring probabilities
### Aliases: ipcw ipcw.none ipcw.marginal ipcw.nonpar ipcw.cox ipcw.aalen
### Keywords: survival

### ** Examples


library(prodlim)
library(rms)
dat=SimSurv(30)

dat <- dat[order(dat$time),]

# using the marginal Kaplan-Meier for the censoring times

WKM=ipcw(Hist(time,status)~X2,
  data=dat,
  method="marginal",
  times=sort(unique(dat$time)),
  subject.times=dat$time,keep=c("fit"))
plot(WKM$fit)
WKM$fit

# using the Cox model for the censoring times given X2
library(survival)
WCox=ipcw(Hist(time=time,event=status)~X2,
  data=dat,
  method="cox",
  times=sort(unique(dat$time)),
  subject.times=dat$time,keep=c("fit"))
WCox$fit

plot(WKM$fit)
lines(sort(unique(dat$time)),
      1-WCox$IPCW.times[1,],
      type="l",
      col=2,
      lty=3,
      lwd=3)
lines(sort(unique(dat$time)),
      1-WCox$IPCW.times[5,],
      type="l",
      col=3,
      lty=3,
      lwd=3)

# using the stratified Kaplan-Meier
# for the censoring times given X2

WKM2=ipcw(Hist(time,status)~X2,
  data=dat,
  method="nonpar",
  times=sort(unique(dat$time)),
  subject.times=dat$time,keep=c("fit"))
plot(WKM2$fit,add=FALSE)





cleanEx()
nameEx("penalizedS3")
### * penalizedS3

flush(stderr()); flush(stdout())

### Name: penalizedS3
### Title: S3-wrapper for S4 function penalized
### Aliases: penalizedS3

### ** Examples

library(prodlim)
## Not run: 
##D ## too slow
##D library(penalized)
##D set.seed(8)
##D d <- sampleData(200,outcome="binary")
##D newd <- sampleData(80,outcome="binary")
##D fitridge <- penalizedS3(Y~X1+X2+pen(7:8), data=d, type="ridge",
##D standardize=TRUE, model="logistic",trace=FALSE)
##D fitlasso <- penalizedS3(Y~X1+X2+pen(7:8), data=d, type="lasso",
##D standardize=TRUE, model="logistic",trace=FALSE)
##D # fitnet <- penalizedS3(Y~X1+X2+pen(7:8), data=d, type="elastic.net",
##D # standardize=TRUE, model="logistic",trace=FALSE)
##D predictRisk(fitridge,newdata=newd)
##D predictRisk(fitlasso,newdata=newd)
##D # predictRisk(fitnet,newdata=newd)
##D Score(list(fitridge),data=newd,formula=Y~1)
##D Score(list(fitridge),data=newd,formula=Y~1,split.method="bootcv",B=2)
## End(Not run)
## Not run: 
##D  data(nki70) ## S4 fit
##D pen <- penalized(Surv(time, event), penalized = nki70[,8:77],
##D                  unpenalized = ~ER+Age+Diam+N+Grade, data = nki70,
##D lambda1 = 1)
##D penS3 <- penalizedS3(Surv(time,event)~ER+Age+Diam+pen(8:77)+N+Grade,
##D                      data=nki70, lambda1=1)
##D ## or
##D penS3 <- penalizedS3(Surv(time,event)~ER+pen(TSPYL5,Contig63649_RC)+pen(10:77)+N+Grade,
##D                      data=nki70, lambda1=1)
##D ## also this works
##D penS3 <- penalizedS3(Surv(time,event)~ER+Age+pen(8:33)+Diam+pen(34:77)+N+Grade,
##D                     data=nki70, lambda1=1)
## End(Not run)



cleanEx()
nameEx("plot.riskRegression")
### * plot.riskRegression

flush(stderr()); flush(stdout())

### Name: plot.riskRegression
### Title: Plotting predicted risk
### Aliases: plot.riskRegression
### Keywords: survival

### ** Examples


library(survival)
library(prodlim)
data(Melanoma)
fit.arr <- ARR(Hist(time,status)~invasion+age+strata(sex),data=Melanoma,cause=1)
plot(fit.arr,xlim=c(500,3000))





cleanEx()
nameEx("plotAUC")
### * plotAUC

flush(stderr()); flush(stdout())

### Name: plotAUC
### Title: Plot of time-dependent AUC curves
### Aliases: plotAUC

### ** Examples

library(survival)
library(prodlim)
d=sampleData(100,outcome="survival")
nd=sampleData(100,outcome="survival")
f1=coxph(Surv(time,event)~X1+X6+X8,data=d,x=TRUE,y=TRUE)
f2=coxph(Surv(time,event)~X2+X5+X9,data=d,x=TRUE,y=TRUE)
xx=Score(list("X1+X6+X8"=f1,"X2+X5+X9"=f2), formula=Surv(time,event)~1,
data=nd, metrics="auc", null.model=FALSE, times=seq(3:10))
aucgraph <- plotAUC(xx)
plotAUC(xx,conf.int=TRUE)
## difference between 
plotAUC(xx,which="contrasts",conf.int=TRUE)





cleanEx()
nameEx("plotBrier")
### * plotBrier

flush(stderr()); flush(stdout())

### Name: plotBrier
### Title: Plot Brier curve
### Aliases: plotBrier

### ** Examples

# survival
library(survival)
library(prodlim)
ds1=sampleData(40,outcome="survival")
ds2=sampleData(40,outcome="survival")
f1 <- coxph(Surv(time,event)~X1+X3+X5+X7+X9,data=ds1,x=TRUE)
f2 <- coxph(Surv(time,event)~X2+X4+X6+X8+X10,data=ds1,x=TRUE)
xscore <- Score(list(f1,f2),formula=Hist(time,event)~1,data=ds2,times=0:12,metrics="brier")
plotBrier(xscore)



cleanEx()
nameEx("plotCalibration")
### * plotCalibration

flush(stderr()); flush(stdout())

### Name: plotCalibration
### Title: Plot Calibration curve
### Aliases: plotCalibration

### ** Examples

library(prodlim)
# binary 
db=sampleData(100,outcome="binary")
fb1=glm(Y~X1+X5+X7,data=db,family="binomial")
fb2=glm(Y~X1+X3+X6+X7,data=db,family="binomial")
xb=Score(list(model1=fb1,model2=fb2),Y~1,data=db,
          plots="cal")
plotCalibration(xb,brier.in.legend=TRUE)
plotCalibration(xb,bars=TRUE,model="model1")
plotCalibration(xb,models=1,bars=TRUE,names.cex=1.3)

# survival
library(survival)
library(prodlim)
dslearn=sampleData(56,outcome="survival")
dstest=sampleData(100,outcome="survival")
fs1=coxph(Surv(time,event)~X1+X5+X7,data=dslearn,x=1)
fs2=coxph(Surv(time,event)~strata(X1)+X3+X6+X7,data=dslearn,x=1)
xs=Score(list(Cox1=fs1,Cox2=fs2),Surv(time,event)~1,data=dstest,
          plots="cal",metrics=NULL)
plotCalibration(xs)
plotCalibration(xs,cens.method="local",pseudo=1)
plotCalibration(xs,method="quantile")


# competing risks

## Not run: 
##D data(Melanoma)
##D f1 <- CSC(Hist(time,status)~age+sex+epicel+ulcer,data=Melanoma)
##D f2 <- CSC(Hist(time,status)~age+sex+logthick+epicel+ulcer,data=Melanoma)
##D x <- Score(list(model1=f1,model2=f2),Hist(time,status)~1,data=Melanoma,
##D            cause= 2,times=5*365.25,plots="cal")
##D plotCalibration(x)
## End(Not run)




cleanEx()
nameEx("plotEffects")
### * plotEffects

flush(stderr()); flush(stdout())

### Name: plotEffects
### Title: Plotting time-varying effects from a risk regression model.
### Aliases: plotEffects
### Keywords: survival

### ** Examples


library(survival)
library(prodlim)
data(Melanoma)

fit.tarr <- ARR(Hist(time,status)~strata(sex),
                data=Melanoma,
                cause=1)
plotEffects(fit.tarr)

fit.tarr <- ARR(Hist(time,status)~strata(sex)+strata(invasion),
                data=Melanoma,
                cause=1,
                times=seq(800,3000,20))
plotEffects(fit.tarr,formula=~sex)
plotEffects(fit.tarr,formula=~invasion)
plotEffects(fit.tarr,
            formula=~invasion,
            level="invasionlevel.1")

## legend arguments are transcluded:
plotEffects(fit.tarr,
            formula=~invasion,
            legend.bty="b",
            legend.cex=1)

## and other smart arguments too:
plotEffects(fit.tarr,
	    formula=~invasion,
	    legend.bty="b",
axis2.las=2,
	    legend.cex=1)





cleanEx()
nameEx("plotPredictRisk")
### * plotPredictRisk

flush(stderr()); flush(stdout())

### Name: plotPredictRisk
### Title: Plotting predicted risks curves.
### Aliases: plotPredictRisk
### Keywords: survival

### ** Examples

library(survival)
# generate survival data
# no effect
set.seed(8)
d <- sampleData(80,outcome="survival",formula = ~f(X6, 0) + f(X7, 0))
d[,table(event)]
f <- coxph(Surv(time,event)~X6+X7,data=d,x=1)
plotPredictRisk(f)

# large effect
set.seed(8)
d <- sampleData(80,outcome="survival",formula = ~f(X6, 0.1) + f(X7, -0.1))
d[,table(event)]
f <- coxph(Surv(time,event)~X6+X7,data=d,x=1)
plotPredictRisk(f)

# generate competing risk data
# small effect
set.seed(8)
d <- sampleData(40,formula = ~f(X6, 0.01) + f(X7, -0.01))
d[,table(event)]
f <- CSC(Hist(time,event)~X5+X6,data=d)
plotPredictRisk(f)

# large effect
set.seed(8)
d <- sampleData(40,formula = ~f(X6, 0.1) + f(X7, -0.1))
d[,table(event)]
f <- CSC(Hist(time,event)~X5+X6,data=d)
plotPredictRisk(f)



cleanEx()
nameEx("plotROC")
### * plotROC

flush(stderr()); flush(stdout())

### Name: plotROC
### Title: Plot ROC curves
### Aliases: plotROC

### ** Examples

## binary
set.seed(18)
library(randomForest)
library(prodlim)
bdl <- sampleData(40,outcome="binary")
bdt <- sampleData(58,outcome="binary")
bdl[,y:=factor(Y)]
bdt[,y:=factor(Y)]
fb1 <- glm(y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=bdl,family="binomial")
fb2 <- randomForest(y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=bdl)
xb <- Score(list("glm"=fb1,"rf"=fb2),y~1,data=bdt,
            plots="roc",metrics=c("auc","brier"))
plotROC(xb,brier.in.legend=1L)

# with cross-validation
## Not run: 
##D xb3 <- Score(list("glm"=fb1,"rf"=fb2),y~1,data=bdl,
##D             plots="roc",B=3,split.method="bootcv",
##D             metrics=c("auc"))
## End(Not run)
## survival
set.seed(18)
library(survival)
sdl <- sampleData(40,outcome="survival")
sdt <- sampleData(58,outcome="survival")
fs1 <- coxph(Surv(time,event)~X3+X5+X6+X7+X8+X10,data=sdl,x=TRUE)
fs2 <- coxph(Surv(time,event)~X1+X2+X9,data=sdl,x=TRUE)
xs <- Score(list(model1=fs1,model2=fs2),Hist(time,event)~1,data=sdt,
            times=5,plots="roc",metrics="auc")
plotROC(xs)
## competing risks
data(Melanoma)
f1 <- CSC(Hist(time,status)~age+sex+epicel+ulcer,data=Melanoma)
f2 <- CSC(Hist(time,status)~age+sex+logthick+epicel+ulcer,data=Melanoma)
x <- Score(list(model1=f1,model2=f2),Hist(time,status)~1,data=Melanoma,
            cause=1,times=5*365.25,plots="roc",metrics="auc")
plotROC(x)



cleanEx()
nameEx("plotRisk")
### * plotRisk

flush(stderr()); flush(stdout())

### Name: plotRisk
### Title: plot predicted risks
### Aliases: plotRisk

### ** Examples

library(prodlim)
## uncensored
learndat = sampleData(40,outcome="binary")
testdat = sampleData(40,outcome="binary")
lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family="binomial")
lr2 = glm(Y~X3+X5+X6,data=learndat,family="binomial")
xb=Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5+X6)"=lr2),formula=Y~1,
         data=testdat,summary="risks",null.model=0L)
plotRisk(xb)
## survival
library(survival)
learndat = sampleData(40,outcome="survival")
testdat = sampleData(40,outcome="survival")
cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=learndat,x=TRUE)
cox2 = coxph(Surv(time,event)~X3+X5+X6,data=learndat,x=TRUE)
xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),formula=Surv(time,event)~1,
         data=testdat,summary="risks",null.model=0L,times=c(3,5,6))
plotRisk(xs,times=5)
## competing risk
## Not run: 
##D library(prodlim)
##D library(survival)
##D set.seed(8)
##D learndat = sampleData(80,outcome="competing.risk")
##D testdat = sampleData(140,outcome="competing.risk")
##D m1 = FGR(Hist(time,event)~X2+X7+X9,data=learndat,cause=1)
##D m2 = CSC(Hist(time,event)~X2+X7+X9,data=learndat,cause=1)
##D xcr=Score(list("FGR"=m1,"CSC"=m2),formula=Hist(time,event)~1,
##D          data=testdat,summary="risks",null.model=0L,times=c(3,5))
##D plotRisk(xcr,times=1)
## End(Not run)



cleanEx()
nameEx("predict.CauseSpecificCox")
### * predict.CauseSpecificCox

flush(stderr()); flush(stdout())

### Name: predict.CauseSpecificCox
### Title: Predicting Absolute Risk from Cause-Specific Cox Models
### Aliases: predict.CauseSpecificCox predictBig.CauseSpecificCox

### ** Examples

library(survival)
library(prodlim)
#### generate data ####
set.seed(5)
d <- sampleData(80,outcome="comp") ## training dataset
nd <- sampleData(4,outcome="comp") ## validation dataset
d$time <- round(d$time,1) ## create tied events
ttt <- sort(sample(x = unique(d$time), size = 10))

## estimate a CSC model based on the coxph function
CSC.fit <- CSC(Hist(time,event)~ X3+X8, data=d, method = "breslow")

## compute the absolute risk of cause 1, in the validation dataset
## at time 1:10
CSC.risk <-  predict(CSC.fit, newdata=nd, times=1:10, cause=1)
CSC.risk

## compute absolute risks with CI for cause 2
## (without displaying the value of the covariates)
predict(CSC.fit,newdata=nd,times=1:10,cause=2,se=TRUE,
        keep.newdata = FALSE)

## other example
library(survival)
CSC.fit.s <- CSC(list(Hist(time,event)~ strata(X1)+X2+X9,
 Hist(time,event)~ X2+strata(X4)+X8+X7),data=d, method = "breslow")
predict(CSC.fit.s,cause=1,times=ttt,se=1L) ## note: absRisk>1 due to small number of observations

## using the cph function instead of coxph
CSC.cph <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow", fitter = "cph")#' 
predict(CSC.cph, newdata = d, cause = 2, times = ttt)

## landmark analysis
T0 <- 1
predCSC_afterT0 <- predict(CSC.fit, newdata = d, cause = 2, times = ttt[ttt>T0], landmark = T0)
predCSC_afterT0



cleanEx()
nameEx("predict.FGR")
### * predict.FGR

flush(stderr()); flush(stdout())

### Name: predict.FGR
### Title: Predict subject specific risks (cumulative incidence) based on
###   Fine-Gray regression model
### Aliases: predict.FGR

### ** Examples

library(prodlim)
library(survival)
set.seed(10)
d <- sampleData(101, outcome = "competing.risk")
tFun<-function(t) {t}
fgr<-FGR(Hist(time, event)~X1+strata(X2)+X6+cov2(X7, tf=tFun),
         data=d, cause=1)
predictRisk(fgr,times=5,newdata=d[1:10])



cleanEx()
nameEx("predict.riskRegression")
### * predict.riskRegression

flush(stderr()); flush(stdout())

### Name: predict.riskRegression
### Title: Predict individual risk.
### Aliases: predict.riskRegression
### Keywords: survival

### ** Examples


data(Melanoma)
library(prodlim)
library(survival)

fit.tarr <- ARR(Hist(time,status)~age+invasion+strata(sex),data=Melanoma,cause=1)
predict(fit.tarr,newdata=data.frame(age=48,
                     invasion=factor("level.1",
                         levels=levels(Melanoma$invasion)),
                     sex=factor("Female",levels=levels(Melanoma$sex))))
predict(fit.tarr,newdata=data.frame(age=48,
                     invasion=factor("level.1",
                         levels=levels(Melanoma$invasion)),
                     sex=factor("Male",levels=levels(Melanoma$sex))))
predict(fit.tarr,newdata=data.frame(age=c(48,58,68),
                     invasion=factor("level.1",
                         levels=levels(Melanoma$invasion)),
                     sex=factor("Male",levels=levels(Melanoma$sex))))
predict(fit.tarr,newdata=Melanoma[1:4,])




cleanEx()
nameEx("predictCox")
### * predictCox

flush(stderr()); flush(stdout())

### Name: predictCox
### Title: Fast computation of survival probabilities, hazards and
###   cumulative hazards from Cox regression models
### Aliases: predictCox

### ** Examples


library(survival)
library(data.table)

#### generate data ####
set.seed(10)
d <- sampleData(40,outcome="survival") ## training dataset
nd <- sampleData(4,outcome="survival") ## validation dataset
d$time <- round(d$time,1) ## create tied events
# table(duplicated(d$time))

#### stratified Cox model ####
fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
             data=d, ties="breslow", x = TRUE, y = TRUE)

## compute the baseline cumulative hazard
fit.haz <- predictCox(fit)
cbind(survival::basehaz(fit), fit.haz$cumhazard)

## compute individual specific cumulative hazard and survival probabilities 
fit.pred <- predictCox(fit, newdata=nd, times=c(3,8), se = TRUE, band = TRUE)
fit.pred

####  other examples ####
# one strata variable
fitS <- coxph(Surv(time,event)~strata(X1)+X2,
              data=d, ties="breslow", x = TRUE, y = TRUE)

predictCox(fitS)
predictCox(fitS, newdata=nd, times = 1)

# two strata variables
set.seed(1)
d$U=sample(letters[1:5],replace=TRUE,size=NROW(d))
d$V=sample(letters[4:10],replace=TRUE,size=NROW(d))
nd$U=sample(letters[1:5],replace=TRUE,size=NROW(nd))
nd$V=sample(letters[4:10],replace=TRUE,size=NROW(nd))
fit2S <- coxph(Surv(time,event)~X1+strata(U)+strata(V)+X2,
              data=d, ties="breslow", x = TRUE, y = TRUE)

cbind(survival::basehaz(fit2S),predictCox(fit2S,type="cumhazard")$cumhazard)
predictCox(fit2S)
predictCox(fitS, newdata=nd, times = 3)

# left truncation
test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8), 
              stop=c(2,3,6,7,8,9,9,9,14,17), 
              event=c(1,1,1,1,1,1,1,0,0,0), 
              x=c(1,0,0,1,0,1,1,1,0,0)) 
m.cph <- coxph(Surv(start, stop, event) ~ 1, test2, x = TRUE)
as.data.table(predictCox(m.cph))

basehaz(m.cph)



cleanEx()
nameEx("predictCoxPL")
### * predictCoxPL

flush(stderr()); flush(stdout())

### Name: predictCoxPL
### Title: Computation of survival probabilities from Cox regression models
###   using the product limit estimator.
### Aliases: predictCoxPL

### ** Examples

library(survival)

#### generate data ####
set.seed(10)
d <- sampleData(40,outcome="survival")
nd <- sampleData(4,outcome="survival")
d$time <- round(d$time,1)

#### Cox model ####
fit <- coxph(Surv(time,event)~ X1 + X2 + X6,
             data=d, ties="breslow", x = TRUE, y = TRUE)

## exponential approximation
predictCox(fit, newdata = d, times = 1:5)

## product limit
predictCoxPL(fit, newdata = d, times = 1:5)

#### stratified Cox model ####
fitS <- coxph(Surv(time,event)~ X1 + strata(X2) + X6,
             data=d, ties="breslow", x = TRUE, y = TRUE)

## exponential approximation
predictCox(fitS, newdata = d, times = 1:5)

## product limit
predictCoxPL(fitS, newdata = d, times = 1:5)

#### fully stratified Cox model ####
fitS <- coxph(Surv(time,event)~ 1,
             data=d, ties="breslow", x = TRUE, y = TRUE)

## product limit
GS <- survfit(Surv(time,event)~1, data = d)
range(predictCoxPL(fitS)$survival - GS$surv)

fitS <- coxph(Surv(time,event)~ strata(X2),
             data=d, ties="breslow", x = TRUE, y = TRUE)

## product limit
GS <- survfit(Surv(time,event)~X2, data = d)
range(predictCoxPL(fitS)$survival - GS$surv)




cleanEx()
nameEx("predictRisk")
### * predictRisk

flush(stderr()); flush(stdout())

### Name: predictRisk
### Title: Extrating predicting risks from regression models
### Aliases: predictRisk predictRisk.CauseSpecificCox
###   predictRisk.riskRegression predictRisk.FGR predictRisk.prodlim
###   predictRisk.rfsrc predictRisk.aalen predictRisk.ARR
###   predictRisk.cox.aalen predictRisk.coxph predictRisk.cph
###   predictRisk.default predictRisk.matrix predictRisk.pecCtree
###   predictRisk.pecCforest predictRisk.psm predictRisk.selectCox
###   predictRisk.survfit predictRisk.randomForest predictRisk.lrm
###   predictRisk.glm predictRisk.rpart predictRisk.gbm
###   predictRisk.flexsurvreg predictRisk.double predictRisk.integer
###   predictRisk.factor predictRisk.numeric predictRisk.formula
###   predictRisk.BinaryTree predictRisk.coxphTD predictRisk.CSCTD
###   predictRisk.coxph.penal predictRisk.ranger predictRisk.penfitS3
###   predictRisk.SuperPredictor
### Keywords: survival

### ** Examples

## binary outcome
library(rms)
set.seed(7)
d <- sampleData(80,outcome="binary")
nd <- sampleData(80,outcome="binary")
fit <- lrm(Y~X1+X8,data=d)
predictRisk(fit,newdata=nd)
## Not run: 
##D library(SuperLearner)
##D set.seed(1)
##D sl = SuperLearner(Y = d$Y, X = d[,-1], family = binomial(),
##D       SL.library = c("SL.mean", "SL.glmnet", "SL.randomForest"))
## End(Not run)

## survival outcome
# generate survival data
library(prodlim)
set.seed(100)
d <- sampleData(100,outcome="survival")
d[,X1:=as.numeric(as.character(X1))]
d[,X2:=as.numeric(as.character(X2))]
# then fit a Cox model
library(rms)
cphmodel <- cph(Surv(time,event)~X1+X2,data=d,surv=TRUE,x=TRUE,y=TRUE)
# or via survival
library(survival)
coxphmodel <- coxph(Surv(time,event)~X1+X2,data=d,x=TRUE,y=TRUE)

# Extract predicted survival probabilities 
# at selected time-points:
ttt <- quantile(d$time)
# for selected predictor values:
ndat <- data.frame(X1=c(0.25,0.25,-0.05,0.05),X2=c(0,1,0,1))
# as follows
predictRisk(cphmodel,newdata=ndat,times=ttt)
predictRisk(coxphmodel,newdata=ndat,times=ttt)

# stratified cox model
sfit <- coxph(Surv(time,event)~strata(X1)+X2,data=d,x=TRUE,y=TRUE)
predictRisk(sfit,newdata=d[1:3,],times=c(1,3,5,10))

## simulate learning and validation data
learndat <- sampleData(100,outcome="survival")
valdat <- sampleData(100,outcome="survival")
## use the learning data to fit a Cox model
library(survival)
fitCox <- coxph(Surv(time,event)~X1+X2,data=learndat,x=TRUE,y=TRUE)
## suppose we want to predict the survival probabilities for all subjects
## in the validation data at the following time points:
## 0, 12, 24, 36, 48, 60
psurv <- predictRisk(fitCox,newdata=valdat,times=seq(0,60,12))
## This is a matrix with event probabilities (1-survival)
## one column for each of the 5 time points
## one row for each validation set individual

# Do the same for a randomSurvivalForest model
# library(randomForestSRC)
# rsfmodel <- rfsrc(Surv(time,event)~X1+X2,data=learndat)
# prsfsurv=predictRisk(rsfmodel,newdata=valdat,times=seq(0,60,12))
# plot(psurv,prsfsurv)

## Cox with ridge option
f1 <- coxph(Surv(time,event)~X1+X2,data=learndat,x=TRUE,y=TRUE)
f2 <- coxph(Surv(time,event)~ridge(X1)+ridge(X2),data=learndat,x=TRUE,y=TRUE)
## Not run: 
##D plot(predictRisk(f1,newdata=valdat,times=10),
##D      riskRegression:::predictRisk.coxph(f2,newdata=valdat,times=10),
##D      xlim=c(0,1),
##D      ylim=c(0,1),
##D      xlab="Unpenalized predicted survival chance at 10",
##D      ylab="Ridge predicted survival chance at 10")
## End(Not run)

## competing risks

library(survival)
library(riskRegression)
library(prodlim)
train <- prodlim::SimCompRisk(100)
test <- prodlim::SimCompRisk(10)
cox.fit  <- CSC(Hist(time,cause)~X1+X2,data=train)
predictRisk(cox.fit,newdata=test,times=seq(1:10),cause=1)

## with strata
cox.fit2  <- CSC(list(Hist(time,cause)~strata(X1)+X2,Hist(time,cause)~X1+X2),data=train)
predictRisk(cox.fit2,newdata=test,times=seq(1:10),cause=1)




cleanEx()
nameEx("riskLevelPlot")
### * riskLevelPlot

flush(stderr()); flush(stdout())

### Name: riskLevelPlot
### Title: Level plots for risk prediction models
### Aliases: riskLevelPlot

### ** Examples


# ---------- logistic regression --------------------
expit <- function(x){exp(x)/(1+exp(x))}
partyData <- function(N){
  Age <- runif(N,.5,15)
  Parasites <- rnorm(N,mean=3.5-0.03*Age)
  Fever <- factor(rbinom(N,1,expit(-3.5-.3*Age+.55*Parasites+0.15*Age*Parasites)))
  data.frame(Fever,Age,Parasites)
}
d <- partyData(100)
f <- glm(Fever~Age+Parasites,data=d,family="binomial")
riskLevelPlot(f,Fever~Age+Parasites,d)
library(randomForest)
rf <- randomForest(Fever~Age+Parasites,data=d)
riskLevelPlot(f,Fever~Age+Parasites,d)
riskLevelPlot(rf,Fever~Age+Parasites,d)

# ---------- survival analysis --------------------

# --simulate an artificial data frame
# with survival response and three predictors

library(survival)
library(prodlim)
set.seed(140515)
sdat <- sampleData(43,outcome="survival")
# -- fit a Cox regression model 
survForm = Surv(time,event) ~ X8 + X9
cox <- coxph(survForm, data = sdat,x=TRUE)

# --choose a time horizon for the predictions and plot the risks
timeHorizon <- floor(median(sdat$time))
riskLevelPlot(cox, survForm, data = sdat, horizon = timeHorizon)

# ---------- competing risks --------------------

# -- simulate an artificial data frame
# with competing cause response and three predictors
library(cmprsk)
library(riskRegression)
set.seed(140515)
crdat <- sampleData(49)

# -- fit a cause-specific Cox regression model
crForm <- Hist(time,event)~X8+X9
csCox  <- CSC(crForm, data=crdat)

# -- choose a time horizon and plot the risk for a given cause
timeHorizon <- floor(median(crdat$time))
riskLevelPlot(csCox, crForm, data = crdat, horizon = timeHorizon, cause = 1)




cleanEx()
nameEx("riskRegression")
### * riskRegression

flush(stderr()); flush(stdout())

### Name: riskRegression
### Title: Risk Regression Fits a regression model for the risk of an event
###   - allowing for competing risks.
### Aliases: riskRegression ARR LRR
### Keywords: survival

### ** Examples


library(prodlim)
data(Melanoma,package="riskRegression")
## tumor thickness on the log-scale
Melanoma$logthick <- log(Melanoma$thick)

# Single binary factor

## absolute risk regression
library(survival)
library(prodlim)
fit.arr <- ARR(Hist(time,status)~sex,data=Melanoma,cause=1)
print(fit.arr)
# show predicted cumulative incidences
plot(fit.arr,col=3:4,newdata=data.frame(sex=c("Female","Male")))

## compare with non-parametric Aalen-Johansen estimate
library(prodlim)
fit.aj <- prodlim(Hist(time,status)~sex,data=Melanoma)
plot(fit.aj,conf.int=FALSE)
plot(fit.arr,add=TRUE,col=3:4,newdata=data.frame(sex=c("Female","Male")))

## with time-dependent effect
fit.tarr <- ARR(Hist(time,status)~strata(sex),data=Melanoma,cause=1)
plot(fit.tarr,newdata=data.frame(sex=c("Female","Male")))

## logistic risk regression
fit.lrr <- LRR(Hist(time,status)~sex,data=Melanoma,cause=1)
summary(fit.lrr)


# Single continuous factor

## tumor thickness on the log-scale
Melanoma$logthick <- log(Melanoma$thick)

## absolute risk regression 
fit2.arr <- ARR(Hist(time,status)~logthick,data=Melanoma,cause=1)
print(fit2.arr)
# show predicted cumulative incidences
plot(fit2.arr,col=1:5,newdata=data.frame(logthick=quantile(Melanoma$logthick)))

## comparison with nearest neighbor non-parametric Aalen-Johansen estimate
library(prodlim)
fit2.aj <- prodlim(Hist(time,status)~logthick,data=Melanoma)
plot(fit2.aj,conf.int=FALSE,newdata=data.frame(logthick=quantile(Melanoma$logthick)))
plot(fit2.arr,add=TRUE,col=1:5,lty=3,newdata=data.frame(logthick=quantile(Melanoma$logthick)))

## logistic risk regression
fit2.lrr <- LRR(Hist(time,status)~logthick,data=Melanoma,cause=1)
summary(fit2.lrr)

## change model for censoring weights
library(rms)
fit2a.lrr <- LRR(Hist(time,status)~logthick,
                 data=Melanoma,
                 cause=1,
                 cens.model="cox",
                 cens.formula=~sex+epicel+ulcer+age+logthick)
summary(fit2a.lrr)

##  compare prediction performance
Score(list(ARR=fit2.arr,AJ=fit2.aj,LRR=fit2.lrr),formula=Hist(time,status)~1,data=Melanoma)


# multiple regression
library(riskRegression)
library(prodlim)
# absolute risk model
multi.arr <- ARR(Hist(time,status)~logthick+sex+age+ulcer,data=Melanoma,cause=1)

# stratified model allowing different baseline risk for the two gender
multi.arr <- ARR(Hist(time,status)~thick+strata(sex)+age+ulcer,data=Melanoma,cause=1)

# stratify by a continuous variable: strata(age)
multi.arr <- ARR(Hist(time,status)~tp(thick,power=0)+strata(age)+sex+ulcer,
                 data=Melanoma,
                 cause=1)

fit.arr2a <- ARR(Hist(time,status)~tp(thick,power=1),data=Melanoma,cause=1)
summary(fit.arr2a)
fit.arr2b <- ARR(Hist(time,status)~timevar(thick),data=Melanoma,cause=1)
summary(fit.arr2b)

## logistic risk model
fit.lrr <- LRR(Hist(time,status)~thick,data=Melanoma,cause=1)
summary(fit.lrr)





## nearest neighbor non-parametric Aalen-Johansen estimate
library(prodlim)
fit.aj <- prodlim(Hist(time,status)~thick,data=Melanoma)
plot(fit.aj,conf.int=FALSE)

# prediction performance
x <- Score(list(fit.arr2a,fit.arr2b,fit.lrr),
             data=Melanoma,
             formula=Hist(time,status)~1,
             cause=1,
             split.method="none")





cleanEx()
nameEx("riskRegression.options")
### * riskRegression.options

flush(stderr()); flush(stdout())

### Name: riskRegression.options
### Title: Global options for 'riskRegression'
### Aliases: riskRegression.options

### ** Examples

options <- riskRegression.options()

## add new method.predictRiskIID
riskRegression.options(method.predictRiskIID = c(options$method.predictRiskIID,"xx"))

riskRegression.options()



cleanEx()
nameEx("rowCenter_cpp")
### * rowCenter_cpp

flush(stderr()); flush(stdout())

### Name: rowCenter_cpp
### Title: Apply - by row
### Aliases: rowCenter_cpp

### ** Examples

x <- matrix(1,6,5)
sweep(x, MARGIN = 2, FUN = "-", STATS = 1:5)
rowCenter_cpp(x, 1:5 )

rowCenter_cpp(x, colMeans(x) )



cleanEx()
nameEx("rowCumProd")
### * rowCumProd

flush(stderr()); flush(stdout())

### Name: rowCumProd
### Title: Apply cumprod in each row
### Aliases: rowCumProd

### ** Examples

x <- matrix(1:8,ncol=2)
rowCumProd(x)



cleanEx()
nameEx("rowCumSum")
### * rowCumSum

flush(stderr()); flush(stdout())

### Name: rowCumSum
### Title: Apply cumsum in each row
### Aliases: rowCumSum

### ** Examples

x <- matrix(1:8,ncol=2)
rowCumSum(x)



cleanEx()
nameEx("rowMultiply_cpp")
### * rowMultiply_cpp

flush(stderr()); flush(stdout())

### Name: rowMultiply_cpp
### Title: Apply * by row
### Aliases: rowMultiply_cpp

### ** Examples

x <- matrix(1,6,5)
sweep(x, MARGIN = 2, FUN = "*", STATS = 1:5)
rowMultiply_cpp(x, 1:5 )

rowMultiply_cpp(x, 1/colMeans(x) )




cleanEx()
nameEx("rowPaste")
### * rowPaste

flush(stderr()); flush(stdout())

### Name: rowPaste
### Title: Collapse Rows of Characters.
### Aliases: rowPaste

### ** Examples

## Not run: 
##D M <- matrix(letters,nrow = 26, ncol = 2)
##D rowPaste(M)
## End(Not run)



cleanEx()
nameEx("rowScale_cpp")
### * rowScale_cpp

flush(stderr()); flush(stdout())

### Name: rowScale_cpp
### Title: Apply / by row
### Aliases: rowScale_cpp

### ** Examples

x <- matrix(1,6,5)
sweep(x, MARGIN = 2, FUN = "/", STATS = 1:5)
rowScale_cpp(x, 1:5 )

rowScale_cpp(x, colMeans(x) )



cleanEx()
nameEx("rowSumsCrossprod")
### * rowSumsCrossprod

flush(stderr()); flush(stdout())

### Name: rowSumsCrossprod
### Title: Apply crossprod and rowSums
### Aliases: rowSumsCrossprod

### ** Examples

x <- matrix(1:10,nrow=5)
y <- matrix(1:20,ncol=4)
rowSumsCrossprod(x,y,0)

x <- matrix(1:10,nrow=5)
y <- matrix(1:20,ncol=5)
rowSumsCrossprod(x,y,1)



cleanEx()
nameEx("sampleData")
### * sampleData

flush(stderr()); flush(stdout())

### Name: sampleData
### Title: Simulate data with binary or time-to-event outcome
### Aliases: sampleData sampleDataTD

### ** Examples

sampleData(10,outcome="binary")
sampleData(10,outcome="survival")
sampleData(10,outcome="competing.risks")



cleanEx()
nameEx("selectCox")
### * selectCox

flush(stderr()); flush(stdout())

### Name: selectCox
### Title: Backward variable selection in the Cox regression model
### Aliases: selectCox
### Keywords: survival

### ** Examples


library(pec)
library(prodlim)
data(GBSG2)
library(survival)
f <- selectCox(Surv(time,cens)~horTh+age+menostat+tsize+tgrade+pnodes+progrec+estrec ,
	       data=GBSG2)




cleanEx()
nameEx("simActiveSurveillance")
### * simActiveSurveillance

flush(stderr()); flush(stdout())

### Name: simActiveSurveillance
### Title: Simulate data of a hypothetical active surveillance prostate
###   cancer study
### Aliases: simActiveSurveillance

### ** Examples

set.seed(71)
simActiveSurveillance(3)



cleanEx()
nameEx("simMelanoma")
### * simMelanoma

flush(stderr()); flush(stdout())

### Name: simMelanoma
### Title: Simulate data alike the Melanoma data
### Aliases: simMelanoma

### ** Examples

set.seed(71)
simMelanoma(3)



cleanEx()
nameEx("sliceMultiply_cpp")
### * sliceMultiply_cpp

flush(stderr()); flush(stdout())

### Name: sliceMultiply_cpp
### Title: Apply * by slice
### Aliases: sliceMultiply_cpp sliceMultiplyPointer_cpp

### ** Examples

x <- array(1, dim = c(2,6,5))
M <- matrix(1:12,2,6)
sweep(x, MARGIN = 1:2, FUN = "*", STATS = M)
sliceMultiply_cpp(x, M) 




cleanEx()
nameEx("sliceScale_cpp")
### * sliceScale_cpp

flush(stderr()); flush(stdout())

### Name: sliceScale_cpp
### Title: Apply / by slice
### Aliases: sliceScale_cpp sliceScalePointer_cpp

### ** Examples

x <- array(1, dim = c(2,6,5))
M <- matrix(1:12,2,6)
sweep(x, MARGIN = 1:2, FUN = "/", STATS = M)
sliceScale_cpp(x, M) 




cleanEx()
nameEx("subjectWeights")
### * subjectWeights

flush(stderr()); flush(stdout())

### Name: subjectWeights
### Title: Estimation of censoring probabilities at subject specific times
### Aliases: subjectWeights subjectWeights.none subjectWeights.km
###   subjectWeights.marginal subjectWeights.nonpar subjectWeights.cox
###   subjectWeights.aalen
### Keywords: survival

### ** Examples


library(prodlim)
library(survival)
dat=SimSurv(300)

dat <- dat[order(dat$time,-dat$status),]

# using the marginal Kaplan-Meier for the censoring times

WKM=subjectWeights(Hist(time,status)~X2,data=dat,method="marginal")
plot(WKM$fit)
WKM$fit
WKM$weights

# using the Cox model for the censoring times given X2

WCox=subjectWeights(Surv(time,status)~X2,data=dat,method="cox")
WCox
plot(WCox$weights,WKM$weights)

# using the stratified Kaplan-Meier for the censoring times given X2

WKM2 <- subjectWeights(Surv(time,status)~X2,data=dat,method="nonpar")
plot(WKM2$fit,add=FALSE)





cleanEx()
nameEx("subsetIndex")
### * subsetIndex

flush(stderr()); flush(stdout())

### Name: subsetIndex
### Title: Extract Specific Elements From An Object
### Aliases: subsetIndex subsetIndex.default subsetIndex.matrix

### ** Examples

M <- matrix(rnorm(50),5,10)
subsetIndex(M, index = c(0,0,1), default = 0)
subsetIndex(M, index = c(0,2,3,NA), default = 0)
subsetIndex(M, index = c(0,NA,2,3,NA), default = 0)

C <- 1:10
subsetIndex(C, index = c(0,0,1,5,NA), default = 0)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
