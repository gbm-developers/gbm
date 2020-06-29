pkgname <- "Plasmode"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Plasmode')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("PlasmodeBin")
### * PlasmodeBin

flush(stderr()); flush(stdout())

### Name: PlasmodeBin
### Title: Performs the plasmode simulation
### Aliases: PlasmodeBin

### ** Examples

{
library(mgcv)
library(nlme)
library(glm2)
library(arm)
library(MASS)
library(lme4)
library(epiDisplay)
library(foreign)
library(nnet)

data("Compaq")
levels(Compaq$stage) <- c(1,2,3,4)
Compaq$stage<-as.numeric(levels(Compaq$stage))[Compaq$stage]
## Creating the binary exposure variable
Compaq$exposure<-ifelse(Compaq$hospital == "Public hospital",1,0)
## Creating binary variables for some confounders
Compaq$ses1<-ifelse(Compaq$ses == "Poor",1,0)
Compaq$ses2<-ifelse(Compaq$ses == "Poor-middle",1,0)
Compaq$ses3<-ifelse(Compaq$ses == "High-middle",1,0)

Compaq$age1<-ifelse(Compaq$agegr == "<40",1,0)
Compaq$age2<-ifelse(Compaq$agegr == "40-49",1,0)
Compaq$age3<-ifelse(Compaq$agegr == "50-59",1,0)


## Creating the formulas for the outcome and the exposure model
form1<- status~ exposure + stage + ses1 + ses2 + ses3 + age1 + age2 + age3
form2<- exposure ~ stage + ses1 + ses2 + ses3 + age1 + age2 + age3

set.seed(111)
Bin_Form1<-PlasmodeBin(formulaOut=form1, objectOut=NULL,formulaExp=form2,
                     objectExp= NULL,data=Compaq,idVar="id",effectOR =1,
                     MMOut=c(1,1,2,1,1,2,1,2),MMExp=c(1,1,1,1,1,1,1),
                     nsim=2, size=nrow(Compaq), eventRate=NULL, exposedPrev=NULL)

Bin_Form2<-PlasmodeBin(formulaOut=form1, objectOut=NULL,formulaExp=NULL,
                      objectExp= NULL,data=Compaq,idVar="id",effectOR =1,
                      MMOut=c(1,1,2,1,1,2,1,2),MMExp=1, nsim=2,
                      size=nrow(Compaq), eventRate=NULL, exposedPrev=NULL)

Bin_Form3<-PlasmodeBin(formulaOut=NULL, objectOut=NULL,formulaExp=form2,
                      objectExp= NULL,data=Compaq,idVar="id",effectOR =1,
                      MMOut=1,MMExp=c(1,1,1,1,1,1,1), nsim=2,
                      size=nrow(Compaq), eventRate=NULL, exposedPrev=NULL)


###################################################################################################
## One can provide the fitted model for the outcome model and the exposure model estimated by
## glm, gam, and bayesglm. The functional form of the fitted model for the outcome variable should
## of the form Outcome ~ Exposure + Confounders. The functional form of the exposure model is,
## Exposure ~ Confounders.
####################################################################################################

Coeff1<- bayesglm(form1, family = "binomial", data=Compaq,control=glm.control(trace=TRUE))
Coeff2<- bayesglm(form2, family = "binomial", data=Compaq,control=glm.control(trace=TRUE))
sizesim<-nrow(model.matrix(Coeff1))
sizesim1<-nrow(model.matrix(Coeff2))

Bin_Obj1<-PlasmodeBin(formulaOut=NULL, objectOut=Coeff1,formulaExp=NULL,
                     objectExp = Coeff2, idVar=Compaq$id,effectOR =1,
                     MMOut=c(1.5,1,2,1,1,1,1,1),MMExp=c(1,1,1,1,1,1,1),
                     nsim=2, size=sizesim, eventRate=NULL, exposedPrev=NULL)

Bin_Obj2<-PlasmodeBin(formulaOut=NULL, objectOut=Coeff1,formulaExp=NULL,
                     objectExp = NULL,idVar=Compaq$id,effectOR =1,
                     MMOut=c(1.5,1,2,1,1,1,1,1),MMExp=1,
                     nsim=2, size=sizesim, eventRate=NULL, exposedPrev=NULL)

Bin_Obj3<-PlasmodeBin(formulaOut=NULL, objectOut=NULL,formulaExp=NULL,
                     objectExp = Coeff2,idVar=Compaq$id,effectOR =1, MMOut=1,
                     MMExp=c(1,1,1,1,1,1,1),
                     nsim=2, size=sizesim1, eventRate=NULL, exposedPrev=NULL)
}



cleanEx()
nameEx("PlasmodeCont")
### * PlasmodeCont

flush(stderr()); flush(stdout())

### Name: PlasmodeCont
### Title: Performs the plasmode simulation
### Aliases: PlasmodeCont

### ** Examples

{
## Example for using the PlasmodeCont
library(twang)
library(gbm)
library(lattice)
library(parallel)
library(survey)
library(grid)
library(Matrix)
library(xtable)
library(latticeExtra)
library(RColorBrewer)
library(arm)
set.seed(1)
data("lalonde")
## Creating the ID variable
lalonde$id <- 1:nrow(lalonde)

str(lalonde)
## Example for PlasmodeCont when the outcome and exposure models formulas are provided.
form1<- re78 ~ treat + age + educ + black + hisp+ nodegr  + married + re74 + re75
form2<- treat ~ age + educ + black + hisp + nodegr + married + re74 + re75
Cont_Form1<-PlasmodeCont(formulaOut=form1, objectOut = NULL,formulaExp=form2,objectExp = NULL,
                        data=lalonde,idVar="id",effectOR =0, MMOut=c(0,1,2,1,1,1,2,2,1),
                        MMExp=c(1,2,1,1,1,2,2,1),nsim=2, size=nrow(lalonde),
                        eventRate=NULL, exposedPrev=NULL)
Cont_Form2<-PlasmodeCont(formulaOut=form1, objectOut = NULL,formulaExp=NULL,objectExp = NULL,
                        data=lalonde,idVar="id",effectOR =0, MMOut=c(0,1,2,1,1,1,2,2,1),MMExp=1,
                        nsim=2, size=nrow(lalonde), eventRate=NULL, exposedPrev=NULL)
Cont_Form3<-PlasmodeCont(formulaOut=NULL, objectOut = NULL,formulaExp=form2,objectExp = NULL,
                        data=lalonde,idVar="id",effectOR =0, MMOut=1,MMExp=c(1,2,1,1,1,2,2,1),
                        nsim=2, size=nrow(lalonde), eventRate=NULL, exposedPrev=NULL)
## Example for PlasmodeCont when the fitted model objects are provided.
###################################################################################################
## One can provide the fitted model for the outcome model and the exposure model estimated by
## glm, gam, and bayesglm. The functional form of the fitted model for the outcome variable should
## of the form Outcome ~ Exposure + Confounders. The functional form of the exposure model is,
## Exposure ~ Confounders.
####################################################################################################
Coeff1c<- bayesglm(form1, family = "gaussian", data=lalonde,control=glm.control(trace=TRUE))
Coeff2c<- bayesglm(form2, family = "binomial", data=lalonde,control=glm.control(trace=TRUE))

sizesim<-nrow(model.matrix(Coeff1c))
sizesim1<-nrow(model.matrix(Coeff2c))

Cont_Obj1<-PlasmodeCont(formulaOut=NULL, objectOut = Coeff1c,formulaExp=NULL,objectExp = Coeff2c,
                       idVar=lalonde$id,effectOR =0, MMOut=c(0,1,2,1,1,1,2,2,1),
                       MMExp=c(1,2,1,1,1,2,2,1),
                       nsim=2, size=nrow(lalonde), eventRate=NULL, exposedPrev=NULL)

Cont_Obj2<-PlasmodeCont(formulaOut=NULL, objectOut = Coeff1c,formulaExp=NULL,objectExp = NULL,
                       idVar=lalonde$id,effectOR =1, MMOut=c(0,1,2,1,1,1,2,2,1),MMExp=1,
                       nsim=2, size=nrow(lalonde), eventRate=NULL, exposedPrev=NULL)

Cont_Obj3<-PlasmodeCont(formulaOut=NULL, objectOut = NULL,formulaExp=NULL,objectExp = Coeff2c,
                       idVar=lalonde$id,effectOR =1, MMOut=c(0,1,2,1,1,1,2,2,1),MMExp=1,
                       nsim=2, size=nrow(lalonde), eventRate=NULL, exposedPrev=NULL)
}



cleanEx()
nameEx("PlasmodeSur")
### * PlasmodeSur

flush(stderr()); flush(stdout())

### Name: PlasmodeSur
### Title: Performs the plasmode simulation
### Aliases: PlasmodeSur

### ** Examples

{
library(survival)
library(splines)
library(glm2)
## Creating data set for simulation
lung <- lung[complete.cases(lung),]
lung$id <- 1:nrow(lung)
lung$meal.cal <- ifelse(lung$meal.cal > 1000, 1, 0)
lung$status <- lung$status - 1

## Formulas for estimating the hazard of outcome event, the hazard of censoring and exposure.

form1<-Surv(lung$time, lung$status)~meal.cal+age+sex+ph.ecog+ph.karno
form2<-Surv(lung$time, !lung$status)~meal.cal+age+sex+ph.ecog+ph.karno
form3<- meal.cal~age+sex+ph.ecog+ph.karno

Sur_Form1<-PlasmodeSur(formulaOut=form1,formulaCen=form2, objectOut=NULL, objectCen = NULL,
            formulaExp=form3,objectExp=NULL,data=lung,idVar="id",effectOR =1, MMOut=c(0.5,2,2,1,3),
            MMExp=c(2,2,2,2), nsim=3, size=nrow(lung), eventRate=NULL, exposedPrev=NULL)

Sur_Form2<-PlasmodeSur(formulaOut=form1,formulaCen=form2, objectOut=NULL, objectCen = NULL,
            formulaExp=NULL,objectExp=NULL,data=lung,idVar="id",effectOR =1, MMOut=c(1,2,2,1,3),
            MMExp=c(1,1,1,1),nsim=3, size=nrow(lung), eventRate=NULL, exposedPrev=NULL)

Sur_Form3<-PlasmodeSur(formulaOut=NULL,formulaCen=NULL, objectOut=NULL, objectCen = NULL,
            formulaExp=form3,objectExp=NULL,data=lung,idVar="id",effectOR =1, MMOut=c(1,2,2,1,3),
            MMExp=c(1,1,1,1),nsim=3, size=nrow(lung), eventRate=NULL, exposedPrev=NULL)

## Objects for the hazard of the outcome event, hazard for censoring and the exposure.

smod1 <- coxph(Surv(lung$time, lung$status)~meal.cal+age+sex+ph.ecog+ph.karno, data = lung,x=TRUE)
smod2 <- coxph(Surv(lung$time, !lung$status)~meal.cal+age+sex+ph.ecog+ph.karno, data = lung,x=TRUE)
pmod1<-glm2(meal.cal~age+sex+ph.ecog+ph.karno, data = lung,family = "binomial",
            control=glm.control(trace=TRUE))

Sur_Obj1<-PlasmodeSur(formulaOut=NULL,formulaCen=NULL, objectOut=smod1,objectCen = smod2,
            formulaExp=NULL,objectExp=pmod1,idVar=lung$id, effectOR =1, MMOut=c(1,2,2,1,3),
            MMExp=1, nsim=3,size=nrow(lung), eventRate=0.5, exposedPrev=NULL)

Sur_Obj2<-PlasmodeSur(formulaOut=NULL,formulaCen=NULL, objectOut=smod1,objectCen = smod2,
            formulaExp=NULL,objectExp=NULL,idVar=lung$id, effectOR =1.5, MMOut=c(1,2,2,1,3),
            MMExp=1, nsim=3,size=nrow(lung), eventRate=0.5, exposedPrev=NULL)

Sur_Obj3<-PlasmodeSur(formulaOut=NULL,formulaCen=NULL, objectOut=NULL,objectCen = NULL,
            formulaExp=NULL,objectExp=pmod1,idVar=lung$id,effectOR =1, MMOut=c(1,2,2,1,3),
            MMExp=1, nsim=3,size=nrow(lung), eventRate=0.5, exposedPrev=NULL)
}



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
