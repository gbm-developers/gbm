pkgname <- "mma"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mma')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("boot.med")
### * boot.med

flush(stderr()); flush(stdout())

### Name: boot.med
### Title: Statistical Inference on Mediation Analysis with Continuous or
###   Binary Predictor
### Aliases: boot.med
### Keywords: Mediation Analysis Continuous Predictor

### ** Examples

data("weight_behavior")
##binary x
#binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
 temp1<-boot.med(data=data.bin,n=2,n2=4)
 temp2<-boot.med(data=data.bin,n=2,n2=4,nu=0.05,nonlinear=TRUE)


##Surv class outcome (survival analysis)
## End(No test)



cleanEx()
nameEx("boot.mod")
### * boot.mod

flush(stderr()); flush(stdout())

### Name: boot.mod
### Title: Statistical Inference on Mediation Analysis with Continuous or
###   Binary Predictor at different level of the moderator
### Aliases: boot.mod
### Keywords: Moderation Effect with Mediation Analysis

### ** Examples




cleanEx()
nameEx("cgd1")
### * cgd1

flush(stderr()); flush(stdout())

### Name: cgd1
### Title: cgd1 Data Set
### Aliases: cgd1
### Keywords: Datasets

### ** Examples

data(cgd1)
names(cgd1)



cleanEx()
nameEx("data.org")
### * data.org

flush(stderr()); flush(stdout())

### Name: data.org
### Title: Data Organization and Identify Potential Mediators
### Aliases: data.org
### Keywords: Mediator Tests

### ** Examples

data("weight_behavior")
#binary predictor
 #binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.b.b.2.1<-data.org(x,y,mediator=5:12,jointm=list(n=1,j1=c(5,7,9)),
                        pred=pred,predref="M", alpha=0.4,alpha2=0.4)
 summary(data.b.b.2.1)
 #Or you can specify the potential mediators and change the reference 
 #group for binary or categorical mediators. In the following code,
 #potential continuous mediators are columns 8,9,10,12, and 13 of x,
 #binary mediators are columns 7 and 11, and categorical mediator is
 #column 6 of x with 1 to be the reference group for all categorical
 #and binary mediators. 
  data.b.b.2<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
   binref=c(1,1),catmed=5,catref=1,jointm=list(n=1,j1=c(5,7,9)),
   predref="M",alpha=0.4,alpha2=0.4) 
  summary(data.b.b.2)
 #use the mediator argument instead of contmet, binmed and catmed
 
 #multivariate predictor



cleanEx()
nameEx("form.interaction")
### * form.interaction

flush(stderr()); flush(stdout())

### Name: form.interaction
### Title: Create interaction terms of predictor(s) and potential
###   moderator(s).
### Aliases: form.interaction

### ** Examples

data("weight_behavior")
pred=data.frame(weight_behavior[,3])
names(pred)="pred"
x=weight_behavior[,c(2,4:14)]
inter=form.interaction(x,pred,inter.cov=c("sports","sweat"),predref=NULL) 
x=cbind(x,inter)
head(x)



cleanEx()
nameEx("med")
### * med

flush(stderr()); flush(stdout())

### Name: med
### Title: Mediation Analysis with Binary or Continuous Predictor
### Aliases: med
### Keywords: Mediation Analysis

### ** Examples

data("weight_behavior")
##binary x
#binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
  binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp1<-med(data=data.bin,n=2)
#or use self-defined final function
temp1<-med(data=data.bin,n=2,custom.function = 
           'glm(responseY~.,data=dataset123,family="quasibinomial",
           weights=weights123)')
temp2<-med(data=data.bin,n=2,nonlinear=TRUE)




cleanEx()
nameEx("mma-package")
### * mma-package

flush(stderr()); flush(stdout())

### Name: mma-package
### Title: Mediation Analysis Package
### Aliases: mma-package
### Keywords: Package

### ** Examples

data("weight_behavior")
#binary predictor
 #binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 temp.b.b.glm<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,n=2,n2=2)



cleanEx()
nameEx("mma")
### * mma

flush(stderr()); flush(stdout())

### Name: mma
### Title: Multiple Mediation Analysis
### Aliases: mma
### Keywords: Mediation Analysis Mediator Tests

### ** Examples

data("weight_behavior")
#binary predictor
 #binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 temp.b.b.glm<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,n=2,n2=2)




cleanEx()
nameEx("moderate")
### * moderate

flush(stderr()); flush(stdout())

### Name: moderate
### Title: Calculate and plot the direct effect of the selected exposure
###   variable at each level of the moderator.
### Aliases: moderate

### ** Examples




cleanEx()
nameEx("plot.med")
### * plot.med

flush(stderr()); flush(stdout())

### Name: plot.med
### Title: Plot the mediation effect on the fitted med object
### Aliases: plot.med

### ** Examples

data("weight_behavior")
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
  binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp1<-med(data=data.bin,n=2)
temp2<-med(data=data.bin,n=2,nonlinear=TRUE)
plot(temp1,data.bin,vari="exercises",xlim=c(0,50))
plot(temp2,data.bin,vari="sports")



cleanEx()
nameEx("plot.mma")
### * plot.mma

flush(stderr()); flush(stdout())

### Name: plot.mma
### Title: Relative effects plot of the fitted mma object
### Aliases: plot.mma

### ** Examples

data("weight_behavior")
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 temp.b.b.glm<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,n=2,n2=2)
plot(temp.b.b.glm,vari="exercises",xlim=c(0,50))
plot(temp.b.b.glm,vari="sports")



cleanEx()
nameEx("plot2.mma")
### * plot2.mma

flush(stderr()); flush(stdout())

### Name: plot2.mma
### Title: Relative effects plot of the fitted mma object with moderator
### Aliases: plot2.mma

### ** Examples

#see boot.mod menu.



cleanEx()
nameEx("print.med")
### * print.med

flush(stderr()); flush(stdout())

### Name: print.med
### Title: Print an med object
### Aliases: print.med

### ** Examples

data("weight_behavior")
##binary x
#binary y
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                    binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
 temp1<-med(data=data.bin,n=2)
 temp2<-med(data=data.bin,n=2,nonlinear=TRUE)
 temp1
 print(temp2,digit=5)



cleanEx()
nameEx("print.mma")
### * print.mma

flush(stderr()); flush(stdout())

### Name: print.mma
### Title: Print a mma object
### Aliases: print.mma

### ** Examples

 data("weight_behavior")
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 temp.b.b.glm<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,n=2,n2=2)
 print(temp.b.b.glm,digit=8)



cleanEx()
nameEx("summary.med_iden")
### * summary.med_iden

flush(stderr()); flush(stdout())

### Name: summary.med_iden
### Title: Summary method for class "med_iden".
### Aliases: summary.med_iden print.summary.med_iden

### ** Examples

data("weight_behavior")
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
 data.b.b.2<-data.org(x,y,mediator=5:12,jointm=list(n=1,j1=c(5,7,9)),
                        pred=pred,predref="M", alpha=0.4,alpha2=0.4)
 summary(data.b.b.2)



cleanEx()
nameEx("summary.mma")
### * summary.mma

flush(stderr()); flush(stdout())

### Name: summary.mma
### Title: Summary of an mma project
### Aliases: summary.mma print.summary.mma

### ** Examples

data("weight_behavior")
 x=weight_behavior[,c(2,4:14)]
 pred=weight_behavior[,3]
 y=weight_behavior[,15]
  temp.b.b.glm<-mma(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),binref=c(1,1),
                    catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4,n=2,n2=2)
 summary(temp.b.b.glm, RE=TRUE, ball.use=FALSE)
 summary(temp.b.b.glm, ball.use=FALSE)



cleanEx()
nameEx("test.moderation")
### * test.moderation

flush(stderr()); flush(stdout())

### Name: test.moderation
### Title: Test for moderation effects.
### Aliases: test.moderation

### ** Examples

data("weight_behavior")
x=weight_behavior[,c(2,4:14)]
pred=weight_behavior[,3]
y=weight_behavior[,15]
data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                   binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp2<-med(data=data.bin,n=2,nonlinear=TRUE)
test.moderation(temp2,c("sports","sweat"),j=1,kx=NULL)

x=cbind(x,form.interaction(x,pred,inter.cov=c("sports","sweat"),predref=NULL))

data.bin<-data.org(x,y,pred=pred,contmed=c(7:9,11:12),binmed=c(6,10),
                   binref=c(1,1),catmed=5,catref=1,predref="M",alpha=0.4,alpha2=0.4)
temp1<-med(data=data.bin,n=2)
test.moderation(temp1,c("sports","sweat"),j=1,kx=NULL)



cleanEx()
nameEx("weight_behavior")
### * weight_behavior

flush(stderr()); flush(stdout())

### Name: weight_behavior
### Title: Weight_Behavior Data Set
### Aliases: weight_behavior
### Keywords: Datasets

### ** Examples

data(weight_behavior)
names(weight_behavior)



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
