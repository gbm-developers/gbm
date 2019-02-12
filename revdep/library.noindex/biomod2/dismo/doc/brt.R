### R code from vignette source 'brt.Rnw'

###################################################
### code chunk number 1: foo
###################################################
options(width = 65)
set.seed(0)


###################################################
### code chunk number 2: brt-0
###################################################
library(dismo)
data(Anguilla_train)
head(Anguilla_train)


###################################################
### code chunk number 3: brt-1
###################################################
angaus.tc5.lr01 <- gbm.step(data=Anguilla_train, gbm.x = 3:13, gbm.y = 2,
                        family = "bernoulli", tree.complexity = 5,
                        learning.rate = 0.01, bag.fraction = 0.5)


###################################################
### code chunk number 4: brt-2
###################################################
names(angaus.tc5.lr01)
summary(angaus.tc5.lr01)


###################################################
### code chunk number 5: brt-3
###################################################
angaus.tc5.lr005 <- gbm.step(data=Anguilla_train, gbm.x = 3:13, gbm.y = 2, 
                      family = "bernoulli", tree.complexity = 5,
                      learning.rate = 0.005, bag.fraction = 0.5)


###################################################
### code chunk number 6: brt-4
###################################################
angaus.simp <- gbm.simplify(angaus.tc5.lr005, n.drops = 5)


###################################################
### code chunk number 7: brt-5
###################################################
angaus.tc5.lr005.simp <- gbm.step(Anguilla_train, 
                   gbm.x=angaus.simp$pred.list[[1]], gbm.y=2,
                   tree.complexity=5, learning.rate=0.005)


###################################################
### code chunk number 8: brt-6
###################################################
gbm.plot(angaus.tc5.lr005, n.plots=11, write.title = FALSE)


###################################################
### code chunk number 9: brt-7
###################################################
gbm.plot.fits(angaus.tc5.lr005)


###################################################
### code chunk number 10: dismo-8
###################################################
find.int <- gbm.interactions(angaus.tc5.lr005)
find.int$interactions
find.int$rank.list


###################################################
### code chunk number 11: brt-9
###################################################
gbm.perspec(angaus.tc5.lr005, 7, 1, y.range=c(15,20), z.range=c(0,0.6))


###################################################
### code chunk number 12: dismo-10
###################################################
data(Anguilla_test)
library(gbm)
preds <- predict.gbm(angaus.tc5.lr005, Anguilla_test,
         n.trees=angaus.tc5.lr005$gbm.call$best.trees, type="response")

calc.deviance(obs=Anguilla_test$Angaus_obs, pred=preds, calc.mean=TRUE)
d <- cbind(Anguilla_test$Angaus_obs, preds)
pres <- d[d[,1]==1, 2]
abs <- d[d[,1]==0, 2]
e <- evaluate(p=pres, a=abs)
e


###################################################
### code chunk number 13: dismo-11
###################################################
angaus.5000 <- gbm.fixed(data=Anguilla_train, gbm.x=3:13, gbm.y=2,
               learning.rate=0.005, tree.complexity=5, n.trees=5000)
tree.list <- seq(100, 5000, by=100)
pred <- predict.gbm(angaus.5000, Anguilla_test, n.trees=tree.list, "response")


###################################################
### code chunk number 14: dismo-12
###################################################
angaus.pred.deviance <- rep(0,50)
for (i in 1:50) {
   angaus.pred.deviance[i] <- calc.deviance(Anguilla_test$Angaus_obs,
                               pred[,i], calc.mean=TRUE)
}


###################################################
### code chunk number 15: brt-12
###################################################
plot(tree.list, angaus.pred.deviance, ylim=c(0.7,1), xlim=c(-100,5000),
     type='l', xlab="number of trees", ylab="predictive deviance",
     cex.lab=1.5) 


###################################################
### code chunk number 16: brt-13
###################################################
data(Anguilla_grids)
plot(Anguilla_grids)


###################################################
### code chunk number 17: brt14
###################################################
Method <- factor('electric', levels = levels(Anguilla_train$Method))
add <- data.frame(Method)
p <- predict(Anguilla_grids, angaus.tc5.lr005, const=add, 
       n.trees=angaus.tc5.lr005$gbm.call$best.trees, type="response")
p <- mask(p, raster(Anguilla_grids, 1))
plot(p, main='Angaus - BRT prediction')


