### R code from vignette source 'replicated.Rnw'

###################################################
### code chunk number 1: replicated.Rnw:29-30
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: replicated.Rnw:35-42
###################################################
library(spatstat)
spatstat.options(image.colfun=function(n) { grey(seq(0,1,length=n)) })
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 3: replicated.Rnw:180-181
###################################################
waterstriders


###################################################
### code chunk number 4: replicated.Rnw:199-200
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(waterstriders, main="")


###################################################
### code chunk number 5: replicated.Rnw:207-208
###################################################
summary(waterstriders)


###################################################
### code chunk number 6: replicated.Rnw:216-217
###################################################
X <- listof(rpoispp(100), rpoispp(100), rpoispp(100))


###################################################
### code chunk number 7: replicated.Rnw:222-224
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(X)
X


###################################################
### code chunk number 8: replicated.Rnw:253-254 (eval = FALSE)
###################################################
## hyperframe(...)


###################################################
### code chunk number 9: replicated.Rnw:279-281
###################################################
H <- hyperframe(X=1:3, Y=list(sin,cos,tan))
H


###################################################
### code chunk number 10: replicated.Rnw:289-294
###################################################
G <- hyperframe(X=1:3, Y=letters[1:3], Z=factor(letters[1:3]),
                W=list(rpoispp(100),rpoispp(100), rpoispp(100)),
                U=42,
                V=rpoispp(100), stringsAsFactors=FALSE)
G


###################################################
### code chunk number 11: replicated.Rnw:324-325
###################################################
simba


###################################################
### code chunk number 12: replicated.Rnw:338-339
###################################################
pyramidal


###################################################
### code chunk number 13: replicated.Rnw:345-346
###################################################
ws <- hyperframe(Striders=waterstriders)


###################################################
### code chunk number 14: replicated.Rnw:353-355
###################################################
H$X
H$Y


###################################################
### code chunk number 15: replicated.Rnw:365-367
###################################################
H$U <- letters[1:3]
H


###################################################
### code chunk number 16: replicated.Rnw:372-376
###################################################
G <- hyperframe()
G$X <- waterstriders
G$Y <- 1:3
G


###################################################
### code chunk number 17: replicated.Rnw:384-388
###################################################
H[,1]
H[2,]
H[2:3, ]
H[1,1]


###################################################
### code chunk number 18: replicated.Rnw:394-397
###################################################
H[,1,drop=TRUE]
H[1,1,drop=TRUE]
H[1,2,drop=TRUE]


###################################################
### code chunk number 19: replicated.Rnw:410-411 (eval = FALSE)
###################################################
## plot.listof(x, ..., main, arrange = TRUE, nrows = NULL, ncols = NULL)


###################################################
### code chunk number 20: replicated.Rnw:426-427
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(waterstriders, pch=16, nrows=1)


###################################################
### code chunk number 21: replicated.Rnw:442-443
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(simba)


###################################################
### code chunk number 22: replicated.Rnw:455-457
###################################################
getOption("SweaveHooks")[["fig"]]()
H <- hyperframe(X=1:3, Y=list(sin,cos,tan))
plot(H$Y)


###################################################
### code chunk number 23: replicated.Rnw:469-470 (eval = FALSE)
###################################################
## plot(h, e)


###################################################
### code chunk number 24: replicated.Rnw:479-480
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(demohyper, quote({ plot(Image, main=""); plot(Points, add=TRUE) }))


###################################################
### code chunk number 25: replicated.Rnw:492-494
###################################################
getOption("SweaveHooks")[["fig"]]()
H <- hyperframe(Bugs=waterstriders)
plot(H, quote(plot(Kest(Bugs))), marsize=1)


###################################################
### code chunk number 26: replicated.Rnw:507-509
###################################################
df <- data.frame(A=1:10, B=10:1)
with(df, A-B)


###################################################
### code chunk number 27: replicated.Rnw:522-523 (eval = FALSE)
###################################################
## with(h,e)


###################################################
### code chunk number 28: replicated.Rnw:533-536
###################################################
H <- hyperframe(Bugs=waterstriders)
with(H, npoints(Bugs))
with(H, distmap(Bugs))


###################################################
### code chunk number 29: replicated.Rnw:559-560
###################################################
with(simba, npoints(Points))


###################################################
### code chunk number 30: replicated.Rnw:567-569
###################################################
H <- hyperframe(Bugs=waterstriders)
K <- with(H, Kest(Bugs))


###################################################
### code chunk number 31: replicated.Rnw:577-578
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(K)


###################################################
### code chunk number 32: replicated.Rnw:583-585
###################################################
H <- hyperframe(Bugs=waterstriders)
with(H, nndist(Bugs))


###################################################
### code chunk number 33: replicated.Rnw:591-592
###################################################
with(H, min(nndist(Bugs)))


###################################################
### code chunk number 34: replicated.Rnw:604-605
###################################################
simba$Dist <- with(simba, distmap(Points))


###################################################
### code chunk number 35: replicated.Rnw:618-622
###################################################
getOption("SweaveHooks")[["fig"]]()
lambda <- rexp(6, rate=1/50)
H <- hyperframe(lambda=lambda)
H$Points <- with(H, rpoispp(lambda))
plot(H, quote(plot(Points, main=paste("lambda=", signif(lambda, 4)))))


###################################################
### code chunk number 36: replicated.Rnw:628-629
###################################################
H$X <- with(H, rpoispp(50))


###################################################
### code chunk number 37: replicated.Rnw:658-659
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(simba, quote(plot(density(Points), main="")), nrows=2)


###################################################
### code chunk number 38: replicated.Rnw:678-680
###################################################
getOption("SweaveHooks")[["fig"]]()
rhos <- with(demohyper, rhohat(Points, Image))
plot(rhos)


###################################################
### code chunk number 39: replicated.Rnw:698-699 (eval = FALSE)
###################################################
## mppm(formula, data, interaction, ...)


###################################################
### code chunk number 40: replicated.Rnw:709-710 (eval = FALSE)
###################################################
## mppm(Points ~ group, simba, Poisson())


###################################################
### code chunk number 41: replicated.Rnw:743-744
###################################################
mppm(Points ~ 1, simba)


###################################################
### code chunk number 42: replicated.Rnw:751-752
###################################################
mppm(Points ~ group, simba)


###################################################
### code chunk number 43: replicated.Rnw:758-759
###################################################
mppm(Points ~ id, simba)


###################################################
### code chunk number 44: replicated.Rnw:769-770
###################################################
mppm(Points ~ Image, data=demohyper)


###################################################
### code chunk number 45: replicated.Rnw:788-789 (eval = FALSE)
###################################################
## mppm(Points ~ offset(log(Image)), data=demohyper)


###################################################
### code chunk number 46: replicated.Rnw:801-802 (eval = FALSE)
###################################################
## mppm(Points ~ log(Image), data=demop)


###################################################
### code chunk number 47: replicated.Rnw:819-820 (eval = FALSE)
###################################################
## mppm(formula, data, interaction, ..., iformula=NULL)


###################################################
### code chunk number 48: replicated.Rnw:870-871
###################################################
radii <- with(simba, mean(nndist(Points)))


###################################################
### code chunk number 49: replicated.Rnw:878-880
###################################################
Rad <- hyperframe(R=radii)
Str <- with(Rad, Strauss(R))


###################################################
### code chunk number 50: replicated.Rnw:885-887
###################################################
Int <- hyperframe(str=Str)
mppm(Points ~ 1, simba, interaction=Int)


###################################################
### code chunk number 51: replicated.Rnw:914-917
###################################################
h <- hyperframe(Y=waterstriders)
g <- hyperframe(po=Poisson(), str4 = Strauss(4), str7= Strauss(7))
mppm(Y ~ 1, data=h, interaction=g, iformula=~str4)


###################################################
### code chunk number 52: replicated.Rnw:928-929
###################################################
fit <- mppm(Points ~ 1, simba, Strauss(0.07), iformula = ~Interaction*group)


###################################################
### code chunk number 53: replicated.Rnw:947-948
###################################################
fit


###################################################
### code chunk number 54: replicated.Rnw:951-953
###################################################
co <- coef(fit)
si <- function(x) { signif(x, 4) }


###################################################
### code chunk number 55: replicated.Rnw:964-965
###################################################
coef(fit)


###################################################
### code chunk number 56: replicated.Rnw:1022-1023 (eval = FALSE)
###################################################
## interaction=hyperframe(po=Poisson(), str=Strauss(0.07))


###################################################
### code chunk number 57: replicated.Rnw:1028-1029 (eval = FALSE)
###################################################
## iformula=~ifelse(group=="control", po, str)


###################################################
### code chunk number 58: replicated.Rnw:1039-1040 (eval = FALSE)
###################################################
## iformula=~I((group=="control")*po) + I((group=="treatment") * str)


###################################################
### code chunk number 59: replicated.Rnw:1050-1055
###################################################
g <- hyperframe(po=Poisson(), str=Strauss(0.07))
fit2 <- mppm(Points ~ 1, simba, g, 
             iformula=~I((group=="control")*po) 
                     + I((group=="treatment") * str))
fit2


###################################################
### code chunk number 60: replicated.Rnw:1178-1180
###################################################
H <- hyperframe(P=waterstriders)
mppm(P ~ 1, H, random=~1|id)


###################################################
### code chunk number 61: replicated.Rnw:1187-1188 (eval = FALSE)
###################################################
## mppm(Neurons ~ AstroIm, random=~AstroIm|WellNumber)


###################################################
### code chunk number 62: replicated.Rnw:1211-1214
###################################################
H <- hyperframe(W=waterstriders)
fit <- mppm(W ~ 1, H)
subfits(fit)


###################################################
### code chunk number 63: replicated.Rnw:1235-1236 (eval = FALSE)
###################################################
## subfits <- subfits.new


###################################################
### code chunk number 64: replicated.Rnw:1248-1250
###################################################
H <- hyperframe(W=waterstriders)
with(H, ppm(W))


###################################################
### code chunk number 65: replicated.Rnw:1273-1275
###################################################
fit <- mppm(P ~ x, hyperframe(P=waterstriders))
res <- residuals(fit)


###################################################
### code chunk number 66: replicated.Rnw:1285-1286
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(res)


###################################################
### code chunk number 67: replicated.Rnw:1291-1293
###################################################
getOption("SweaveHooks")[["fig"]]()
smor <- with(hyperframe(res=res), Smooth(res, sigma=4))
plot(smor)


###################################################
### code chunk number 68: replicated.Rnw:1305-1308
###################################################
fit <- mppm(P ~ x, hyperframe(P=waterstriders))
res <- residuals(fit)
totres <- sapply(res, integral.msr)


###################################################
### code chunk number 69: replicated.Rnw:1314-1321
###################################################
getOption("SweaveHooks")[["fig"]]()
fit <- mppm(Points~Image, data=demohyper)
resids <- residuals(fit, type="Pearson")
totres <- sapply(resids, integral.msr)
areas <- with(demohyper, area.owin(as.owin(Points)))
df <- as.data.frame(demohyper[, "Group"])
df$resids <- totres/areas
plot(resids~Group, df)


###################################################
### code chunk number 70: replicated.Rnw:1342-1345
###################################################
getOption("SweaveHooks")[["fig"]]()
fit <- mppm(P ~ 1, hyperframe(P=waterstriders))
sub <- hyperframe(Model=subfits(fit))
plot(sub, quote(diagnose.ppm(Model)))


###################################################
### code chunk number 71: replicated.Rnw:1358-1366
###################################################
H <- hyperframe(P = waterstriders)
fitall <- mppm(P ~ 1, H)
together <- subfits(fitall)
separate <- with(H, ppm(P))
Fits <- hyperframe(Together=together, Separate=separate)
dr <- with(Fits, unlist(coef(Separate)) - unlist(coef(Together)))
dr
exp(dr)


###################################################
### code chunk number 72: replicated.Rnw:1383-1392
###################################################
H <- hyperframe(X=waterstriders)

# Poisson with constant intensity for all patterns
fit1 <- mppm(X~1, H)
quadrat.test(fit1, nx=2)

# uniform Poisson with different intensity for each pattern
fit2 <- mppm(X ~ id, H)
quadrat.test(fit2, nx=2)


###################################################
### code chunk number 73: replicated.Rnw:1421-1422 (eval = FALSE)
###################################################
## kstest.mppm(model, covariate)


