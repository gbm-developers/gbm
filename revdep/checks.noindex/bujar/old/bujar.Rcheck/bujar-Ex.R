pkgname <- "bujar"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('bujar')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bujar")
### * bujar

flush(stderr()); flush(stdout())

### Name: bujar
### Title: Buckley-James Regression
### Aliases: bujar print.bujar plot.bujar summary.bujar

### ** Examples

data("wpbc", package = "TH.data")
wpbc2 <- wpbc[, 1:12]
wpbc2$status <- as.numeric(wpbc2$status) - 1
fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x= wpbc2[, -(1:2)])
print(fit)
coef(fit)
pr <- predict(fit)
plot(fit)
fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x= wpbc2[, -(1:2)], tuning = TRUE)
## Not run: 
##D fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x=wpbc2[, -(1:2)], learner="pspline")
##D fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x=wpbc2[, -(1:2)], 
##D  learner="tree", degree=2)
##D ### select tuning parameter for "enet"
##D tmp <- gcv.enet(y=log(wpbc2$time), cens=wpbc2$status, x=wpbc2[, -(1:2)])
##D fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x=wpbc2[, -(1:2)], learner="enet", 
##D lamb = tmp$lambda, s=tmp$s)
##D 
##D fit <- bujar(y=log(wpbc2$time),cens=wpbc2$status, x=wpbc2[, -(1:2)], learner="mars", 
##D degree=2)
##D summary(fit)
## End(Not run)



cleanEx()
nameEx("chop")
### * chop

flush(stderr()); flush(stdout())

### Name: chop
### Title: Survival of CHOP for diffuse large B cell lymphoma
### Aliases: chop
### Keywords: datasets

### ** Examples

data(chop)
str(chop)



cleanEx()
nameEx("rchop")
### * rchop

flush(stderr()); flush(stdout())

### Name: rchop
### Title: Survival of R-CHOP for diffuse large B cell lymphoma
### Aliases: rchop
### Keywords: datasets

### ** Examples

data(rchop)
str(rchop)



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
