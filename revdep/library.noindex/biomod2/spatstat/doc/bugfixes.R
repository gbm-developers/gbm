### R code from vignette source 'bugfixes.Rnw'

###################################################
### code chunk number 1: bugfixes.Rnw:20-24
###################################################
library(spatstat)
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 2: bugfixes.Rnw:37-41
###################################################
nbugs <- nrow(news(grepl("^BUG", Category), 
                   package="spatstat"))
nbugssince <- nrow(news(Version > "1.42-0" & grepl("^BUG", Category), 
                   package="spatstat"))


###################################################
### code chunk number 3: bugfixes.Rnw:57-58 (eval = FALSE)
###################################################
## bugfixes


###################################################
### code chunk number 4: bugfixes.Rnw:64-65 (eval = FALSE)
###################################################
## bugfixes(sincedate="2017-06-30")


###################################################
### code chunk number 5: bugfixes.Rnw:68-69 (eval = FALSE)
###################################################
## bugfixes(sinceversion="1.42-0")


