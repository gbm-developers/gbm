pkgname <- "twang"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('twang')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("raceprofiling")
### * raceprofiling

flush(stderr()); flush(stdout())

### Name: raceprofiling
### Title: Traffic stop data
### Aliases: raceprofiling
### Keywords: datasets

### ** Examples

data(raceprofiling)

# the first five lines of the dataset
raceprofiling[1:5,]



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
