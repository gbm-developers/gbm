pkgname <- "mob"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mob')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("arb_bin")
### * arb_bin

flush(stderr()); flush(stdout())

### Name: arb_bin
### Title: Monotonic binning based on the decision tree
### Aliases: arb_bin

### ** Examples

data(hmeq)
arb_bin(hmeq, BAD, DEROG)



cleanEx()
nameEx("bad_bin")
### * bad_bin

flush(stderr()); flush(stdout())

### Name: bad_bin
### Title: Monotonic binning based on bads only
### Aliases: bad_bin

### ** Examples

data(hmeq)
bad_bin(hmeq, BAD, DEROG)



cleanEx()
nameEx("cal_woe")
### * cal_woe

flush(stderr()); flush(stdout())

### Name: cal_woe
### Title: Perform WoE transformation of a numeric variable
### Aliases: cal_woe

### ** Examples

data(hmeq)
bin_out <- qtl_bin(hmeq, BAD, DEROG)
cal_woe(hmeq, "DEROG", bin_out$df)



cleanEx()
nameEx("gbm_bin")
### * gbm_bin

flush(stderr()); flush(stdout())

### Name: gbm_bin
### Title: Monotonic binning based on the generalized boosted regression
###   model (GBM)
### Aliases: gbm_bin

### ** Examples

data(hmeq)
gbm_bin(hmeq, BAD, DEROG)



cleanEx()
nameEx("iso_bin")
### * iso_bin

flush(stderr()); flush(stdout())

### Name: iso_bin
### Title: Monotonic binning based on isotonic regression
### Aliases: iso_bin

### ** Examples

data(hmeq)
iso_bin(hmeq, BAD, DEROG)



cleanEx()
nameEx("manual_bin")
### * manual_bin

flush(stderr()); flush(stdout())

### Name: manual_bin
### Title: Binning based on cut points
### Aliases: manual_bin

### ** Examples

data(hmeq)
manual_bin(hmeq, "BAD", "DEROG", c(0, 1, 2))



cleanEx()
nameEx("qtl_bin")
### * qtl_bin

flush(stderr()); flush(stdout())

### Name: qtl_bin
### Title: Monotonic binning by quantile
### Aliases: qtl_bin

### ** Examples

data(hmeq)
qtl_bin(hmeq, BAD, DEROG)



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
