pkgname <- "SRGnet"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SRGnet')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("SRGnet")
### * SRGnet

flush(stderr()); flush(stdout())

### Name: SRGnet
### Title: Synergistic response to gene mutations specific network
### Aliases: SRGnet

### ** Examples

{
data(Differentially_expressed_genes)
data(Transcriptomics)
data(PLCRG)
SRGnet("F") #Fast run  
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
