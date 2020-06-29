pkgname <- "crispRdesignR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('crispRdesignR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("crispRdesignRUI")
### * crispRdesignRUI

flush(stderr()); flush(stdout())

### Name: crispRdesignRUI
### Title: UI caller for crispRdesignR
### Aliases: crispRdesignRUI

### ** Examples

requireNamespace("gbm", quietly = TRUE)
requireNamespace("Biostrings", quietly = TRUE)
if (interactive()) {
  crispRdesignRUI()
  }



cleanEx()
nameEx("getofftargetdata")
### * getofftargetdata

flush(stderr()); flush(stdout())

### Name: getofftargetdata
### Title: Off Target Data Frame Creation
### Aliases: getofftargetdata

### ** Examples


## Quick example without off-target searching or annotation
## First generate data with the sgRNA_Design Function
testseq <- "GGCAGAGCTTCGTATGTCGGCGATTCATCTCAAGTAGAAGATCCTGGTGCAGTAGG"
usergenome <- "placeholder"
gtfname <- "placeholder"
alldata <- sgRNA_design(testseq, usergenome, gtfname, calloffs = FALSE)
## Then separate and format the off-target data with getofftargetdata()
final_data <- getofftargetdata(alldata)




cleanEx()
nameEx("getsgRNAdata")
### * getsgRNAdata

flush(stderr()); flush(stdout())

### Name: getsgRNAdata
### Title: sgRNA Data Frame Creation
### Aliases: getsgRNAdata

### ** Examples


## Quick example without off-target searching or annotation
## First generate data with the sgRNA_Design Function
testseq <- "GGCAGAGCTTCGTATGTCGGCGATTCATCTCAAGTAGAAGATCCTGGTGCAGTAGG"
usergenome <- "placeholder"
gtfname <- "placeholder"
alldata <- sgRNA_design(testseq, usergenome, gtfname, calloffs = FALSE)
## Then separate and format the sgRNA data with getsgRNAdata()
final_data <- getsgRNAdata(alldata)




cleanEx()
nameEx("sgRNA_design")
### * sgRNA_design

flush(stderr()); flush(stdout())

### Name: sgRNA_design
### Title: sgRNA Target Design
### Aliases: sgRNA_design

### ** Examples

## Quick example without off-target searching or annotation
testseq <- "GGCAGAGCTTCGTATGTCGGCGATTCATCTCAAGTAGAAGATCCTGGTGCAGTAGG"
usergenome <- "placeholder"
gtfname <- "placeholder"
alldata <- sgRNA_design(testseq, usergenome, gtfname, calloffs = FALSE)




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
