pkgname <- "corrgrapher"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('corrgrapher')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("calculate_cors")
### * calculate_cors

flush(stderr()); flush(stdout())

### Name: calculate_cors
### Title: Calculate correlation coefficients
### Aliases: calculate_cors calculate_cors.explainer calculate_cors.matrix
###   calculate_cors.table calculate_cors.default

### ** Examples


data(mtcars)
# Make sure, that categorical variables are factors
mtcars$vs <- factor(mtcars$vs, labels = c('V-shaped', 'straight'))
mtcars$am <- factor(mtcars$am, labels = c('automatic', 'manual'))
calculate_cors(mtcars)

# For a table:
data(HairEyeColor)
calculate_cors(HairEyeColor)

# Custom functions:
num_mtcars <- mtcars[,-which(colnames(mtcars) %in% c('vs', 'am'))]
my_f <- function(x,y) cor.test(x, y, method = 'spearman', exact=FALSE)$estimate
calculate_cors(num_mtcars, num_num_f = my_f, max_cor = 1)




cleanEx()
nameEx("corrgrapher")
### * corrgrapher

flush(stderr()); flush(stdout())

### Name: corrgrapher
### Title: Create a 'corrgrapher' object
### Aliases: corrgrapher corrgrapher.explainer corrgrapher.matrix
###   corrgrapher.default

### ** Examples

# convert the category variable
df <- as.data.frame(datasets::Seatbelts)
df$law <- factor(df$law) 
cgr <- corrgrapher(df)



cleanEx()
nameEx("plot.corrgrapher")
### * plot.corrgrapher

flush(stderr()); flush(stdout())

### Name: plot.corrgrapher
### Title: Visualize correlations in a corrgrapher object
### Aliases: plot.corrgrapher

### ** Examples

df <- as.data.frame(datasets::Seatbelts)[,1:7] # drop the binary target variable
cgr <- corrgrapher(df)
plot(cgr)



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
