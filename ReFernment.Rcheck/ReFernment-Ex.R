pkgname <- "ReFernment"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ReFernment')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ReFernment")
### * ReFernment

flush(stderr()); flush(stdout())

### Name: ReFernment
### Title: Annotate RNA editing in the GB files of plastomes
### Aliases: ReFernment

### ** Examples

genomes <- c("Asplenium_pek", "Woodwardia_uni")
gbFolder <- "../examples/GB/"
gffFolder <- "../examples/GFF/"
outputFolder <- "../examples"
ReFernment(gbFolder, gffFolder, fastaFolder, outputFolder, genomes)




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
