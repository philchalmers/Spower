# Vignettes that take too long are precompiled

library(knitr)
setwd('vignettes')
knitr::knit("original/gpower_examples.Rmd", "gpower_examples.Rmd")
