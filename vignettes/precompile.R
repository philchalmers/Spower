# Vignettes that take too long are precompiled

library(knitr)
setwd('vignettes')
knit("original/gpower_examples.Rmd", "gpower_examples.Rmd")
