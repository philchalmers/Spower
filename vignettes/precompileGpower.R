# Vignettes that take too long are precompiled

library(knitr)
setwd('vignettes')
rmarkdown::render("original/gpower_examples.Rmd", output_file="gpower_examples.Rmd")
