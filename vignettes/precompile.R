# Vignettes that take too long are precompiled

library(knitr)
setwd('vignettes')
knit("original/gpower_examples.Rmd", "gpower_examples.Rmd")

Sys.setenv(SPOWER_EVAL = TRUE)
rmarkdown::render("SpowerIntro.Rmd")
rmarkdown::render("SpowerIntro_logicals.Rmd")
SimDesign::SimClean("SpowerIntro.html", "SpowerIntro_logicals.html")
