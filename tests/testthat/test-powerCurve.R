context('powerCurve')

test_that('powerCurve', {

	expect_class <- function(x, class) expect_true(inherits(x, class))
    library(Spower)

	p_t.test(n=NA, d=.2) |>
		powerCurve(n=c(30, 90, 270), verbose=FALSE, replications=10, plotly = FALSE) -> out
	expect_class(out, 'gg')

	p_t.test(d=.2) |>
		powerCurve(n=c(30, 90, 270), verbose=FALSE, replications=10, plotly = FALSE) -> out
	expect_class(out, 'gg')

	p_t.test(n=NA, d=0.2) |>
	   powerCurve(power=c(.2, .4), interval=c(10, 1000), verbose=F) -> out
	expect_class(out, 'gg')

	p_t.test() |>
		powerCurve(n=c(30, 90, 270, 550),
				   d=c(.2, .5, .8), replications=10, verbose=F) -> out
	expect_class(out, 'gg')

	p_t.test(n=NA) |>
		powerCurve(d=c(.2, .5, .8), power=.80, interval=c(10, 1000), verbose=F) -> out
	expect_class(out, 'gg')



})


