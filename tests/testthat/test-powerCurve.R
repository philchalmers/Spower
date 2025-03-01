context('SpowerCurve')

test_that('SpowerCurve', {

	expect_class <- function(x, class) expect_true(inherits(x, class))
    library(Spower)

	p_t.test(n=NA, d=.2) |>
		SpowerCurve(n=c(30, 90, 270), verbose=FALSE, replications=10, plotly = FALSE) -> out
	expect_class(out, 'gg')

	p_t.test(d=.2) |>
		SpowerCurve(n=c(30, 90, 270), verbose=FALSE, replications=10, plotly = FALSE) -> out
	expect_class(out, 'gg')

	p_t.test(n=NA, d=0.2) |>
	   SpowerCurve(power=c(.2, .4), interval=c(10, 1000), verbose=F, plotly = FALSE) -> out
	expect_class(out, 'gg')

	p_t.test() |>
		SpowerCurve(n=c(30, 90, 270, 550),
				   d=c(.2, .5, .8), replications=10, verbose=F, plotly = FALSE) -> out
	expect_class(out, 'gg')

	p_t.test(n=NA) |>
		SpowerCurve(d=c(.2, .5, .8), power=.80, interval=c(10, 1000), verbose=F, plotly = FALSE) -> out
	expect_class(out, 'gg')

	# p_t.test(n=250) |>
	# 	SpowerCurve(d=c(.2, .5, .6), power=.80, sig.level=NA, verbose=F) -> out
	# expect_class(out, 'gg')



})

test_that('scope', {

	mygen_fun <- function(n, n2_n1, d, df, ...){
		group1 <- rchisq(n, df=df)
		group1 <-  (group1 - df) / sqrt(2*df)   # Adjusted mean to 0, sd = 1
		group2 <- rnorm(n*n2_n1, mean=d)
		dat <- data.frame(group = factor(rep(c('G1', 'G2'),
											 times = c(n, n*n2_n1))),
						  DV = c(group1, group2))
		dat
	}

	p_my_t.test <- function(n, d, var.equal=FALSE, n2_n1=1, df=10, ...){
		dat <- mygen_fun(n=n, n2_n1=n2_n1, d=d, df=df, ...)
		obj <- t.test(DV ~ group, dat, var.equal=var.equal)

		# p-value must be first element when using default summarise()
		with(obj, c(p=p.value,
					mean_diff=unname(estimate[2] - estimate[1]),
					SE=stderr))
	}

	# Solve N to get .80 power (a priori power analysis), using defaults
	gg <- p_my_t.test(d = .5) |>
		SpowerCurve(n=c(30, 60, 90), verbose=F, replications=10, plotly=F)
	expect_is(gg, 'gg')

	gg <- p_my_t.test(d = .5) |>
		SpowerCurve(n=c(30, 60, 90), verbose=F, replications=50, plotly=F,
				   parallel=TRUE, ncores=2)
	expect_is(gg, 'gg')


})
