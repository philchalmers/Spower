context('Spower')

test_that('Spower', {

	expect_class <- function(x, class) expect_true(inherits(x, class))

    library(Spower)

	set.seed(1234)
	out1 <- Spower(p_t.test(n = 50, d = .5), replications=10, verbose=FALSE)
	expect_class(out1, 'Spower')
	expect_equal(out1$power, .6)
	outu <- update(out1, sig.level = .10)
	expect_equal(outu$power, .8)

	set.seed(4321)
	out2 <- Spower(p_t.test(n = NA, d = .5), power=.5, interval=c(10, 100),
				   maxiter = 40, verbose=FALSE)
	expect_class(out2, 'Spower')
	expect_equal(out2$n, 32.15, tolerance=.01)

	set.seed(42)
	out <- Spower(p_t.test(n = 50, d = .5), beta_alpha = 2,
				  verbose=FALSE, replications=1000)
	expect_class(out, 'Spower')
	expect_equal(out$power, .813, tolerance=.01)
	expect_equal(out$sig.level, .093, tolerance=.01)

	out2 <- update(out, beta_alpha=4)
	expect_equal(out2$power, .751, tolerance=.01)
	expect_equal(out2$sig.level, .062, tolerance=.01)

	# same, but in parallel
	set.seed(1234)
	out1 <- Spower(p_t.test(n = 50, d = .5), replications=10, verbose=FALSE,
				   parallel=TRUE, ncores=2)
	expect_class(out1, 'Spower')
	expect_equal(out1$power, .7)
	outu <- update(out1, sig.level = .10)
	expect_equal(outu$power, .8)

	set.seed(4321)
	out2 <- Spower(p_t.test(n = NA, d = .5), power=.5, interval=c(10, 100),
				   maxiter = 40, verbose=FALSE,  parallel=TRUE, ncores=2)
	expect_equal(out2$n, 31.00, tolerance=.01)

	set.seed(42)
	out <- Spower(p_t.test(n = 50, d = .5), beta_alpha = 2,
				  verbose=FALSE, replications=1000,
				  parallel=TRUE, ncores=2)
	expect_equal(out$power, 0.790, tolerance=.01)
	expect_equal(out$sig.level, 0.104, tolerance=.01)

	out2 <- update(out, beta_alpha=4)
	expect_equal(out2$power, 0.726, tolerance=.01)
	expect_equal(out2$sig.level, .068, tolerance=.01)


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

	p_my_t.test <- function(n, d, var.equal=FALSE, n2_n1=1, df=10,
							gen_fun=mygen_fun, ...){
		dat <- gen_fun(n=n, n2_n1=n2_n1, d=d, df=df, ...)
		obj <- t.test(DV ~ group, dat, var.equal=var.equal)

		# p-value must be first element when using default summarise()
		with(obj, c(p=p.value,
					mean_diff=unname(estimate[2] - estimate[1]),
					SE=stderr))
	}

	# Solve N to get .80 power (a priori power analysis), using defaults
	set.seed(1234)
	out <- Spower(p_my_t.test(n = NA, d = .5), power=.8,
				  interval=c(2,500), maxiter=40, verbose=FALSE)
	expect_equal(out$n, 64.42, tolerance=.01)

	out2 <- Spower(p_my_t.test(n = 100, d = .5),
				   replications=1000, verbose=FALSE)
	expect_equal(out2$power, 0.928, tolerance=.01)

	# in parallel
	set.seed(1234)
	out <- Spower(p_my_t.test(n = NA, d = .5), power=.8,
				  interval=c(2,500), maxiter=40, verbose=FALSE,
				  parallel=TRUE, ncores = 2)
	expect_equal(out$n, 66.01, tolerance=.01)

	out2 <- Spower(p_my_t.test(n = 100, d = .5),
				   verbose=FALSE, replications=1000,
				   parallel=TRUE, ncores = 2)
	expect_equal(out2$power, 0.934, tolerance=.01)

	# constants
	set.seed(4321)
	n <- 206
	n2_n1 <- 51/n
	p_2r(n=n, r.ab1=.75, r.ab2=.88, n2_n1=n2_n1) |>
		Spower(replications=100, verbose=FALSE) -> out
	expect_equal(out$power, .72)
	p_2r(n=n, r.ab1=.75, r.ab2=.88, n2_n1=n2_n1) |>
		Spower(replications=100, verbose=FALSE,
			   parallel=TRUE, ncores=2) -> out2
	expect_equal(out2$power, .68)

	# conflicting constants
	set.seed(123)
	PI <- pi <- .65
	g <- .15
	p <- PI + g
	p_prop.test(n=20, prop=p, pi=PI, two.tailed=FALSE) |>
		Spower(replications=10, verbose=FALSE) -> out
	expect_equal(out$power, .6)

	# global pi
	p_prop.test(n=20, prop=p, pi=pi, two.tailed=FALSE) |>
		Spower(replications=10, verbose=FALSE) -> out2
	expect_equal(out2$power, .2)

	# unnamed inputs
	p_t.test(100, .5) |> Spower(replications=10, verbose=FALSE) -> out
	expect_is(out, 'Spower')
})

