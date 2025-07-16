context('Spower')

expect_class <- function(x, class) expect_true(inherits(x, class))

test_that('Spower', {

    library(Spower)

	set.seed(1234)
	out1 <- Spower(p_t.test(n = 50, d = .5), replications=10, verbose=FALSE)
	expect_class(out1, 'Spower')
	expect_equal(out1$power, .6)
	outu <- update(out1, sig.level = .10)
	expect_equal(outu$power, .8)

	# selection test
	ind.t.test <- function(n, d){
		g1 <- rnorm(n)
		g2 <- rnorm(n, mean=d)
		out <- t.test(g2, g1, var.equal=TRUE)
		c(p=out$p.value, xbar=out$estimate[1], se=out$stderr)
	}
	out1 <- ind.t.test(n = 50, d = .5) |>
		Spower(replications=10, select='p', verbose=FALSE)
	expect_class(out1, 'Spower')
	expect_equal(ncol(out1), 8)
	expect_equal(ncol(SimResults(out1)), 6)

	set.seed(4321)
	out2 <- Spower(p_t.test(n = NA, d = .5), power=.5, interval=c(10, 100),
				   maxiter = 40, verbose=FALSE)
	expect_class(out2, 'Spower')
	expect_equal(out2$n, 31.3, tolerance=.01)

	set.seed(42)
	out <- Spower(p_t.test(n = 50, d = .5), beta_alpha = 2,
				  verbose=FALSE, replications=1000)
	expect_class(out, 'Spower')
	expect_equal(out$power, .813, tolerance=.01)
	expect_equal(out$sig.level, .093, tolerance=.01)

	out2 <- update(out, beta_alpha=4)
	expect_equal(out2$power, .751, tolerance=.01)
	expect_equal(out2$sig.level, .062, tolerance=.01)

	set.seed(90210)
	out <- p_t.test(n = 50, d = .5) |>
		Spower(power=.90, sig.level=NA, verbose=FALSE, maxiter = 40)
	expect_class(out, 'Spower')
	expect_equal(out$power, .9, tolerance=.01)
	expect_equal(out$sig.level, .224, tolerance=.01)


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

test_that('multi', {

	p_my_t.test <- function(n, d){
		g1 <- rnorm(n)
		g2 <- rnorm(n, mean=d)
		p1 <- t.test(g1, g2, var.equal=FALSE)$p.value
		p2 <- t.test(g1, g2, var.equal=TRUE)$p.value
		c(welch=p1, ind=p2)
	}

	set.seed(90210)
	p_my_t.test(n=100, d=.2) |>
		Spower(replications=100, verbose=FALSE) -> sim
	expect_class(sim, 'Spower')
	expect_equal(sim$power.welch, .31)
	new_sim <- update(sim, sig.level=.01)
	expect_equal(new_sim$power.welch, .1)
	new_sim2 <- update(sim, beta_alpha=4)
	expect_equal(new_sim2$power, .494, tolerance=1e-2)

})

test_that('rerun', {

	set.seed(1234321)
	out <- p_t.test(n = 50, d = .5) |> Spower(replications=100, verbose=FALSE,
											  beta_alpha = 4)
	expect_equal(out$power, .740, tol=1e-2)
	p_t.test(n = 50, d = .5) |>
		Spower(replications=100, lastSpower=out,
			   verbose=FALSE, beta_alpha = 4) -> out2
	expect_equal(out2$REPLICATIONS, 200)
	expect_equal(out2$power, .759, tol=1e-2)


	set.seed(90210)
	out <- p_t.test(n = 50, d = .5) |> Spower(replications=100, verbose=FALSE)
	p_t.test(n = 50, d = .5) |>
		Spower(replications=100, lastSpower=out, verbose=FALSE) -> out2
	expect_equal(out2$REPLICATIONS, 200)
	expect_equal(attr(out2, 'extra_info')$SEED_history, c(1471380155, 1411519934))

	out <- p_t.test(n = NA, d = .5) |> Spower(power=.8, interval=c(10, 100),
											  maxiter=40, verbose=FALSE)
	expect_equal(out$n, 62.04, tolerance=1e-4)
	expect_equal(unname(summary(out)$predCIs_root), c(60.20545, 64.10115), tolerance=1e-4)

	p_t.test(n = NA, d = .5) |>
		Spower(power=.8, interval=c(10, 100), lastSpower=out, maxiter=70, verbose=FALSE) -> out2
	expect_equal(out2$n, 63.06, tolerance=1e-2)
	expect_equal(unname(summary(out2)$predCIs_root), c(62.37223, 63.73347), tolerance=1e-4)

	# multi
	p_my_t.test <- function(n, d){
		g1 <- rnorm(n)
		g2 <- rnorm(n, mean=d)
		p1 <- t.test(g1, g2, var.equal=FALSE)$p.value
		p2 <- t.test(g1, g2, var.equal=TRUE)$p.value
		c(welch=p1, ind=p2)
	}

	set.seed(54321)
	out <- p_my_t.test(n = 30, d = .5) |>
		Spower(replications=100, verbose=FALSE)
	p_my_t.test(n = 30, d = .5) |>
		Spower(replications=900, lastSpower=out, verbose=FALSE) -> out2
	expect_equal(c(out2$power.welch, out2$power.ind), c(.494, .497))
	expect_equal(out2$REPLICATIONS, 1000)
	expect_equal(attr(out2, 'extra_info')$SEED_history, c(1500916448, 1842577306))


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
		p <- obj$p.value
		p
	}

	p_my_t.test.defaults <- function(n = 30, d = .5,
									 var.equal=FALSE, n2_n1=1, df=10,
							gen_fun=mygen_fun, ...){
		dat <- gen_fun(n=n, n2_n1=n2_n1, d=d, df=df, ...)
		obj <- t.test(DV ~ group, dat, var.equal=var.equal)

		# p-value must be first element when using default summarise()
		p <- obj$p.value
		p
	}

	# Solve N to get .80 power (a priori power analysis), using defaults
	set.seed(1234)
	out <- Spower(p_my_t.test(n = NA, d = .5), power=.8,
				  interval=c(2,500), maxiter=40, verbose=FALSE)
	expect_equal(out$n, 66.2, tolerance=.01)

	out2 <- Spower(p_my_t.test(n = 100, d = .5),
				   replications=1000, verbose=FALSE)
	expect_equal(out2$power, 0.928, tolerance=.01)

	# in parallel
	set.seed(1234)
	out <- Spower(p_my_t.test(n = NA, d = .5), power=.8,
				  interval=c(2,500), maxiter=40, verbose=FALSE,
				  parallel=TRUE, ncores = 2)
	expect_equal(out$n, 67.4, tolerance=.01)

	out2 <- Spower(p_my_t.test(n = 100, d = .5),
				   verbose=FALSE, replications=1000,
				   parallel=TRUE, ncores = 2)
	expect_equal(out2$power, 0.934, tolerance=.01)

	# with assigned defaults
	set.seed(1234)
	out <- Spower(p_my_t.test.defaults(n=NA), power=.8,
				  interval=c(2,500), maxiter=40, verbose=FALSE)
	expect_equal(out$n, 66.2, tolerance=.01)

	out2 <- Spower(p_my_t.test.defaults(),
				   replications=1000, verbose=FALSE)
	expect_equal(out2$power, 0.491, tolerance=.01)

	# in parallel
	set.seed(1234)
	out <- Spower(p_my_t.test(n = NA, d = .5), power=.8,
				  interval=c(2,500), maxiter=40, verbose=FALSE,
				  parallel=TRUE, ncores = 2)
	expect_equal(out$n, 67.4, tolerance=.01)

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

