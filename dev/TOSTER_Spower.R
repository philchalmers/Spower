# Equivalence test example using TOSTER and Spower
library(Spower)
library(TOSTER)

equiv.t.test <- function(n, m1, sd1, m2, sd2, LB, UB, sig.level = .05){
	g1 <- rnorm(n, m1, sd1)
	g2 <- rnorm(n, m2, sd2)
	out <- tsum_TOST(m1=mean(g1), sd1=sd(g1),
			  m2=mean(g2), sd2=sd(g2),
			  n1=n, n2=n,
			  low_eqbound=LB,
			  high_eqbound=UB,
			  eqbound_type = "raw",
			  paired = FALSE) # This is quite slow due to extra computations
	out$TOST[2,]$p.value < sig.level/2 &
		out$TOST[3,]$p.value < sig.level/2
}

# These are equivalent
equiv.t.test(n=100, m1=12, m2=11, sd1=2.5, sd2=2.5, LB=-2.5, UB=2.5) |>
	Spower(replications=1000)

power_t_TOST(n = 100,
			 delta = 1,
			 sd = 2.5,
			 eqb = 2.5,
			 alpha = .025,
			 power = NULL,
			 type = "two.sample")

# solve sample size
equiv.t.test(n=NA, m1=12, m2=11, sd1=2.5, sd2=2.5, LB=-2.5, UB=2.5) |>
	Spower(power=.95, interval=c(50, 200))

power_t_TOST(n = NULL,
			 delta = 1,
			 sd = 2.5,
			 eqb = 2.5,
			 alpha = .025,
			 power = .95,
			 type = "two.sample")


######################

# Using CI approach (much faster)
CI.equiv.t.test <- function(n, m1, sd1, m2, sd2, LB, UB, conf.level=.95){
	g1 <- rnorm(n, m1, sd1)
	g2 <- rnorm(n, m2, sd2)
	CI <- t.test(g2, g1, var.equal=TRUE,
				 conf.level=conf.level)$conf.int
	# LB < CI[1] && CI[2] < UB             # manually
	is.CI_within(CI, interval = c(LB, UB)) # this is equivalent
}

CI.equiv.t.test(n=100, m1=12, m2=11, sd1=2.5, sd2=2.5, LB=-2.5, UB=2.5) |>
	Spower()

# solve sample size
CI.equiv.t.test(n=NA, m1=12, m2=11, sd1=2.5, sd2=2.5, LB=-2.5, UB=2.5) |>
	Spower(power=.95, interval=c(50, 200))
