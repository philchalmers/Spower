# power for multi-level model
# https://stats.stackexchange.com/questions/446803/power-analysis-for-glmer-using-simr

simulate_binary <- function (n) {
	K <- 8 # number of measurements per subject
	t_max <- 15 # maximum follow-up time

	# we constuct a data frame with the design:
	# everyone has a baseline measurement, and then measurements at random follow-up times
	DF <- data.frame(id = rep(seq_len(n), each = K),
					 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
					 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

	# design matrices for the fixed and random effects
	X <- model.matrix(~ sex * time, data = DF)
	Z <- model.matrix(~ time, data = DF)

	betas <- c(-2.13, -0.25, 0.24, -0.05) # fixed effects coefficients
	D11 <- 0.48 # variance of random intercepts
	D22 <- 0.1 # variance of random slopes

	# we simulate random effects
	b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
	# linear predictor
	eta_y <- drop(X %*% betas + rowSums(Z * b[DF$id, ]))
	# we simulate binary longitudinal data
	DF$y <- rbinom(n * K, 1, plogis(eta_y))
	DF
}

###################################################################

if(FALSE){
	library("GLMMadaptive")
	M <- 1000 # number of simulations to estimate power
	p_values <- numeric(M)
	for (m in seq_len(M)) {
		DF_m <- simulate_binary(n = 100)
		fm_m <- mixed_model(y ~ sex * time, random = ~ time | id,
							data = DF_m, family = binomial())
		p_values[m] <- coef(summary(fm_m))["sexfemale:time", "p-value"]
	}
	# assuming a significance level of 5%, the power will be
	mean(p_values < 0.05)
}


###############
# for Spower
p_mixed_model <- function(n){
	DF_m <- simulate_binary(n = n)
	fm_m <- mixed_model(y ~ sex * time, random = ~ time | id,
						data = DF_m, family = binomial())
	p <- coef(summary(fm_m))["sexfemale:time", "p-value"]
	p
}

# estimate power given n
Spower(p_mixed_model, n=100, replications=1000,
	   packages="GLMMadaptive")

# estimate n to achieve 80% power
Spower(p_mixed_model, n=NA, power=.8,
	   interval=c(100, 1000), packages="GLMMadaptive")

# estimate power given n
Spower(p_mixed_model, n=100, replications=1000, parallel=TRUE,
	   packages="GLMMadaptive")

# estimate n to achieve 80% power
Spower(p_mixed_model, n=NA, power=.8,
	   interval=c(100, 1000), parallel=TRUE, packages="GLMMadaptive")
