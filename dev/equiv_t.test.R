# Parameters
n_sim <- 10000    # Number of simulations
n_per_group <- 30 # Sample size per group
mu1 <- 0          # Mean of group 1
mu2 <- 0.5        # Mean of group 2
sigma <- 1        # Standard deviation (assumed equal)
alpha <- 0.05     # Significance level
margin <- 0.1     # Non-inferiority margin

# Function to perform a single simulation
simulate_t_test <- function() {
	group1 <- rnorm(n_per_group, mean = mu1, sd = sigma)
	group2 <- rnorm(n_per_group, mean = mu2, sd = sigma)
	t_test <- t.test(group2, group1, alternative = "greater", var.equal = TRUE)
	return(t_test$conf.int[1] > -margin)  # Check if lower bound is above -margin
}

# Run simulations
results <- replicate(n_sim, simulate_t_test())

# Calculate power
power <- mean(results)

cat("Estimated power:", power, "\n")


##################

library(Spower)

nequiv_t.test <- function(n, d, equivL, equivU=Inf, ...) {
	group1 <- rnorm(n)
	group2 <- rnorm(n, mean=d)
	t_test <- t.test(group2, group1, var.equal = TRUE, conf.level=.9) # .9 for non-inferiority, .95 for equivalence?
	CI <- t_test$conf.int
	nequiv <- CI[1] < equivL || CI[2] >  equivU
	nequiv
}

# non-inferiority (difference greater than .4)
nequiv_t.test(1000, .5, equivL=.4)

# Non-inferiority power
nequiv_t.test(1000, .5, equivL=.4) |> Spower(replications=1000)
nequiv_t.test(1000, .5, equivL=.35) |> Spower(replications=1000)

# Equivalence test (between .4 and .6)
nequiv_t.test(1000, .5, equivL=.4, equivU=.6)

# Equivalence power
nequiv_t.test(1000, .5, equivL=.4, equivU=.6) |> Spower(replications=1000)
nequiv_t.test(1000, .5, equivL=.3, equivU=.7) |> Spower(replications=1000)

# search for non-inferiority bound that would give 80% power (breaks? Issue with SimSolve?)
p_t.test(n=1000, d=NA) |> Spower(power=.95, interval=c(.01, .5))
nequiv_t.test(n=1000, d=.5, equivL=NA) |> Spower(power=.80, interval=c(.1, .5))

