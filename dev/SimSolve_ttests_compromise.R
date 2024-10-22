#' ---
#' title: "Compromise analysis example"
#' output:
#'   html_document:
#'     theme: readable
#'     code_download: true
#' ---

library(SimDesign)

Design <- createDesign(crit = NA,
                       N = 100,
					   d = .3)
Design    # solve for NA's

#~~~~~~~~~~~~~~~~~~~~~~~~
#### Step 2 --- Define generate, analyse, and summarise functions

Generate <- function(condition, fixed_objects = NULL) {
	Attach(condition)
	group1 <- rnorm(N)
	group2 <- rnorm(N, mean=d)
	dat <- data.frame(DV.a = group1,
	                  DV.b = group2)
	dat
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
	t.a <- t.test(dat$DV.a)$statistic
	t.b <- t.test(dat$DV.b)$statistic
	with(condition, c(a=t.a > crit, b=t.b < crit))
}

Summarise <- function(condition, results, fixed_objects = NULL) {
    means <- colMeans(results)
    q <- means[2] / means[1]
    q
}

#~~~~~~~~~~~~~~~~~~~~~~~~
#### Step 3 --- Optimize N over the rows in design

# In this example, b = desired ratio
solved <- SimSolve(design=Design, b=1, interval=c(1, 2),
				   generate=Generate, analyse=Analyse,
				   summarise=Summarise, integer = FALSE,
				   family='gaussian',
				   control=list(summarise.reg_data=TRUE))
solved
summary(solved)
plot(solved, 1)


# also can plot median history and estimate precision
plot(solved, 1, type = 'history')

# check using better functions
fun <- \(alpha, N=100, d=.3, target_ratio = 1) {
    beta <- alpha * target_ratio
    abs(qt(1 - alpha, ncp = 0, df = N - 1)  # crit for alpha
        -
            qt(beta, ncp = d * sqrt(N), df = N - 1)  # crit for beta
    )
}
alpha <- optim(0.1, fun, lower = 1e-09, upper = 1 -
                   1e-09, method = "L-BFGS-B")$par
N <- 100
crit <- qt(1 - alpha, ncp = 0, df = N - 1)
crit



#####################
# solve for alpha cutoff rather than critical value 
library(SimDesign)

Design <- createDesign(N = 100,
                       d = .3, 
                       alpha = .05)
Design    # solve for NA's

Attach(Design, RStudio_flags = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~
#### Step 2 --- Define generate, analyse, and summarise functions

Generate <- function(condition, fixed_objects = NULL) {
    Attach(condition)
    group1 <- rnorm(N)
    group2 <- rnorm(N, mean=d)
    dat <- data.frame(DV = c(group1, group2),
                      group = rep(c('G1', 'G2'), each=N))
    dat
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
    p <- t.test(DV ~ group, data=dat)$p.value
    p
}

Summarise <- function(condition, results, fixed_objects = NULL) {
    rate <- EDR(results, alpha=condition$alpha)
    ret <- c(beta_alpha = unname((1-rate) / condition$alpha))
    ret
}

#~~~~~~~~~~~~~~~~~~~~~~~~
#### Step 3 --- Optimize N over the rows in design

sim <- runSimulation(design=Design, replications=10000, 
                     generate=Generate, analyse=Analyse,
                     summarise=Summarise, store_results=TRUE)
sim

# compute beta/alpha ratio given different alpha
compromise <- function(alpha, sim, Design){
    Design$alpha <- alpha
    out <- reSummarise(Summarise, results=sim, Design=Design)    
    out$beta_alpha
}

compromise(.3, sim=sim, Design=Design)
compromise(.001, sim=sim, Design=Design)

compromise_root <- function(alpha, beta_alpha, ...)
    compromise(alpha, ...) - beta_alpha


# equal beta/alpha trade-off  
uniroot(compromise_root, c(.001, .3), beta_alpha=1, 
        sim=sim, Design=Design)

# beta/alpha trade-off 4 times worse (beta = 4 * alpha)
uniroot(compromise_root, c(.001, .3), beta_alpha=4, 
        sim=sim, Design=Design)


#####################
# same as above, however if Type I error not nominal then may wish to use 
# empirical Type I error estimate instead
library(SimDesign)

Design <- createDesign(N = 100,
                       d = .3, 
                       alpha = .05)
Design    # solve for NA's

#~~~~~~~~~~~~~~~~~~~~~~~~
#### Step 2 --- Define generate, analyse, and summarise functions

Generate <- function(condition, fixed_objects = NULL) {
    Attach(condition)
    group1 <- rnorm(N)
    group2 <- rnorm(N, mean=d)
    group3 <- rnorm(N)    # for H0 test
    dat <- data.frame(DV = c(group1, group2),
                      DV.null = c(group1, group3),
                      group = rep(c('G1', 'G2'), each=N))
    dat
}

Analyse <- function(condition, dat, fixed_objects = NULL) {
    p.beta <- t.test(DV ~ group, data=dat)$p.value # (1-beta)
    p.alpha <- t.test(DV.null ~ group, data=dat)$p.value  
    nc(p.alpha, p.beta)
}

Summarise <- function(condition, results, fixed_objects = NULL) {
    rate <- EDR(results, alpha=condition$alpha)
    ret <- c(beta_alpha = unname((1-rate["p.beta"]) / rate["p.alpha"]))
    ret
}

#~~~~~~~~~~~~~~~~~~~~~~~~
#### Step 3 

sim <- runSimulation(design=Design, replications=10000, 
                     generate=Generate, analyse=Analyse,
                     summarise=Summarise, store_results=TRUE)
sim

# compute beta/alpha ratio given different alpha
compromise <- function(alpha, sim, Design){
    Design$alpha <- alpha
    out <- reSummarise(Summarise, results=sim, Design=Design)    
    out$beta_alpha
}

compromise(.3, sim=sim, Design=Design)
compromise(.001, sim=sim, Design=Design)

compromise_root <- function(alpha, beta_alpha, ...)
    compromise(alpha, ...) - beta_alpha


# equal beta/alpha trade-off  
uniroot(compromise_root, c(.001, .3), beta_alpha=1, 
        sim=sim, Design=Design)

# beta/alpha trade-off 4 times worse (beta = 4 * alpha)
uniroot(compromise_root, c(.001, .3), beta_alpha=4, 
        sim=sim, Design=Design)
