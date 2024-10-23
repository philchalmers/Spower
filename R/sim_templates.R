#' @export
p_t.test <- function(N, d) {
	group1 <- rnorm(N)
	group2 <- rnorm(N, mean=d)
	dat <- data.frame(group = gl(2, N, labels=c('G1', 'G2')),
					  DV = c(group1, group2))
	p <- t.test(DV ~ group, dat, var.equal=TRUE)$p.value
	p
}
