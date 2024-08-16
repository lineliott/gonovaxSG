
sample_uptake <- function(offered, accepted, n) {
  rbeta(n, accepted, offered - accepted)
}

sample_gamma <- function(mean, n, se) {
  var <- (se * mean) ^ 2
  scale <- var / mean
  shape <- mean / scale
  rgamma(n = n, shape = shape, scale = scale)
}

fit_beta <- function(mean, l, u, ci = 0.95) {
  a <- (1 - ci) / 2
  p <- c(a, 1 - a)
  
  f <- function(alpha) {
    beta <- alpha * (1 - mean) / mean
    x <- qbeta(p = p, shape1 = alpha, shape2 = beta)
    sum((x - c(l, u)) ^ 2)
  }
  
  alpha <- optimise(f = f, interval = c(0, 200), maximum = FALSE)$minimum
  beta <- alpha * (1 - mean) / mean
  qs <- qbeta(p, shape1 = alpha, shape2 = beta)
  
  print(sprintf("fitted qs = (%.3f, %.3f)", qs[1], qs[2]))
  print(sprintf("target qs = (%f, %f)", l, u))
  
  list(alpha = alpha, beta = beta[[1]])
  
}

plot_param <- function(params, ...) {
  for (i in seq_along(params)) {
    x <- params[[i]]
    hist(x, main = "", freq = FALSE,
         xlab = names(params)[i], ...)
    abline(v = mean_ci(x), lty = c(1, 2, 2), col = "darkred")
  }
}

