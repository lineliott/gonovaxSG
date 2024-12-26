# calibration support was adapted from gonovax pacakge uk for singapore setting
version_check <- function(package, version) {
  if (packageVersion(package) < version) {
    stop(sprintf(
      paste("Please update %s with:",
            "devtools::install_github('mrc-ide/%s', force = TRUE)"),
      package, package))
  }
  
}

# demograohic info for singapore [adapted from gonovax package]
demographic_params_sg <- list(    
  # msm
  N0 = 139000,
  enr= 2780,
  male_enr = 19782.5,

  q = c(0.85, 0.15), # proportion in group L
  p = c(0.6, 15.6), # partner change rate in group L/H

  exr = 1 / 50
)

mcmc_control <- function(short_run) {
  control <- list(
    n_mcmc   = 5e5,
    #n_mcmc = 2e5,
    n_chains = 8,
    n_sample = 1e3,
    burnin   = 1e5)
  if (short_run) {
    control$n_mcmc <- 1e4
    control$n_sample <- 1e3
    #control$burnin <- 50
    control$burnin <- 100
  }
  control
}

transform_to_normal_distribution <- function(posterior_estimate) {
  posterior_mean <- posterior_estimate[1]  # Extracting the posterior mean
  credible_interval <- posterior_estimate[-1]  # Extracting the credible interval
  
  # Calculating the standard deviation based on the credible interval
  half_width <- (credible_interval[2] - credible_interval[1]) / 2
  z_score <- qnorm(0.975)  # Fowr a 95% credible interval
  std_dev <- half_width / z_score
  
  # Returning the parameters for the normal distribution
  gaussian_mean <- posterior_mean
  gaussian_sd <- std_dev
  
  return(list(gaussian_mean = gaussian_mean, gaussian_sd = gaussian_sd))
}



prepare_priors <- function(prior) {
  lapply(split(prior, prior$name), make_prior)
}

prepare_parameters <- function(info, priors, proposal) {
  
  pars <- Map(mcstate::pmcmc_parameter,
              name     = info$name,
              initial  = info$initial,
              min      = info$min,
              max      = info$max,
              # eliott
              integer = info$discrete,
              #discrete = info$discrete,
              prior    = priors[info$name])
  
  proposal <- as.matrix(proposal)[info$name, info$name]
  
  pars <- mcstate::pmcmc_parameters$new(pars, proposal)
  
}


make_prior <- function(d) {
  if (d$type == "gamma") {
    
    shape <- d$gamma_shape
    scale <- d$gamma_scale
    
    if (d$as_duration) {
      function(p) {
        dgamma(1 / p, shape = shape, scale = scale, log = TRUE)
      }
    } else {
      function(p) {
        dgamma(p, shape = shape, scale = scale, log = TRUE)
      }
    }
  } 
  else if (d$type == "beta") {
    shape1 <- d$beta_shape1
    shape2 <- d$beta_shape2
    function(p) {
      dbeta(p, shape1 = shape1, shape2 = shape2, log = TRUE)
    }
  } 
  else if (d$type == "lnorm") {
    function(p) {
      dlnorm(p, meanlog = d$lnorm_meanlog, sdlog = d$lnorm_sdlog, log = TRUE)
    }
  }
  # eliott
  else if (d$type == "gaussain") {
    mean <- d$gaussian_mean
    sd <- d$gaussian_sd
    function(p) {
      dnorm(p, mean = mean, sd = sd, log = TRUE)
    }
  }
  else if (d$type == "null") {
    NULL
  } 
  else {
    stop(d, "Unknown prior type")
  }
}

calculate_parameters_info <- function(results, parameters_info) {
  i_mpe <- which.max(results$probabilities[, "log_posterior"])
  pars_mpe <- results$pars[i_mpe, ]
  parameters_info$initial <- pars_mpe
  parameters_info
}

calculate_mcmc_metrics <- function(mcmc_results) {
  
  results <- lapply(mcmc_results, function(x) coda::as.mcmc(x$pars))
  results <- coda::as.mcmc.list(results)
  ret <- list(gr  = coda::gelman.diag(results, multivariate = FALSE),
              ess = coda::effectiveSize(results),
              acceptance_rate = 1 - coda::rejectionRate(results))
  ret
}


create_vax_map <- function(n_vax, v, i_u, i_v) {
  
  # ensure vaccine input is of correct length
  n_group <- 2
  
  stopifnot(length(v) == n_group)
  stopifnot(all(v %in% c(0, 1)))
  if (length(i_v) > 0) {
    stopifnot(max(i_u, i_v) <= n_vax)
  }
  
  # set up vaccination matrix
  vax_map <- array(0, dim = c(n_group, n_vax, n_vax))
  
  for (i in seq_along(i_u)) {
    vax_map[, i_u[i], i_u[i]] <-  v
    vax_map[, i_v[i], i_u[i]] <- -v
  }
  
  vax_map
}

vax_params0 <- function(n_diag_rec = 1) {
  n_group <- 2
  n_vax <- n_diag_rec
  v <- array(0, dim = c(n_group, n_vax, n_vax))
  list(n_vax = n_diag_rec,
       willing = 1,
       u_vbe = 0,
       u = v,
       vbe = v,
       vos = v,
       vod = v,
       vea = 0,
       ved = 0,
       ves = 0,
       vei = 0,
       w = as.matrix(0),
       D = 0,
       vax_t = c(0, 99),
       vax_y = c(0, 0))
}

initial_params <- function(pars, n_vax = 1, coverage = 1) {
  
  stopifnot(length(coverage) == n_vax)
  stopifnot(sum(coverage) == 1)
  
  U0 <- I0 <- A0 <- S0 <- T0 <- array(0, c(2, n_vax))
  
  # separate into 1:low and 2:high activity groups and by coverage
  N0 <- pars$N0 * outer(pars$q, coverage)
  # set initial asymptomatic prevalence in each group (unvaccinated only)
  A0[, 1] <- round(N0[, 1] * c(pars$prev_Asl, pars$prev_Ash))
  
  # set initial uninfecteds
  U0 <- round(N0) - A0
  
  list(U0 = U0, I0 = I0, A0 = A0, S0 = S0, T0 = T0)
}


assert_is <- function(x, what, name = deparse(substitute(x))) {
        if (!inherits(x, what)) {
          stop(sprintf("'%s' must be a %s", name,
                       paste(what, collapse = " / ")), call. = FALSE)
        }
        invisible(x)
      }
      
assert_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1L) {
    stop(sprintf("'%s' must be a scalar", name), call. = FALSE)
  }
  invisible(x)
}

assert_integer <- function(x, name = deparse(substitute(x)),
                           what = "integer") {
  if (!(is.integer(x))) {
    eps <- sqrt(.Machine$double.eps)
    usable_as_integer <- is.numeric(x) && (max(abs(round(x) - x)) < eps)
    if (!usable_as_integer) {
      stop(sprintf("'%s' must be an %s", name, what), call. = FALSE)
    }
    x <- as.integer(round(x))
  }
  invisible(x)
}
assert_scalar_positive <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  if (x <= 0) {
    stop(sprintf("'%s' must be greater than 0", name), call. = FALSE)
  }
  invisible(x)
}

assert_scalar_positive_integer <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  x <- assert_integer(x, name)
  if (x < 1L) {
    stop(sprintf("'%s' must be at least 1", name), call. = FALSE)
  }
  invisible(x)
}

assert_scalar_unit_interval <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  if (x < 0 || x > 1) {
    stop(sprintf("'%s' must be between 0 and 1", name), call. = FALSE)
  }
  invisible(x)
}

check_gono_params <- function(pars) {
  with(pars, {
    assert_scalar_unit_interval(psi)
    assert_scalar_unit_interval(prev_Asl)
    assert_scalar_unit_interval(prev_Ash)
    assert_scalar_unit_interval(epsilon)
    assert_scalar_positive(sigma)
    assert_scalar_positive(nu)
    assert_scalar_positive(mu)
    assert_scalar_positive(rho)
  })
}

gonovax_year_eliott <- function(year) {
  years_after_2004 <- as.numeric(year - 2004)
  if (any(years_after_2004 < 0)) {
    stop("Negative dates, gonovax_year likely applied twice")
  }
  years_after_2004
}

transform_eliott <- function(pars, fix_par_t = TRUE) {
  # reformat time-varying pars
  pars <- as.list(pars)
  check_gono_params(pars)
  with(pars, {
    assert_scalar_positive(beta_2004)
    assert_scalar_positive(beta_2018)
    assert_scalar_positive(eta_h)
    assert_scalar_unit_interval(omega)
  })


  t0 <- 2004
  t1 <- 2018
  t_max <- 1e2


  t_fix <- ifelse(fix_par_t, t1, t0 + t_max)
  pars$tt <- gonovax_year_eliott(pmin(seq(t0, t0 + t_max), t_fix))

  pars$beta_t <-  ( ((pars$beta_2018 - pars$beta_2004)/(t1-t0)) * pars$tt) + pars$beta_2004
  pars$eta_l_t <- pars$eta_h * pars$omega * (1 + 0 * pars$tt)
  pars$eta_h_t <- pars$eta_h * (1 + 0 * pars$tt)

  pars$tt <- seq(0, t_max)
  pars
}




##' Run a mcmc sampler [adapted from gnovax package]
##'
##' This is a basic Metropolis-Hastings MCMC sampler.  The
##' `model` is run with a set of parameters to evaluate the
##' likelihood. A new set of parameters is proposed, and these
##' likelihoods are compared, jumping with probability equal to their
##' ratio. This is repeated for `n_steps` proposals.
##'
##' This function is adapted from the `pmcmc` in the mcstate package
##' https://github.com/mrc-ide/mcstate
##'
##' @title Run a mcmc sampler
##'
##' @param pars A [`mcmc_parameters`] object containing
##'   information about parameters ( parameter ranges, priors, proposal kernel,
##'   observation functions).
##'
##' @param n_steps Number of MCMC steps to run
##'
##' @param compare likelihood function to compare data to epidemic trajectory
##' should return a single value representing the log-likelihood. Default is
##' compare_basic()
##'
##' @param progress Logical, indicating if a progress bar should be
##'   displayed, using [`progress::progress_bar`].
##'
##' @param n_chains Optional integer, indicating the number of chains
##'   to run. If more than one then we run a series of chains and
##'   merge them with [mcmc_combine()]. Chains are run in series,
##'   with the same model.
##'
##' @return A `gonovax_mcmc` object containing `pars`
##'   (sampled parameters) and `probabilities` (log prior, log
##'   likelihood and log posterior values for these
##'   probabilities).
##'
##' @export
# ##' @importFrom stats runif dnorm
mcmc_eliott <- function(pars, n_steps, compare = NULL, progress = FALSE,
                        n_chains = 1) {
  
  assert_is(pars, "pmcmc_parameters")
  #assert_scalar_positive_integer(n_steps)
  assert_scalar_positive_integer(control$n_chains)
  assert_scalar_positive_integer(control$n_mcmc)
  
  compare <- compare %||% compare_basic
  
  if (n_chains == 1) {
    mcmc_single_chain(pars, n_steps, compare, progress)
  } else {
    samples <- vector("list", n_chains)
    for (i in seq_along(samples)) {
      if (progress) {
        message(sprintf("Running chain %d / %d", i, n_chains))
      }
      samples[[i]] <- mcmc_single_chain(pars, n_steps, compare, progress)
    }
    gonovax::mcmc_combine(samples = samples)
  }
}



mcmc_single_chain <- function(pars, n_steps, compare, progress) {

  history_pars <- history_collector(n_steps)
  history_probabilities <- history_collector(n_steps)

  curr_pars <- pars$initial()
  curr_lprior <- pars$prior(curr_pars)
  curr_llik <- compare(curr_pars)
  curr_lpost <- curr_lprior + curr_llik

  history_pars$add(curr_pars)
  history_probabilities$add(c(curr_lprior, curr_llik, curr_lpost))

  #tick <- mcmc_progress(n = n_steps[1], show = progress)
  tick <- mcmc_progress(n_steps, progress)

  for (i in seq_len(n_steps)) {
    tick()

    # prop_pars <- pars$propose(curr_pars)
    # prop_lprior <- pars$prior(prop_pars)
    # prop_llik <- compare(prop_pars)
    prop_llik <- NaN  # Initialize prop_llik to NaN

    while (is.nan(prop_llik)) {
      prop_pars <- pars$propose(curr_pars)
      prop_lprior <- pars$prior(prop_pars)
      prop_llik <- compare(prop_pars)
    }
    # if (any(is.nan(prop_llik))) {
    #   print(paste0("nan prop_llik produced"))
    # }

    prop_lpost <- prop_lprior + prop_llik

    if (runif(1) < exp(prop_lpost - curr_lpost)) {
      curr_pars <- prop_pars
      curr_lprior <- prop_lprior
      curr_llik <- prop_llik
      curr_lpost <- prop_lpost
    }

    history_pars$add(curr_pars)
    history_probabilities$add(c(curr_lprior, curr_llik, curr_lpost))
  }
  pars_matrix <- do.call(rbind, history_pars$get())
  probabilities <- do.call(rbind, history_probabilities$get())
  colnames(probabilities) <- c("log_prior", "log_likelihood", "log_posterior")

  gonovax_mcmc(pars_matrix, probabilities)
}



## Generic history collector, collects anything at all into a list
##
## This would be more nicely done as a simple R6 class but it's a bit
## slow in testing; this version speeds up the total mcmc runtime by a
## factor of ~3x (0.4s/1000 iterations to 0.13s/1000) mostly by
## reducing the number of garbage collections considerably.
history_collector <- function(n) {
  data <- vector("list", length(n) + 1L)
  i <- 0L
  add <- function(value) {
    i <<- i + 1L
    data[[i]] <<- value
  }
  
  get <- function() {
    data
  }
  
  list(add = add, get = get)
}


##' @title Calculate the log likelihood of the data given the parameters
##' diagnoses and attendances lhoods are negative binomial
##' p_symp lhood is betabinomial
##' @param pars A named vector of parameters
##' @param transform the transform function to use in the comparison
##'
##' @return a single log likelihood
##' @importFrom stats dnbinom
##' @export

compare_eliott <- function(pars, a) {
  
  
  # run odin model
  y <- gonovax::run(tt = c(0, sg_data$gonovax_year), 
                    gono_params = transform_eliott(pars), 
                    demographic_params = demographic_params_sg,
                    transform = FALSE)
  
  # output total treated and total attended per year from odin model
  diagnosed <- diff(y[, "tot_treated"])
  
  # ## check if diagnosed [1:16] values are all positive
  # if (any(diagnosed[1:15] <= 0)) {
  #   message("Some diagnosed values are non-positive. Further inspection required.")
  #   print(paste0("NaNs diagnosed ", toString(diagnosed)))
  #   print(paste0("NaNs beta_2004 ", toString(transform_eliott(pars)[['beta_2004']])))
  #   print(paste0("NaNs beta_2018 ", toString(transform_eliott(pars)[['beta_2018']])))
  #   print(paste0("NaNs beta_t ", toString(transform_eliott(pars)[['beta_t']][1:16])))
  #   print(paste0("NaNs eta_l_t ", toString(transform_eliott(pars)[['eta_l_t']][1])))
  #   print(paste0("NaNs eta_h_t ", toString(transform_eliott(pars)[['eta_h_t']][1])))
  #   # pars$beta_t <-  (((pars$beta_2018 - pars$beta_2004)/(t1-t0)) * pars$tt) + pars$beta_2004
  #   # pars$eta_l_t <- pars$eta_h * pars$omega * (1 + 0 * pars$tt)
  #   # pars$eta_h_t <- pars$eta_h * (1 + 0 * pars$tt)
  #   #
  #   # pars$tt <- seq(0, t_max)
  # 
  # }
  
  ## compare to data
  ## in dnbinom size = shape param of gamma dist, mu = mean
  
  size <- 1 / pars[["k_surveillance"]]
  lltreated <- dnbinom(sg_data$diagnosed, size = size, mu = diagnosed, log = TRUE)
  
  # ## check if lltreated is NaN
  # if (any(is.nan(lltreated))) {
  #   print(paste0("NaNs beta_2004 ", toString(pars[['beta_2004']])))
  #   print(paste0("NaNs beta_2018 ", toString(pars[['beta_2018']])))
  #   print(paste0("NaNs diagnosed ", toString(diagnosed)))
  #   print(paste0("NaNs size ", toString(size)))
  # }
  # 
  ## output log likelihood
  if (any(is.nan(lltreated))) {
    NaN
  } else {
    sum(lltreated, na.rm = TRUE)
  }
  
}


# this is imported from mcstate

mcmc_progress <- function(n, show , force = FALSE) {
  if (show) {
    fmt <- "Step :current / :total [:bar] ETA :eta | :elapsedfull so far"
    t0 <- Sys.time()
    callback <- function(p) {
      message(sprintf("Finished %d steps in %s",
                      n, format(Sys.time() - t0, digits = 1)))
    }
    p <- progress::progress_bar$new(fmt, n, callback = callback, force = force)
    p$tick(0)
    p$tick
  } else {
    function() NULL
  }
}

gonovax_mcmc <- function(pars, probabilities, chain = NULL,
                         iteration = NULL) {
  
  ret <- list(chain = chain,
              iteration = iteration %||% seq.int(0, length.out = nrow(pars)),
              pars = pars,
              probabilities = probabilities)
  class(ret) <- "gonovax_mcmc"
  ret
}

assert_is <- function(x, what, name = deparse(substitute(x))) {
  if (!inherits(x, what)) {
    stop(sprintf("'%s' must be a %s", name,
                 paste(what, collapse = " / ")), call. = FALSE)
  }
  invisible(x)
}


assert_named <- function(x, unique = FALSE, name = deparse(substitute(x))) {
  if (is.null(names(x))) {
    stop(sprintf("'%s' must be named", name), call. = FALSE)
  }
  if (unique && any(duplicated(names(x)))) {
    stop(sprintf("'%s' must have unique names", name), call. = FALSE)
  }
}


assert_integer <- function(x, name = deparse(substitute(x)),
                           what = "integer") {
  if (!(is.integer(x))) {
    eps <- sqrt(.Machine$double.eps)
    usable_as_integer <- is.numeric(x) && (max(abs(round(x) - x)) < eps)
    if (!usable_as_integer) {
      stop(sprintf("'%s' must be an %s", name, what), call. = FALSE)
    }
    x <- as.integer(round(x))
  }
  invisible(x)
}

assert_positive_integer <- function(x, name = deparse(substitute(x))) {
  force(name)
  x <- assert_integer(x, name)
  if (any(x < 1L)) {
    stop(sprintf("'%s' must be at least 1", name), call. = FALSE)
  }
  invisible(x)
}


assert_character <- function(x, name = deparse(substitute(x))) {
  if (!(is.character(x))) {
    stop(sprintf("'%s' must be a character", name), call. = FALSE)
  }
  invisible(x)
}


assert_strictly_increasing <- function(x, name = deparse(substitute(x))) {
  if (any(diff(x) <= 0)) {
    stop(sprintf("'%s' must be strictly increasing", name), call. = FALSE)
  }
  invisible(x)
}


assert_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1L) {
    stop(sprintf("'%s' must be a scalar", name), call. = FALSE)
  }
  invisible(x)
}


assert_scalar_positive_integer <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  x <- assert_integer(x, name)
  if (x < 1L) {
    stop(sprintf("'%s' must be at least 1", name), call. = FALSE)
  }
  invisible(x)
}

assert_scalar_positive <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  if (x <= 0) {
    stop(sprintf("'%s' must be greater than 0", name), call. = FALSE)
  }
  invisible(x)
}

assert_scalar_unit_interval <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  if (x < 0 || x > 1) {
    stop(sprintf("'%s' must be between 0 and 1", name), call. = FALSE)
  }
  invisible(x)
}





`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}



read_csv <- function(path, ...) {
  read.csv(path, ..., stringsAsFactors = FALSE, check.names = FALSE)
}

mean_ci <- function(x, qs = c(0.025, 0.975), na.rm = TRUE) {
  c(mean = mean(x, na.rm = na.rm),
    quantile(x, probs = qs, na.rm = na.rm))
}

read_rds <- function(dir, combine = FALSE, id = NULL) {
  filepaths <- list.files(dir)
  filenames <- gsub(".rds", "", filepaths)
  results <- lapply(sprintf("%s/%s", dir, filepaths), readRDS)
  names(results) <- filenames
  
  if (combine) {
    results <- dplyr::bind_rows(results, .id = id)
  }
  
  results
}

named_list <- function(nms) {
  x <- vector("list", length(nms))
  names(x) <- nms
  x
}

