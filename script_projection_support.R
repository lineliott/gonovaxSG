# These are functions adapted or cited from gonovax package in support of running script_projection.R 
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


create_diagnosis_waning_map <- function(n_vax, z, n_diag_rec = 1) {

  stopifnot(z > 0)
  stopifnot(n_vax %% n_diag_rec == 0)

  # set up waning map
  wd <- array(0, dim = c(n_vax, n_vax))

  #different base number of vaccine statuses (e.g. if X, V, W, then ntype = 3)
  ntype <- n_vax / n_diag_rec

  if (n_diag_rec >= 2) {
    for (k in 1:(ntype)) {
      for (j in 1:(n_diag_rec - 1)) {
        wd[(k - 1) * n_diag_rec + j, (k - 1) * n_diag_rec + j + 1] <- z
        wd[(k - 1) * n_diag_rec + j + 1, (k - 1) * n_diag_rec + j + 1] <- -z
      }

    }
  }

  wd
}


create_waning_map <- function(n_vax, i_v, i_w, z, n_diag_rec = 1) {

  stopifnot(z > 0)
  stopifnot(length(z) %in% c(1, length(i_v)))
  stopifnot(length(i_w) == n_diag_rec)
  # set up waning map
  w <- array(0, dim = c(n_vax, n_vax))

  for (i in seq_along(i_v)) {
    for (j in 1:n_diag_rec){

      w[i_w[j], i_v[(i - 1) * n_diag_rec + j]] <-
        ifelse(length(z) == 1,  z, z[i])

      w[i_v[(i - 1) * n_diag_rec + j], i_v[(i - 1) * n_diag_rec + j]] <-
        -w[i_w[j], i_v[(i - 1) * n_diag_rec + j]]
    }
  }

  w
}


create_uptake_map_xvwv <- function(n_group, n_vax, primary_uptake,
                                   booster_uptake, idx, n_diag_rec = 1,
                                   screening_or_diagnosis) {

  # set up uptake matrix rows = groups, columns = vaccine strata
  u <- array(0, dim = c(n_group, n_vax, n_vax))

  for (i in 1:n_diag_rec){

    if (screening_or_diagnosis == "screening") {
      temp <- i
    } else if (screening_or_diagnosis == "diagnosis") {

      if (i < n_diag_rec) {
        temp <- i + 1
      } else {
        temp <- i
      }
    } else {
      print("uptake map type not specified.")
    }

    u[, i, i] <- primary_uptake
    u[, idx$V[temp], i] <- primary_uptake

    u[, idx$W[i], idx$W[i]] <- booster_uptake
    u[, idx$V[temp], idx$W[i]] <- booster_uptake


  }

  u
}

set_strategy <- function(strategy = NULL, include_vbe = FALSE) {
  # switch for vaccination in group (L, H)
  novax  <- c(0, 0)
  vax_lh <- c(1, 1)
  vax_h  <- c(0, 1)

  if (is.null(strategy)) {
    vos <- vod <- novax

  } else if (strategy == "VoD") {
    vod <- vax_lh
    vos <- novax

  } else if (strategy == "VoA") {
    vod <- vos <- vax_lh

  } else if (strategy == "VoD(H)") {
    vod <- vax_h
    vos <- novax

  } else if (strategy == "VoA(H)") {
    vod <- vos <- vax_h

  } else if (strategy == "VoD+VoA(H)") {
    vod <- vax_lh
    vos <- vax_h

  } else if (strategy == "VoS") {
    vod <- novax
    vos <- vax_lh

  } else if (strategy == "VaH") {
    vod <- vax_lh
    vos <- vax_lh

  } else if (strategy == "VaHonly") {
    vod <- novax
    vos <- vax_lh

  } else {
    stop("strategy not recognised")
  }

  if (include_vbe) {
    vbe <- vax_lh
  } else {
    vbe <- novax
  }

  list(vod = vod, vos = vos, vbe = vbe)
}

create_vax_map_branching <- function(n_vax, v, i_e, i_p, set_vbe = FALSE, idx) {

  # ensure vaccine input is of correct length
  n_group <- 2
  n_vax <- idx$n_vax

  stopifnot(length(v) == n_group)
  stopifnot(all(v %in% c(0, 1)))


  if (length(i_e) > 0) {
    stopifnot(max(i_e, i_p) <= n_vax)
  }

  # set up vaccination matrix
  vax_map <- array(0, dim = c(n_group, n_vax, n_vax))

  if (set_vbe == TRUE) {

    vax_map[, idx$X[1], idx$X[1]] <-  v
    vax_map[, idx$V[1], idx$X[1]] <- -v

  } else {

    #repeat over stratum 1 column 1 for ease
    for (i in seq_along(i_e)) {

      vax_map[, i_e[i], i_e[i]] <-  v
      vax_map[, i_p[i], i_e[i]] <- -v

    }
  }

  vax_map
}


stratum_index_xvwv <- function(n_erlang = 1, n_diag_rec = 1, strategy = NULL) {

  # for an n_erlang of 3, and n_diag_rec of 2, the list of indexes returned
  # will be in the following order, where roman numerals refer to n_diag_rec,
  # and arabic numerals refer to erlang:
  # X.I, X.II, V1.I, V1.II, V2.I, V2.II, V3.I, V3.II, W.I, W.II

  ret <- list(X = seq_len(n_diag_rec))
  ret$V <- max(ret$X) + seq_len(n_erlang * n_diag_rec)
  ret$W <- max(ret$V) + seq_len(n_diag_rec)
  ret$n_vax <- max(ret$W)
  n_vax <- ret$n_vax

  ret$V1 <- ret$V[1:n_diag_rec]

  # strata people are diagnosed from
  ret$diagnosedfrom <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 0]

  # strata people are diagnosed to
  ret$diagnosedto <- seq_len(n_vax)[seq_len(n_vax) %% n_diag_rec  != 1]

  # strata people are vaccinated (before entry) from and to
  ret$vaccinatedfrom_vbe <- ret$X
  ret$vaccinatedto_vbe <- ret$V1

  # strata people are vaccinated (on screening) from and to
  if (!is.null(strategy)) {
    if (!is.null(strategy) && (strategy == "VaH" || strategy ==  "VaHonly")) {
      ret$vaccinatedfrom_vos <- c(ret$X[-1], ret$W[-1])
      ret$vaccinatedto_vos <- c(ret$V1[-1], ret$V1[-1])
    } else {
      ret$vaccinatedfrom_vos <- c(ret$X, ret$W)
      ret$vaccinatedto_vos <- c(ret$V1, ret$V1)
    }
  }

  #strata people are vaccinated (on diagnosis) from and to
  ret$vaccinatedfrom_vod <- c(ret$X, ret$W)
  ret$vaccinatedto_vod <- c(ret$V, ret$V)
  ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] <-
    ret$vaccinatedto_vod[ret$vaccinatedto_vod %% n_diag_rec != 0] + 1

  ret
}


vax_params_xvwv <- function(vea = 0, vei = 0, ved = 0, ves = 0,
                            dur = 1e3, uptake = 0, strategy = NULL,
                            vbe = 0, t_stop = 99, n_diag_rec = 1) {

  assert_scalar_unit_interval(vea)
  assert_scalar_unit_interval(vei)
  assert_scalar_unit_interval(ved)
  assert_scalar_unit_interval(ves)
  assert_scalar_positive(dur)
  assert_scalar_unit_interval(uptake)
  assert_scalar_unit_interval(vbe)
  assert_scalar_positive(t_stop)

  # waned vaccinees move to own stratum, but are eligible for re-vaccination
  # unvaccinated 1:n_diag_rec (x)
  # vaccinated from n_diag_rec+1 to 2*n_diag_rec (v)
  # waned from  2*n_diag_rec+1 to 3*n_diag_rec (w)

  # generate indices for all strata and
  idx <- stratum_index_xvwv(1, n_diag_rec = n_diag_rec, strategy = strategy)

  n_vax <- idx$n_vax

  i_v <- idx$V
  i_w <- idx$W

  n_group <- 2

  # create diagnosis history mapping
  diag_rec <- create_vax_map_branching(idx$n_vax, c(1, 1), idx$diagnosedfrom,
                                       idx$diagnosedto, set_vbe = FALSE, idx)

  # compartments to which vaccine efficacy applies
  ve <- ifelse(seq_len(n_vax) %in% idx$V, 1, 0)
  ved <- min(ved, 1 - 1e-10) # ensure duration is not divided by 0

  # If uptake of VbE > 0 consider that all adolescents are offered vaccine
  p <- set_strategy(strategy, vbe > 0)

  # set up uptake matrix rows = groups, columns = vaccine strata
  
  
  u_s <- create_uptake_map_xvwv(n_group = n_group, n_vax = n_vax,
                                primary_uptake = uptake,
                                booster_uptake = uptake,
                                idx, n_diag_rec = n_diag_rec,
                                screening_or_diagnosis = "screening")

  u_d <- create_uptake_map_xvwv(n_group = n_group, n_vax = n_vax,
                                primary_uptake = uptake,
                                booster_uptake = uptake,
                                n_diag_rec = n_diag_rec,
                                idx, screening_or_diagnosis = "diagnosis")
  # #eliott
  # u <- u_d
  
  if (sum(p$vod) > 0) {
    #vaccination on diagnosis occuring, so need to scale down diag_rec
    diag_rec[, idx$X, ] <- (1 - uptake) * diag_rec[, idx$X, ]
    diag_rec[, idx$W, ] <- (1 - uptake) * diag_rec[, idx$W, ]
  }

  willing <- rep(0, n_vax)
  willing[1] <- 1

  list(n_vax   = n_vax,
       willing = willing,
       # #eliott
       # u = u,
       u_s = u_s,
       u_d = u_d,
       u_vbe = vbe,
       vbe     = create_vax_map(n_vax, p$vbe, idx$vaccinatedfrom_vbe,
                                idx$vaccinatedto_vbe),
       vod     = create_vax_map(n_vax, p$vod, idx$vaccinatedfrom_vod,
                                idx$vaccinatedto_vod),
       vos     = create_vax_map(n_vax, p$vos, idx$vaccinatedfrom_vos,
                                idx$vaccinatedto_vos),
       vea = vea * ve,
       vei = vei * ve,
       ved = ved * ve,
       ves = ves * ve,
       w = create_waning_map(n_vax, i_v, i_w, 1 / dur, n_diag_rec),
       wd = create_diagnosis_waning_map(n_vax, 1, n_diag_rec),
       vax_t = c(0, t_stop),
       vax_y = c(1, 0),
       diag_rec = diag_rec
  )
}



run_onevax_xvwv <- function (tt, gono_params, demographic_params, init_params = NULL, dur = 1000, vea = 0,
                              vei = 0, ved = 0, ves = 0, vbe = 0, uptake = 0, strategy = NULL,
                              t_stop = 99)
  {
   if(!is.null(gono_params)) {stopifnot(all(lengths(list(uptake, vea, vei, ved, ves, dur)) %in%
                    c(1, length(gono_params))))}
    vax_params <- Map(vax_params_xvwv, uptake = uptake, dur = dur,
                      vea = vea, vei = vei, ved = ved, ves = ves, MoreArgs = list(strategy = strategy,
                                                                                  t_stop = t_stop, vbe = vbe))
    if (is.null(init_params)) {
      ret <- Map(gonovax::run, 
                 gono_params = gono_params, 
                 vax_params = vax_params,
                 MoreArgs = list(tt = tt,
                                 demographic_params = demographic_params_sg))
    }
    else {
      ret <- Map(gonovax::run, 
                 gono_params = gono_params,
                 init_params = init_params,
                 vax_params = vax_params,
                 MoreArgs = list(tt = tt,
                                 demographic_params = demographic_params_sg))
    }
    ret <- lapply(ret, name_outputs, c("X", "V", "W"))
    ret
}

name_outputs <- function(res, strata_names) {
  
  group_names <- c("L", "H")
  state_names <- c("U", "I", "A", "S", "T", "N",
                   "cum_incid", "cum_diag_a", "cum_diag_s",
                   "cum_treated", "cum_screened", "cum_vaccinated",
                   "cum_offered", "cum_vbe")
  
  for (nm in state_names) {
    dimnames(res[[nm]]) <- list(NULL, group_names, strata_names)
  }
  dimnames(res$lambda) <- list(NULL, group_names)
  dimnames(res$eta) <- list(NULL, group_names)
  
  res
}

run_grid_xvwv <- function (gono_params, init_params, cost_params, baseline, model, 
                          dur,uptake_total = 0, uptake_second_dose = uptake_total,
                          eff, vbe = 0, strategy = NULL, 
                          t_stop = 99, full_output = TRUE, disc_rate = 0)
{
  tt <- seq.int(init_params[[1]]$t, length.out = nrow((baseline[[1]][[1]])) +
                  1)
  
  # run model (set up for onevax_xvwv)
  l <- expand.grid(vea = eff, dur = dur)
  res <- furrr::future_pmap(.l = l, .f = model, tt = tt, gono_params = gono_params,
                            init_params = init_params, vbe = vbe, uptake = uptake_total,
                            strategy = strategy, t_stop = t_stop)
  # compare to baseline
  ret <- furrr::future_pmap(.l = list(y = res, baseline = baseline),
                            .f = gonovax::compare_baseline, cost_params = cost_params, uptake_first_dose = uptake_total/uptake_second_dose,
                            uptake_second_dose = uptake_second_dose, disc_rate = disc_rate)
  names(ret) <- sprintf("eff%.2f_dur%05.2f", l$vea, l$dur)
  # prepare results
  out <- list(inputs = list(t = tt, vbe = vbe, strategy = strategy,
                            grid = l, gono_params = gono_params, init_params = init_params,
                            cost_params = cost_params, disc_rate = disc_rate, uptake = list(total = uptake_total,
                                                                                            second_dose = uptake_second_dose)), results = ret)

  # compare to baseline
  ret <- furrr::future_pmap(.l = list(y = res, baseline = baseline),
                            .f = gonovax::compare_baseline, cost_params = cost_params, uptake_first_dose = uptake_total/uptake_second_dose,
                            uptake_second_dose = uptake_second_dose, disc_rate = disc_rate)
  names(ret) <- sprintf("eff%.2f_dur%05.2f", l$vea, l$dur)
  # prepare results
  out <- list(inputs = list(t = tt, vbe = vbe, strategy = strategy,
                            grid = l, gono_params = gono_params, init_params = init_params,
                            cost_params = cost_params, disc_rate = disc_rate, 
                            uptake = list(total = uptake_total,
                                          second_dose = uptake_second_dose)), 
              results = ret)
  if (full_output)
    out$full_results <- res
  class(out) <- "gonovax_grid"
  out
}


select_vaccine_profiles <- function(formation = "all") {
  efficacy <- c("low" = 0.2, "mid" = 0.4, "high" = 0.8)
  duration <- c("low" = 2, "mid" = 4, "high" = 8)
  
  if (formation == "plus") {
    vax_profiles <- data.frame(efficacy = c("mid", "mid", "mid", "low", "high"),
                               duration = c("mid", "low", "high", "mid", "mid"))
  } else if (formation == "cross") {
    vax_profiles <- data.frame(efficacy = c("mid", "low", "low", "high", "high"),
                               duration = c("mid", "low", "high", "low", "high"))
  } else if(formation == "all") {
    #vax_profiles <- expand_grid(duration = names(duration),
    vax_profiles <- tidyr::expand_grid(duration = names(duration),
                                efficacy = names(efficacy))
  }
  
  vax_profiles$name <- sprintf("eff%.2f_dur%02d",
                               efficacy[vax_profiles$efficacy],
                               duration[vax_profiles$duration])
  vax_profiles
}

gonovax_year_as_year_eliott <- function(gonovax_year) {
  assert_positive_integer(gonovax_year)
  2004 + gonovax_year
}

## moved run_grid() into package

tidy_output1 <- function(x, tt) {

  as.data.frame(x) %>% 
    dplyr::mutate(
      #year = gonovax::gonovax_year_as_year(tt),
                  year = gonovax_year_as_year_eliott(tt),
                  t = year - min(year)) %>% 
    #tidyr::pivot_longer(-c(year, t), "parameter_set") %>% 
    tidyr::pivot_longer(-c(year, t), names_to =  "parameter_set") %>% 
    dplyr::mutate(parameter_set = as.numeric(gsub("V", "", parameter_set)))
}

tidy_output_vaccine_profile <- function(vaccine_profile, tt, what = NULL) {

  what <- what %||% names(vaccine_profile)
  ret <- lapply(vaccine_profile[what], tidy_output1, tt = tt)
  dplyr::bind_rows(ret, .id = "what")
}

tidy_output <- function(res, what = NULL) {
  if (is.null(what)) {
    what <- c("treated", # annual cases treated (no reference baseline)
              "inc_cum_treated", # cumulative cases averted compared to baseline
              "inc_doses", # annual doses compared to baseline
              "inc_cum_doses", # cumulative doses compared to baseline
              "cet_20k",
              "cet_30k")
  }
  tt <- res$inputs$t[-1]
  ret <- lapply(res$results, tidy_output_vaccine_profile, tt, what)
  dplyr::bind_rows(ret, .id = "vaccine_profile")
}

tidy_summarise <- function(tidy_results) {
  tidy_results %>% 
    dplyr::group_by(dplyr::across(-c(parameter_set, value))) %>%
    dplyr::summarise(mean = mean(value, na.rm = TRUE),
                     median = median(value, na.rm = TRUE),
                     q2.5 = quantile(value, 0.025, na.rm = TRUE),
                     q5 = quantile(value, 0.05, na.rm = TRUE),
                     q10 = quantile(value, 0.10, na.rm = TRUE),
                     q97.5 = quantile(value, 0.975, na.rm = TRUE))
}

tidy_combine <- function(x) {
  dplyr::bind_rows(x, .id = "strategy")
}

rank_strategies <- function(forecasts) {

  ret <- forecasts %>%
    dplyr::group_by(dplyr::across(-c(value, strategy))) %>% 
    dplyr::mutate(rank = rank(value)) %>% 
    dplyr::select(-value) %>%  
    dplyr::ungroup() %>% 
    dplyr::group_by(dplyr::across(-c(parameter_set, strategy))) %>% 
    dplyr::count(strategy) %>% 
    tidyr::pivot_wider(names_from = c(rank, what),
                       names_prefix = "rank_",
                       values_from = n) %>% 
    dplyr::mutate(efficacy = as.numeric(substr(vaccine_profile, 4, 7)),
                  duration = as.numeric(substr(vaccine_profile, 12, 13)),
                  .after = vaccine_profile) %>%
    dplyr::mutate(strategy = factor(strategy, levels = unique(strategy))) %>% 
    dplyr::arrange(t, efficacy, duration, strategy) %>%
    dplyr::group_by(strategy, .add = TRUE)
  
  ret[is.na(ret)] <- 0
  
  ret
}



stat_boxplot_custom <- function(mapping = NULL, data = NULL,
                                geom = "boxplot", position = "dodge",
                                ...,
                                qs = c(0.025, 0.025, 0.5, 0.975, 0.975),
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatBoxplotCustom,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      qs = qs,
      ...
    )
  )
}


# boxplot customization
StatBoxplotCustom <- ggproto("StatBoxplotCustom", Stat,
                             required_aes = c("x", "y"),
                             non_missing_aes = "weight",
                             setup_params = function(data, params) {
                               params$width <- ggplot2:::"%||%"(
                                 params$width, (resolution(data$x) * 0.75)
                               )
                               
                               if (is.double(data$x) && !ggplot2:::has_groups(data) && any(data$x != data$x[1L])) {
                                 warning(
                                   "Continuous x aesthetic -- did you forget aes(group=...)?",
                                   call. = FALSE
                                 )
                               }
                               
                               params
                             },
                             
                             compute_group = function(data, scales, width = NULL, na.rm = FALSE, qs = c(.05, .25, 0.5, 0.75, 0.95)) {
                               
                               if (!is.null(data$weight)) {
                                 mod <- quantreg::rq(y ~ 1, weights = weight, data = data, tau = qs)
                                 stats <- as.numeric(stats::coef(mod))
                               } else {
                                 stats <- as.numeric(stats::quantile(data$y, qs))
                                 stats[3] <- mean(data$y) # replace median by mean
                               }
                               names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
                               iqr <- diff(stats[c(2, 4)])
                               
                               outliers <- (data$y < stats[1]) | (data$y > stats[5])
                               
                               if (length(unique(data$x)) > 1)
                                 width <- diff(range(data$x)) * 0.9
                               
                               
                               df <- as.data.frame(as.list(stats))
                               df$outliers <- list(data$y[outliers])
                               
                               if (is.null(data$weight)) {
                                 n <- sum(!is.na(data$y))
                               } else {
                                 # Sum up weights for non-NA positions of y and weight
                                 n <- sum(data$weight[!is.na(data$y) & !is.na(data$weight)])
                               }
                               
                               df$notchupper <- df$middle + 1.58 * iqr / sqrt(n)
                               df$notchlower <- df$middle - 1.58 * iqr / sqrt(n)
                               
                               df$x <- if (is.factor(data$x)) data$x[1] else mean(range(data$x))
                               df$width <- width
                               df$relvarwidth <- sqrt(n)
                               df
                             }
)

