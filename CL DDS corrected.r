############################################################
# CL-DDS : Chain-Ladder Diagnostic Decision-Support scheme
# Corrected version of the original diagnostic_global() code.
# Every change is marked [CORR n]; see the list at the end of the file.
############################################################

library(ChainLadder)
library(lmtest)
# [CORR 0] packages tseries, car and dplyr removed: they were not used.

############################################################
# PREPARE DEVELOPMENT DATA
############################################################
# [CORR 1] Uses only the observed pairs (C_{i,j}, C_{i,j+1}), so that the
#          function also works on truncated triangles (rolling backtest).
prepare_dev_data <- function(triangle, j){
  tri  <- unclass(triangle)
  rows <- which(!is.na(tri[, j]) & !is.na(tri[, j + 1]))
  data.frame(
    AY       = rows,
    Dev      = j,
    Calendar = rows + j - 1,
    Ci_j     = tri[cbind(rows, j)],
    Ci_j1    = tri[cbind(rows, j + 1)]
  )
}

############################################################
# CLASSICAL CHAIN-LADDER FACTORS (volume-weighted)
############################################################
# [CORR 2] Computed directly as sum(C_{i,j+1}) / sum(C_{i,j}); identical to
#          MackChainLadder()$f but also valid on truncated triangles.
cl_factors <- function(triangle){
  tri <- unclass(triangle)
  sapply(1:(ncol(tri) - 1), function(j){
    d <- prepare_dev_data(tri, j)
    sum(d$Ci_j1) / sum(d$Ci_j)
  })
}

############################################################
# MAIN DIAGNOSTIC FUNCTION
############################################################
diagnostic_global <- function(triangle,
                              diagnostic_prop        = 0.5,
                              alpha                  = 0.05,
                              min_obs                = 5,
                              min_ages_structural    = 4,   # calibrated by simulation_study.R
                              share_structural       = 0.5,
                              calendar_extrapolation = c("last", "trend"),
                              verbose                = TRUE){

  # [CORR 3] Unused argument global_threshold removed; the decision-rule
  #          parameters (6 ages, 50%) are now explicit arguments.
  calendar_extrapolation <- match.arg(calendar_extrapolation)

  tri <- unclass(triangle)
  n   <- nrow(tri)
  m   <- ncol(tri)

  ##########################################################
  # DIAGNOSTIC ZONE
  ##########################################################
  K <- floor(diagnostic_prop * (m - 1))
  # [CORR 4] Each tested age must have at least min_obs observations
  #          (binds only for very short triangles, e.g. in the backtest).
  while (K > 0 && nrow(prepare_dev_data(tri, K)) < min_obs) K <- K - 1

  if (verbose) {
    cat("\n====================================================\n")
    cat(" CL-DDS : diagnostic zone K =", K, "of", m - 1, "development ages\n")
    cat("====================================================\n")
  }

  summary_results <- data.frame()
  models          <- list()
  t_max           <- numeric(0)

  ##########################################################
  # LOOP ON DEVELOPMENT AGES
  ##########################################################
  for (j in seq_len(K)) {

    df <- prepare_dev_data(tri, j)
    t_max[j] <- max(df$Calendar)

    #######################################################
    # H1 : CONDITIONAL MEAN (RESET)
    #######################################################
    # [CORR 5] The Chain-Ladder regression is a WLS regression through the
    #          origin with weights 1/C_{i,j}; its coefficient is exactly the
    #          Chain-Ladder factor. The original code used OLS (no weights).
    model_cl <- lm(Ci_j1 ~ 0 + Ci_j, data = df, weights = 1 / Ci_j)

    # [CORR 6] RESET computed on the weighted regression (lmtest::resettest
    #          ignores weights). With a single regressor, powers of the fitted
    #          values span the same space as powers of C_{i,j}; C_{i,j} is
    #          rescaled for numerical stability.
    df$z2 <- (df$Ci_j / mean(df$Ci_j))^2
    df$z3 <- (df$Ci_j / mean(df$Ci_j))^3
    model_reset <- lm(Ci_j1 ~ 0 + Ci_j + z2 + z3, data = df, weights = 1 / Ci_j)
    p_reset <- anova(model_cl, model_reset)$`Pr(>F)`[2]

    positive <- all(df$Ci_j > 0 & df$Ci_j1 > 0)
    if (p_reset < alpha && positive) {
      current_model <- lm(log(Ci_j1) ~ log(Ci_j), data = df)   # specification (7')
      model_type    <- "LOG"
    } else {
      current_model <- model_cl
      model_type    <- "CLASSICAL"
    }
    h1_rej <- p_reset < alpha

    #######################################################
    # H2 : CALENDAR-YEAR EFFECT (F-test, nested models)
    #######################################################
    if (model_type == "CLASSICAL") {
      # [CORR 7] Multiplicative calendar term, consistent with eq. (9):
      #          E(C_{i,j+1}) = (f_j + gamma_j t) C_{i,j}. The original code
      #          added Calendar additively and without weights.
      model_exog <- lm(Ci_j1 ~ 0 + Ci_j + Ci_j:Calendar, data = df, weights = 1 / Ci_j)
    } else {
      model_exog <- lm(log(Ci_j1) ~ log(Ci_j) + Calendar, data = df)
    }
    p_exog <- anova(current_model, model_exog)$`Pr(>F)`[2]
    h2_rej <- p_exog < alpha
    if (h2_rej) current_model <- model_exog

    #######################################################
    # STANDARDIZED RESIDUALS ON THE ORIGINAL SCALE
    #######################################################
    # [CORR 8] r = (C_{i,j+1} - fitted) / sqrt(C_{i,j}). Under (M2) these
    #          residuals have constant variance; the raw residuals do not.
    if (model_type == "LOG") {
      s2           <- summary(current_model)$sigma^2
      fitted_level <- exp(fitted(current_model) + s2 / 2)   # lognormal mean
    } else {
      fitted_level <- fitted(current_model)
    }
    dfr <- data.frame(r = (df$Ci_j1 - fitted_level) / sqrt(df$Ci_j), Ci_j = df$Ci_j)

    #######################################################
    # H3 : VARIANCE STRUCTURE (Breusch-Pagan, studentized)
    #######################################################
    # [CORR 9] Applied to the standardized residuals, variance regressed on
    #          C_{i,j}. The original code applied bptest() to the auxiliary
    #          regression lm(res2 ~ Ci_j) built on raw OLS residuals, which
    #          tests homoscedasticity, not assumption (M2).
    bp_test <- bptest(r ~ 1, varformula = ~ Ci_j, data = dfr)
    h3_rej  <- bp_test$p.value < alpha

    #######################################################
    # H4 : INDEPENDENCE (Durbin-Watson)
    #######################################################
    # [CORR 10] Applied to the standardized residuals, ordered by accident
    #           year; alternative: positive first-order autocorrelation.
    dw_test <- dwtest(r ~ 1, data = dfr)
    h4_rej  <- dw_test$p.value < alpha

    if (verbose) {
      cat(sprintf("Age %d : %-9s  H1 p=%.4f  H2 p=%.4f  H3 p=%.4f  H4 p=%.4f\n",
                  j, model_type, p_reset, p_exog, bp_test$p.value, dw_test$p.value))
    }

    summary_results <- rbind(summary_results, data.frame(
      Development = j,
      Model       = model_type,
      H1_p = round(p_reset, 4), H2_p = round(p_exog, 4),
      H3_p = round(bp_test$p.value, 4), H4_p = round(dw_test$p.value, 4),
      H1_rej = h1_rej, H2_rej = h2_rej, H3_rej = h3_rej, H4_rej = h4_rej
    ))
    models[[j]] <- current_model
  }

  ##########################################################
  # DECISION RULE (Section 5.3)
  ##########################################################
  rownames(summary_results) <- NULL
  n_H1 <- if (K > 0) sum(summary_results$H1_rej) else 0
  n_H2 <- if (K > 0) sum(summary_results$H2_rej) else 0
  structural <- (K >= min_ages_structural) &&
                (max(n_H1, n_H2) >= ceiling(share_structural * K))
  decision <- if (structural) "STRUCTURAL" else "LOCAL"

  rejection_summary <- data.frame(
    Hypothesis     = c("H1", "H2", "H3", "H4"),
    Rejection_rate = if (K > 0) colMeans(summary_results[, c("H1_rej", "H2_rej", "H3_rej", "H4_rej")]) else rep(NA, 4)
  )
  rownames(rejection_summary) <- NULL

  if (verbose) {
    print(rejection_summary)
    cat("Decision :", decision, "\n")
  }

  ##########################################################
  # PROJECTION
  ##########################################################
  f_cl <- cl_factors(tri)

  # [CORR 11] Projection uses the model retained at each tested age, cell by
  #           cell. In the original code:
  #           - the log model selected under H1 was never used (factor reset
  #             to the classical one), so H1 had no effect on the reserve;
  #           - an OLS coefficient replaced the Chain-Ladder factor even when
  #             no hypothesis was rejected (e.g. RAA age 4: 1.1717 -> 1.1620);
  #           - with a calendar effect, the coefficient of C_{i,j} in a model
  #             with an additive Calendar term was used as a factor
  #             (e.g. RAA age 2: 1.6235 -> 1.1396).
  # [CORR 12] Calendar effect in the future: "last" (default) keeps it at its
  #           last observed level for that age; "trend" extrapolates the
  #           linear trend. The choice must be stated in the paper.
  predict_next <- function(j, c_val, t_val){
    if (j > K) return(f_cl[j] * c_val)
    mod <- models[[j]]
    if (calendar_extrapolation == "last") t_val <- min(t_val, t_max[j])
    nd <- data.frame(Ci_j = c_val, Calendar = t_val)
    if (summary_results$Model[j] == "CLASSICAL") {
      as.numeric(predict(mod, newdata = nd))
    } else {
      s2 <- summary(mod)$sigma^2
      as.numeric(exp(predict(mod, newdata = nd) + s2 / 2))
    }
  }

  proj <- tri
  for (i in 1:n) {
    for (j in 1:(m - 1)) {
      if (is.na(proj[i, j + 1]) && !is.na(proj[i, j])) {
        proj[i, j + 1] <- predict_next(j, proj[i, j], i + j - 1)
      }
    }
  }

  latest_observed <- apply(tri, 1, function(x) tail(na.omit(x), 1))
  ultimates       <- proj[, m]
  ibnr            <- ultimates - latest_observed

  # Implied factor for the first future cell of each age (for reporting)
  implied <- sapply(1:(m - 1), function(j){
    i  <- max(which(!is.na(tri[, j])))
    c0 <- tri[i, j]
    predict_next(j, c0, i + j - 1) / c0
  })
  factor_summary <- data.frame(Development = 1:(m - 1),
                               Classical   = round(f_cl, 4),
                               Corrected   = round(implied, 4))

  if (verbose) {
    print(factor_summary)
    cat("Total ultimate :", round(sum(ultimates), 2), "\n")
    cat("Total reserve  :", round(sum(ibnr), 2), "\n")
  }

  list(diagnostic_zone   = K,
       summary_tests     = summary_results,
       rejection_rates   = rejection_summary,
       final_decision    = decision,
       models            = models,
       f_cl              = f_cl,
       predict_next      = predict_next,
       factor_summary    = factor_summary,
       projected_triangle = proj,
       ultimates         = ultimates,
       ibnr              = ibnr)
}

############################################################
# ROLLING-ORIGIN BACKTEST
############################################################
rolling_origin_backtest <- function(triangle, min_train = 6, ...){
  tri <- unclass(triangle)
  n   <- nrow(tri)
  out <- data.frame()

  for (k in min_train:(n - 1)) {
    # [CORR 13] Training data = diagonals 1..k only. The original code used
    #           triangle[1:k, 1:k], which contains cells of later calendar
    #           years, including the cells used as test values (data leakage).
    train <- tri[1:k, 1:k]
    for (i in 1:k) for (j in 1:k) if (i + j - 1 > k) train[i, j] <- NA

    fit <- diagnostic_global(train, verbose = FALSE, ...)

    # Test cells: one step ahead, on diagonal k + 1
    for (r in 2:k) {
      jj <- k + 1 - r                     # C[r, jj] lies on diagonal k
      if (jj + 1 <= ncol(tri) && !is.na(tri[r, jj + 1])) {
        true_value <- tri[r, jj + 1]
        pred_cl    <- fit$f_cl[jj] * tri[r, jj]
        pred_corr  <- fit$predict_next(jj, tri[r, jj], r + jj - 1)
        out <- rbind(out, data.frame(
          Iteration = k, AY = r, Dev = jj, True = true_value,
          Classical = pred_cl, Corrected = pred_corr,
          Error_Classical = abs(pred_cl   - true_value) / true_value,
          Error_Corrected = abs(pred_corr - true_value) / true_value,
          Decision  = fit$final_decision))
      }
    }
  }

  differ <- abs(out$Classical - out$Corrected) > 1e-8
  res <- list(details              = out,
              n_cells              = nrow(out),
              n_cells_corrected    = sum(differ),
              mean_classical_error = mean(out$Error_Classical),
              mean_corrected_error = mean(out$Error_Corrected),
              rmse_classical       = sqrt(mean(out$Error_Classical^2)),
              rmse_corrected       = sqrt(mean(out$Error_Corrected^2)))
  cat("\nROLLING ORIGIN :", res$n_cells, "cells,", res$n_cells_corrected, "corrected\n")
  cat("Mean abs. error  CL:", round(100 * res$mean_classical_error, 2), "%  CL-DDS:",
      round(100 * res$mean_corrected_error, 2), "%\n")
  cat("RMSE             CL:", round(res$rmse_classical, 4), "   CL-DDS:",
      round(res$rmse_corrected, 4), "\n")
  res
}

############################################################
# SENSITIVITY OF THE DECISION TO THE DIAGNOSTIC PROPORTION
############################################################
# For Section 11: reports the decision (not only the CV) for each proportion.
sensitivity_prop <- function(triangle, props = seq(0.3, 0.7, by = 0.1), ...){
  do.call(rbind, lapply(props, function(p){
    fit <- diagnostic_global(triangle, diagnostic_prop = p, verbose = FALSE, ...)
    data.frame(prop     = p,
               K        = fit$diagnostic_zone,
               n_H1     = sum(fit$summary_tests$H1_rej),
               n_H2     = sum(fit$summary_tests$H2_rej),
               decision = fit$final_decision,
               ultimate = round(sum(fit$ultimates), 0))
  }))
}

############################################################
# BOOTSTRAP OF THE CL-DDS RESERVE
############################################################
# Conditional residual bootstrap in the spirit of the Mack bootstrap:
#  - the specification retained by the CL-DDS at each age (classical,
#    classical + calendar, log, log + calendar) is kept FIXED;
#  - standardized residuals e = (C_{i,j+1} - fitted) / (sigma_j sqrt(C_{i,j})),
#    corrected for degrees of freedom, are pooled over all ages (they all
#    have unit variance, unlike the raw ratios of the original code);
#  - parameter error: pseudo-values C*_{i,j+1} = fitted + sigma_j sqrt(C_{i,j}) e*
#    (conditional on the observed C_{i,j}); each age is re-estimated;
#  - process error: future cells drawn from a gamma distribution with mean
#    given by the re-estimated model and variance sigma_j^2 C_{i,j} (M2);
#  - calendar effect kept at its last observed level ("last").
# Check: with force_classical = TRUE the bootstrap standard error must be
# close to the Mack standard error (RAA: Mack S.E. = 26 909).
# reselect = TRUE: in each replication the tests of H1 and H2 are run
#   again on the pseudo-data of each tested age, and the specification is
#   chosen by the same rules as diagnostic_global(). The bootstrap then
#   includes the uncertainty of the model choice. (The structural/local
#   decision is not re-evaluated: it is LOCAL for the three triangles and
#   far from the threshold.) The frequency with which each specification
#   is selected is returned in $selection.
bootstrap_cl_dds <- function(triangle, B = 1000, seed = 123,
                             force_classical = FALSE, reselect = FALSE,
                             alpha = 0.05, ...){
  set.seed(seed)
  tri <- unclass(triangle); n <- nrow(tri); A <- ncol(tri) - 1

  dds <- diagnostic_global(tri, verbose = FALSE, alpha = alpha, ...)
  K   <- dds$diagnostic_zone
  spec <- rep("CL", A)
  if (!force_classical && K > 0) for (j in 1:K) {
    spec[j] <- paste0(if (dds$summary_tests$Model[j] == "LOG") "LOG" else "CL",
                      if (dds$summary_tests$H2_rej[j]) "CAL" else "")
  }
  npar <- c(CL = 1, CLCAL = 2, LOG = 2, LOGCAL = 3)
  dat  <- lapply(1:A, function(j) prepare_dev_data(tri, j))
  tmax <- sapply(dat, function(d) max(d$Calendar))

  fit_age <- function(s, d){
    switch(s,
      CL     = list(s = s, f = sum(d$Ci_j1) / sum(d$Ci_j)),
      CLCAL  = list(s = s, b = coef(lm(Ci_j1 ~ 0 + Ci_j + Ci_j:Calendar,
                                       data = d, weights = 1 / Ci_j))),
      LOG    = { md <- lm(log(Ci_j1) ~ log(Ci_j), data = d)
                 list(s = s, b = coef(md), s2 = summary(md)$sigma^2) },
      LOGCAL = { md <- lm(log(Ci_j1) ~ log(Ci_j) + Calendar, data = d)
                 list(s = s, b = coef(md), s2 = summary(md)$sigma^2) })
  }
  mu_age <- function(p, c, t){
    lc <- log(pmax(c, 1e-12))
    switch(p$s,
      CL     = p$f * c,
      CLCAL  = (p$b[1] + p$b[2] * t) * c,
      LOG    = ifelse(c > 0, exp(p$b[1] + p$b[2] * lc + p$s2 / 2), 0),
      LOGCAL = ifelse(c > 0, exp(p$b[1] + p$b[2] * lc + p$b[3] * t + p$s2 / 2), 0))
  }

  # Same selection rules as diagnostic_global() (RESET, then calendar F-test)
  select_spec <- function(d){
    m_cl <- lm(Ci_j1 ~ 0 + Ci_j, data = d, weights = 1 / Ci_j)
    d$z2 <- (d$Ci_j / mean(d$Ci_j))^2
    d$z3 <- (d$Ci_j / mean(d$Ci_j))^3
    m_re <- lm(Ci_j1 ~ 0 + Ci_j + z2 + z3, data = d, weights = 1 / Ci_j)
    p1   <- anova(m_cl, m_re)$`Pr(>F)`[2]
    use_log <- isTRUE(p1 < alpha) && all(d$Ci_j > 0 & d$Ci_j1 > 0)
    if (use_log) {
      m0 <- lm(log(Ci_j1) ~ log(Ci_j), data = d)
      m1 <- lm(log(Ci_j1) ~ log(Ci_j) + Calendar, data = d)
    } else {
      m0 <- m_cl
      m1 <- lm(Ci_j1 ~ 0 + Ci_j + Ci_j:Calendar, data = d, weights = 1 / Ci_j)
    }
    p2 <- anova(m0, m1)$`Pr(>F)`[2]
    paste0(if (use_log) "LOG" else "CL", if (isTRUE(p2 < alpha)) "CAL" else "")
  }

  # Fit, sigma_j^2 and pooled standardized residuals
  pars <- lapply(1:A, function(j) fit_age(spec[j], dat[[j]]))
  sig2 <- rep(NA_real_, A); pool <- c()
  for (j in 1:A) {
    d <- dat[[j]]
    e <- (d$Ci_j1 - mu_age(pars[[j]], d$Ci_j, d$Calendar)) / sqrt(d$Ci_j)
    dfree <- nrow(d) - npar[[spec[j]]]
    if (dfree > 0) {
      sig2[j] <- sum(e^2) / dfree
      pool    <- c(pool, e / sqrt(sig2[j]) * sqrt(nrow(d) / dfree))
    }
  }
  for (j in 1:A) if (is.na(sig2[j]))        # Mack (1993) extrapolation
    sig2[j] <- min(sig2[j - 1]^2 / sig2[j - 2], sig2[j - 2], sig2[j - 1])
  pool <- pool - mean(pool)

  latest <- apply(tri, 1, function(x) tail(x[!is.na(x)], 1))
  lastj  <- apply(tri, 1, function(x) max(which(!is.na(x))))

  project <- function(pp, process){
    tot <- 0
    for (i in 1:n) {
      cv <- latest[i]
      if (lastj[i] <= A) for (j in lastj[i]:A) {
        mval <- mu_age(pp[[j]], cv, min(i + j - 1, tmax[j]))
        if (!process)      cv <- mval
        else if (mval <= 0) cv <- 0
        else {
          v  <- sig2[j] * cv
          cv <- if (v > 0) rgamma(1, shape = mval^2 / v, scale = v / mval) else mval
        }
      }
      tot <- tot + cv
    }
    tot - sum(latest)
  }

  R0   <- project(pars, FALSE)
  sims <- numeric(B)
  chosen <- matrix(NA_character_, B, max(K, 1))
  for (b in 1:B) {
    pp <- vector("list", A)
    for (j in 1:A) {
      d   <- dat[[j]]
      mu0 <- mu_age(pars[[j]], d$Ci_j, d$Calendar)
      for (k in 1:100) {
        ys <- mu0 + sqrt(sig2[j] * d$Ci_j) * sample(pool, nrow(d), replace = TRUE)
        if (!grepl("LOG", spec[j]) || all(ys > 0)) break
      }
      ds <- d; ds$Ci_j1 <- ys
      s_b <- spec[j]
      if (reselect && !force_classical && j <= K) {
        s_b <- tryCatch(select_spec(ds), error = function(e) spec[j])
        chosen[b, j] <- s_b
      }
      pp[[j]] <- fit_age(s_b, ds)
    }
    sims[b] <- project(pp, TRUE)
  }

  selection <- NULL
  if (reselect && !force_classical && K > 0)
    selection <- t(sapply(1:K, function(j)
      table(factor(chosen[, j], levels = c("CL", "CLCAL", "LOG", "LOGCAL"))) / B))

  list(specification = spec,
       selection     = selection,
       reserve       = R0,
       boot_mean     = mean(sims),
       se            = sd(sims),
       cv            = sd(sims) / R0,
       quantiles     = quantile(sims, c(0.75, 0.95, 0.995)),
       distribution  = sims)
}

############################################################
# RUN
############################################################
data(RAA); data(GenIns); data(MW2014)

res_RAA <- diagnostic_global(RAA)
res_RAA$summary_tests
res_RAA$factor_summary

res_GenIns <- diagnostic_global(GenIns)
res_MW2014 <- diagnostic_global(MW2014)

bt_RAA    <- rolling_origin_backtest(RAA,    min_train = 6)
bt_GenIns <- rolling_origin_backtest(GenIns, min_train = 6)
bt_MW2014 <- rolling_origin_backtest(MW2014, min_train = 8)

sensitivity_prop(RAA)

# Bootstrap: validation first (must be close to Mack S.E. = 26 909 for RAA)
bs_check <- bootstrap_cl_dds(RAA, B = 1000, force_classical = TRUE)
c(reserve = bs_check$reserve, se = bs_check$se)
bs_RAA <- bootstrap_cl_dds(RAA, B = 1000)
bs_RAA[c("specification", "reserve", "boot_mean", "se", "cv", "quantiles")]
# with re-selection of the model in each replication
bs_RAA_sel <- bootstrap_cl_dds(RAA, B = 1000, reselect = TRUE)
bs_RAA_sel[c("selection", "reserve", "boot_mean", "se", "cv", "quantiles")]

############################################################
# REFERENCE VALUES (independent Python computation, RAA, default settings)
# Compare with res_RAA before using any result:
#   K = 4
#   Age 1 : CLASSICAL  H1 p = 0.0947   H2 p = 0.7607   H3 p = 0.4438
#   Age 2 : CLASSICAL  H1 p = 0.3982   H2 p = 0.0319   H3 p = 0.2860
#   Age 3 : LOG        H1 p = 4.8e-05  H2 p = 0.8095   H3 p = 0.5381
#   Age 4 : CLASSICAL  H1 p = 0.3091   H2 p = 0.9747   H3 p = 0.7565
#   H4 (Durbin-Watson) statistics d = 2.697, 2.513, 1.849, 2.109
#   Decision : LOCAL
#   Classical factors : 2.9994 1.6235 1.2709 1.1717 1.1134 1.0419 1.0333 1.0169 1.0092
#   Corrected (first future cell) : 2.9994 2.0847 1.2487 1.1717 then classical
#   Total ultimate : 222 979 ("last")   225 585 ("trend")   CL : 213 122
#   Bootstrap (B = 2000-3000, 3 seeds, Python):
#     force_classical : reserve 52 135, s.e. 26 900 - 28 000 (Mack 26 909)
#     CL-DDS          : reserve 61 992, boot mean 58 500 - 59 100,
#                       s.e. 22 600 - 24 600, CV 0.37 - 0.40
#     CL-DDS, reselect = TRUE (B = 1000, 2 seeds):
#                       boot mean 58 400 - 59 200, s.e. 30 800 - 31 800,
#                       CV 0.50 - 0.51; calendar term at age 2 selected in
#                       about 63-64% of replications, log model at age 3 in
#                       about 13-15%
#   (different random generator: R values will differ slightly)
############################################################

############################################################
# LIST OF CORRECTIONS
#  0  unused packages removed
#  1  prepare_dev_data() uses observed pairs only
#  2  classical factors computed directly (volume-weighted)
#  3  unused global_threshold removed; decision parameters explicit
#  4  at least min_obs observations per tested age
#  5  Chain-Ladder regression estimated by WLS (weights 1/C)
#  6  RESET computed on the weighted regression
#  7  calendar term multiplicative and weighted (eq. 9)
#  8  standardized residuals on the original scale
#  9  Breusch-Pagan on standardized residuals
# 10  Durbin-Watson on standardized residuals
# 11  projection uses the retained model (log model and calendar term)
# 12  explicit choice for the future calendar effect ("last" / "trend")
# 13  backtest without data leakage (training = diagonals 1..k)
############################################################
