############################################################
# bootstrap_cl_dds() - requires diagnostic_global() and prepare_dev_data()
# from CL_DDS_corrected.R (run that file first)
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
# EXAMPLE
############################################################
# data(RAA)
# bs_check <- bootstrap_cl_dds(RAA, B = 1000, force_classical = TRUE)
# c(reserve = bs_check$reserve, se = bs_check$se)   # se close to Mack 26 909
# bs_RAA <- bootstrap_cl_dds(RAA, B = 1000)                     # fixed model
# bs_RAA_sel <- bootstrap_cl_dds(RAA, B = 1000, reselect = TRUE) # re-selection
# bs_RAA_sel[c("selection", "reserve", "boot_mean", "se", "cv", "quantiles")]
