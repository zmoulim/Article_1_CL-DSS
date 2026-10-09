############################################################
# SECTION 7 - COMPARISON OF ULTIMATES BY ACCIDENT YEAR
# (point estimates only, no bootstrap)
#
# Run CL_DDS_corrected.R first in the same R session.
############################################################

library(ChainLadder)
library(clmplus)      # install.packages("clmplus")

if (!exists("diagnostic_global"))
  stop("Run CL_DDS_corrected.R first (diagnostic_global() not found).")

data(RAA); data(GenIns); data(MW2014)

############################################################
# Corrections with respect to the original script
#
# [1] "England & Verrall (GLM)": the original model was a Gamma GLM
#     (log link) fitted to CUMULATIVE amounts. This is not the model of
#     England & Verrall (2002). It gives 182 352 for RAA (-14.4%), the
#     value of the original table. The England & Verrall model is an
#     over-dispersed Poisson GLM fitted to INCREMENTAL amounts with
#     accident-year and development-year factors; its ultimates equal
#     the Chain-Ladder ultimates exactly (Renshaw & Verrall 1998).
#     odp_ultimates() fits this model (glmReserve() fails on RAA).
#
# [2] "Pittarello et al. (APC)": the original model was mgcv::gam with
#     smooth functions of age, calendar year and accident year on
#     cumulative amounts, which is not their model; the calendar spline
#     is extrapolated outside the observed years. Their model is
#     implemented in the package clmplus (hazard.model = "apc").
#     hazard.model = "a" must reproduce the Chain-Ladder (check).
#
# [3] "Mack": same ultimates as the Chain-Ladder by construction
#     (Mack's model only adds the prediction error). The column is
#     merged with the Chain-Ladder column.
#
# [4] "MDM-CL": replaced by the corrected CL-DDS.
#
# [5] "Maciak et al.": the original model, gam(Cum ~ s(Dev) + factor(AY)),
#     is not the functional profile method of Maciak, Mizera & Pesta
#     (2022), and it was fitted to cumulative amounts. It is replaced by
#     a GAM with smoothed development, NOT attributed to Maciak et al.:
#     over-dispersed Poisson on INCREMENTAL amounts, accident-year factor
#     and a smooth function of development age (the ODP GLM of [1] with
#     the development factors replaced by a spline). If the triangle
#     contains negative incremental amounts, the quasi-Poisson deviance
#     is not defined; these cells are then set to 0 for this model only
#     and their number is reported.
############################################################

library(mgcv)

# Triangle with negative incremental amounts set to 0 (cumulative amounts
# made non-decreasing). Used only for clmplus, which rejects recoveries.
nonneg_triangle <- function(tri){
  tri <- unclass(tri); m <- ncol(tri)
  inc <- tri; inc[, -1] <- tri[, -1] - tri[, -m]
  n_neg <- sum(inc < 0, na.rm = TRUE)
  inc[!is.na(inc) & inc < 0] <- 0
  out <- t(apply(inc, 1, function(x){ y <- cumsum(ifelse(is.na(x), 0, x)); y[is.na(x)] <- NA; y }))
  dimnames(out) <- dimnames(tri)
  attr(out, "n_negative") <- n_neg
  out
}

# ODP GLM fitted by IRLS (quasi-Poisson score equations). glmReserve()
# and glm(family = quasipoisson) fail when the triangle contains negative
# incremental amounts (RAA), because their starting values and deviance
# use log(y). The quasi-likelihood estimating equations remain valid with
# negative values, so the model is fitted directly.
odp_ultimates <- function(tri){
  tri <- unclass(tri); n <- nrow(tri); m <- ncol(tri)
  inc <- tri; inc[, -1] <- tri[, -1] - tri[, -m]
  d   <- data.frame(AY  = factor(rep(1:n, m)),
                    Dev = factor(rep(1:m, each = n)),
                    y   = as.vector(inc))
  X   <- model.matrix(~ AY + Dev, d)
  obs <- !is.na(d$y)
  Xo  <- X[obs, , drop = FALSE]; y <- d$y[obs]
  mu  <- rep(mean(y), length(y)); eta <- log(mu)
  for (it in 1:200) {
    fit     <- lm.wfit(Xo, eta + (y - mu) / mu, w = mu)
    eta_new <- drop(Xo %*% fit$coefficients)
    done    <- max(abs(eta_new - eta)) < 1e-10
    eta <- eta_new; mu <- exp(eta)
    if (done) break
  }
  pred   <- exp(drop(X %*% fit$coefficients))
  future <- matrix(ifelse(obs, 0, pred), n, m)
  latest <- apply(tri, 1, function(x) tail(x[!is.na(x)], 1))
  as.numeric(latest + rowSums(future))
}

gam_smoothed_ultimates <- function(tri){
  tri <- unclass(tri); n <- nrow(tri); m <- ncol(tri)
  inc <- tri; inc[, -1] <- tri[, -1] - tri[, -m]
  d <- data.frame(AY  = factor(rep(1:n, m)),
                  Dev = rep(1:m, each = n),
                  Inc = as.vector(inc))
  d <- d[!is.na(d$Inc), ]
  n_neg <- sum(d$Inc < 0)
  d$Inc <- pmax(d$Inc, 0)
  fit <- gam(Inc ~ AY + s(Dev, k = min(8, m - 1)),
             family = quasipoisson(link = "log"), data = d, method = "REML")
  fut <- expand.grid(AY = factor(1:n, levels = levels(d$AY)), Dev = 1:m)
  fut <- fut[is.na(inc[cbind(as.integer(fut$AY), fut$Dev)]), ]
  fut$pred <- predict(fit, newdata = fut, type = "response")
  latest <- apply(tri, 1, function(x) tail(x[!is.na(x)], 1))
  future <- tapply(fut$pred, factor(fut$AY, levels = 1:n), sum)
  future[is.na(future)] <- 0
  list(ultimates = as.numeric(latest + future), n_negative = n_neg,
       edf = sum(fit$edf[-(1:n)]))
}
############################################################

ultimates_by_origin <- function(tri, name){

  tri    <- as.triangle(unclass(tri))
  n      <- nrow(tri)
  latest <- as.numeric(getLatestCumulative(tri))

  # Chain-Ladder (= Mack point estimate)
  ult_cl <- as.numeric(summary(MackChainLadder(tri))$ByOrigin$Ultimate)

  # ODP GLM (England & Verrall)
  ult_odp <- odp_ultimates(tri)

  # APC (clmplus, Pittarello, Hiabu & Villegas)
  # clmplus does not accept decreasing cumulative amounts (recoveries).
  # As for the GAM, negative incremental amounts are set to 0 for this
  # model only; the reserve is computed on the adjusted triangle and
  # added to the latest observed amounts of the ORIGINAL triangle.
  tri_pos <- nonneg_triangle(tri)
  pp      <- AggregateDataPP(cumulative.payments.triangle = tri_pos, eta = 1/2)
  ult_a   <- tryCatch(
    latest + as.numeric(predict(clmplus(pp, hazard.model = "a"))$reserve),
    error = function(e){ message("clmplus 'a' failed for ", name, ": ",
                                 conditionMessage(e)); rep(NA, n) })
  ult_apc <- tryCatch(
    latest + as.numeric(predict(clmplus(pp, hazard.model = "apc"),
                                gk.fc.model = "a", ckj.fc.model = "a",
                                gk.order = c(1, 1, 0),
                                ckj.order = c(0, 1, 0))$reserve),
    error = function(e){ message("APC failed for ", name, ": ",
                                 conditionMessage(e)); rep(NA, n) })

  # GAM with smoothed development
  gs      <- gam_smoothed_ultimates(tri)
  ult_gam <- gs$ultimates

  # CL-DDS
  dds     <- diagnostic_global(unclass(tri), verbose = FALSE)
  ult_dds <- as.numeric(dds$ultimates)

  tab <- data.frame(AY            = 1:n,
                    Chain_Ladder  = ult_cl,
                    ODP_GLM       = ult_odp,
                    APC_clmplus   = ult_apc,
                    GAM_smoothed  = ult_gam,
                    CL_DDS        = ult_dds)
  tot  <- colSums(tab[, -1])
  diff <- round(100 * (tot / tot["Chain_Ladder"] - 1), 2)

  cat("\n==============================", name, "==============================\n")
  cat("Checks (must be ~0): ODP - CL =", round(tot["ODP_GLM"] - tot["Chain_Ladder"], 2),
      " | clmplus 'a' - CL =", round(sum(ult_a) - tot["Chain_Ladder"], 2), "\n")
  cat("Negative incremental cells (set to 0 for APC and GAM only) =",
      attr(tri_pos, "n_negative"), "\n")
  cat("GAM: negative incremental cells set to 0 =", gs$n_negative,
      " | effective df of the development spline =", round(gs$edf, 2), "\n")
  cat("CL-DDS decision:", dds$final_decision, " (K =", dds$diagnostic_zone, ")\n")
  print(dds$summary_tests)
  print(round(tab))
  cat("Total ultimate:\n");          print(round(tot))
  cat("Difference vs CL (%):\n");    print(diff)
  cat("Reserve (total ultimate - latest):\n"); print(round(tot - sum(latest)))

  invisible(list(table = tab, total = tot, diff = diff, dds = dds))
}

u_RAA    <- ultimates_by_origin(RAA,    "RAA")
u_GenIns <- ultimates_by_origin(GenIns, "GenIns")
u_MW2014 <- ultimates_by_origin(MW2014, "MW2014")
