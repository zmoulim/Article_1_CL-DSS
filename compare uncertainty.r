############################################################
# SECTION 8 - BOOTSTRAP / PREDICTION UNCERTAINTY
# Run CL_DDS_corrected.R first (diagnostic_global, bootstrap_cl_dds).
#
# Corrections with respect to the original code:
# [1] "Mack bootstrap" = BootChainLadder() is in fact the bootstrap of
#     England & Verrall (ODP GLM). It is now labelled "ODP bootstrap".
#     A separate "GLM bootstrap" is therefore not needed: the ODP GLM
#     reproduces the Chain-Ladder and BootChainLadder bootstraps it.
# [2] The original "GLM bootstrap" (Gamma GLM on cumulative amounts) and
#     "APC bootstrap" (mgcv splines) used models that are not those of
#     England & Verrall and Pittarello et al.; they are removed.
#     clmplus provides point estimates only (no bootstrap).
# [3] CV is now computed on the RESERVE (standard practice), not on the
#     ultimate, which contains the paid amounts and has no uncertainty.
# [4] Fixed seed, same B for all methods, and a validation run of the
#     CL-DDS bootstrap with the classical specification at every age,
#     whose standard error must be close to the Mack standard error.
############################################################

library(ChainLadder)

if (!exists("bootstrap_cl_dds"))
  stop("Run CL_DDS_corrected.R first (bootstrap_cl_dds() not found).")

data(RAA); data(GenIns); data(MW2014)

uncertainty_table <- function(tri, name, B = 1000, seed = 123){

  tri <- as.triangle(unclass(tri))

  # Mack (1993), analytical
  mack <- MackChainLadder(tri, est.sigma = "Mack")
  tot  <- summary(mack)$Totals
  R_cl <- tot["IBNR:", 1]; se_mack <- tot["Mack S.E.:", 1]

  # ODP bootstrap (England & Verrall)
  set.seed(seed)
  bt  <- BootChainLadder(tri, R = B, process.distr = "od.pois")
  ib  <- as.numeric(bt$IBNR.Totals)

  # CL-DDS bootstrap: validation (classical everywhere) and CL-DDS
  chk <- bootstrap_cl_dds(unclass(tri), B = B, seed = seed, force_classical = TRUE)
  dds <- bootstrap_cl_dds(unclass(tri), B = B, seed = seed)
  sel <- bootstrap_cl_dds(unclass(tri), B = B, seed = seed, reselect = TRUE)

  q <- function(x) quantile(x, c(0.95, 0.995))
  res <- rbind(
    data.frame(Method = "Mack (analytical)", Reserve = R_cl, Boot_mean = NA,
               SE = se_mack, Q95 = NA, Q995 = NA),
    data.frame(Method = "ODP GLM bootstrap (E&V)", Reserve = R_cl, Boot_mean = mean(ib),
               SE = sd(ib), Q95 = q(ib)[1], Q995 = q(ib)[2]),
    data.frame(Method = "CL-DDS bootstrap, classical (check)", Reserve = chk$reserve,
               Boot_mean = chk$boot_mean, SE = chk$se,
               Q95 = chk$quantiles[2], Q995 = chk$quantiles[3]),
    data.frame(Method = "CL-DDS bootstrap, fixed model", Reserve = dds$reserve,
               Boot_mean = dds$boot_mean, SE = dds$se,
               Q95 = dds$quantiles[2], Q995 = dds$quantiles[3]),
    data.frame(Method = "CL-DDS bootstrap, re-selection", Reserve = sel$reserve,
               Boot_mean = sel$boot_mean, SE = sel$se,
               Q95 = sel$quantiles[2], Q995 = sel$quantiles[3]))
  res$CV <- round(res$SE / res$Reserve, 4)
  res$Triangle <- name
  rownames(res) <- NULL

  cat("\n==============================", name, "==============================\n")
  cat("CL-DDS specification by age:", paste(dds$specification, collapse = " "), "\n")
  print(cbind(res[, 1, drop = FALSE], round(res[, 2:6]), CV = res$CV), row.names = FALSE)
  cat("Frequency of the specification selected at each tested age (re-selection):\n")
  print(round(sel$selection, 3))
  invisible(res)
}

unc_RAA    <- uncertainty_table(RAA,    "RAA")
unc_GenIns <- uncertainty_table(GenIns, "GenIns")
unc_MW2014 <- uncertainty_table(MW2014, "MW2014")
