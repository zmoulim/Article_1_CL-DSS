############################################################
# SECTION 9 - EXTENSION TO ALL TRIANGLES OF THE ChainLadder PACKAGE
#
# Run first, in the same R session:
#   1. CL_DDS_corrected.R            (diagnostic_global, bootstrap_cl_dds)
#   2. compare_ultimates_corrected.R (odp_ultimates, nonneg_triangle,
#                                     gam_smoothed_ultimates)
#
# Corrections with respect to test_compl_boots.R:
# [1] names(triangles) = names(triangles)[-c(13,14)] assigned a SHORTER
#     vector of names to the list: the names were shifted / set to NA and
#     results could be attached to the wrong triangle. Triangles are now
#     listed explicitly by name.
# [2] Old diagnostic_global (OLS, wrong projection, raw residuals) replaced
#     by the corrected version.
# [3] Gamma GLM on cumulative amounts and mgcv "APC" replaced by the ODP GLM
#     (= Chain-Ladder) and clmplus, as in Section 7.
# [4] "Mack bootstrap" (BootChainLadder) relabelled ODP GLM bootstrap; the
#     MDM-CL bootstrap (log residuals, no parameter error, buggy factors)
#     replaced by bootstrap_cl_dds(reselect = TRUE).
# [5] All comparisons on the RESERVE, CV = SE / reserve.
# Each step is wrapped in tryCatch: a failure is reported, not hidden.
############################################################

library(ChainLadder)
library(clmplus)
library(mgcv)

needed <- c("diagnostic_global", "bootstrap_cl_dds", "odp_ultimates",
            "nonneg_triangle", "gam_smoothed_ultimates")
miss <- needed[!sapply(needed, exists)]
if (length(miss)) stop("Run CL_DDS_corrected.R and compare_ultimates_corrected.R first. Missing: ",
                       paste(miss, collapse = ", "))

tri_names <- c("RAA", "GenIns", "MW2014",
               "ABC", "M3IR5", "MCLpaid", "MCLincurred", "MW2008",
               "Mortgage", "UKMotor", "USAApaid", "USAAincurred")

B_boot <- 1000   # reduce to 500 if the run is too long

# Triangles stored as INCREMENTAL amounts in the package (see ?M3IR5:
# "Run off triangle of simulated incremental claims data"): converted to
# cumulative amounts before any analysis.
incremental_data <- c("M3IR5")

safe <- function(expr, what, name){
  tryCatch(expr, error = function(e){
    message("  [", name, "] ", what, " failed: ", conditionMessage(e)); NA })
}

one_triangle <- function(name){
  cat("\n######## ", name, " ########\n")
  data(list = name, package = "ChainLadder", envir = environment())
  tri <- unclass(as.triangle(get(name)))
  if (name %in% incremental_data) tri <- unclass(incr2cum(as.triangle(tri)))
  n <- nrow(tri); m <- ncol(tri)
  latest <- apply(tri, 1, function(x) tail(x[!is.na(x)], 1))
  inc <- tri; inc[, -1] <- tri[, -1] - tri[, -m]
  n_neg  <- sum(inc < 0, na.rm = TRUE)
  n_zero <- sum(tri[, -m] <= 0, na.rm = TRUE)

  # Diagnostics
  dds <- safe(diagnostic_global(tri, verbose = FALSE), "CL-DDS", name)
  ok  <- is.list(dds)
  st  <- if (ok) dds$summary_tests else NULL
  K   <- if (ok) dds$diagnostic_zone else NA
  R   <- if (ok && K > 0) colSums(st[, c("H1_rej","H2_rej","H3_rej","H4_rej")]) else rep(NA, 4)
  if (ok) { cat("K =", K, " decision:", dds$final_decision, "\n"); print(st, row.names = FALSE) }

  # Reserves
  mack  <- safe(MackChainLadder(as.triangle(tri), est.sigma = "Mack"), "Mack", name)
  R_cl  <- if (is.list(mack)) summary(mack)$Totals["IBNR:", 1] else NA
  se_mk <- if (is.list(mack)) summary(mack)$Totals["Mack S.E.:", 1] else NA
  R_odp <- safe(sum(odp_ultimates(tri)) - sum(latest), "ODP GLM", name)
  pp    <- safe(AggregateDataPP(cumulative.payments.triangle = nonneg_triangle(tri), eta = 1/2),
                "clmplus data", name)
  R_apc <- if (is.list(pp)) safe(sum(predict(clmplus(pp, hazard.model = "apc"),
                                gk.fc.model = "a", ckj.fc.model = "a",
                                gk.order = c(1, 1, 0), ckj.order = c(0, 1, 0))$reserve),
                                "APC", name) else NA
  R_gam <- safe(sum(gam_smoothed_ultimates(tri)$ultimates) - sum(latest), "GAM", name)
  R_dds <- if (ok) sum(dds$ibnr) else NA

  # Uncertainty
  bt <- safe(as.numeric(BootChainLadder(as.triangle(tri), R = B_boot,
                                        process.distr = "od.pois")$IBNR.Totals),
             "ODP bootstrap", name)
  bs <- if (ok) safe(bootstrap_cl_dds(tri, B = B_boot, reselect = TRUE),
                     "CL-DDS bootstrap", name) else NA

  pct <- function(x) round(100 * (x / R_cl - 1), 1)
  data.frame(
    Triangle = name, n = n, K = K,
    R1 = R[1], R2 = R[2], R3 = R[3], R4 = R[4],
    Decision = if (ok) dds$final_decision else NA,
    Neg_incr = n_neg, Nonpos_cells = n_zero,
    Reserve_CL = round(R_cl),
    ODP_pct = pct(R_odp), APC_pct = pct(R_apc), GAM_pct = pct(R_gam), DDS_pct = pct(R_dds),
    CV_Mack = round(se_mk / R_cl, 4),
    CV_ODPboot = if (is.numeric(bt) && length(bt) > 1) round(sd(bt) / R_cl, 4) else NA,
    CV_DDSboot = if (is.list(bs)) round(bs$se / R_dds, 4) else NA,
    DDSboot_mean_pct = if (is.list(bs)) round(100 * (bs$boot_mean / R_dds - 1), 1) else NA
  )
}

set.seed(123)
all_res <- do.call(rbind, lapply(tri_names, one_triangle))
rownames(all_res) <- NULL

cat("\n================ DIAGNOSTICS ================\n")
print(all_res[, c("Triangle","n","K","R1","R2","R3","R4","Decision","Neg_incr","Nonpos_cells")])
cat("\n================ RESERVES (difference vs CL, %) ================\n")
print(all_res[, c("Triangle","Reserve_CL","ODP_pct","APC_pct","GAM_pct","DDS_pct")])
cat("\n================ UNCERTAINTY (CV on reserve) ================\n")
print(all_res[, c("Triangle","CV_Mack","CV_ODPboot","CV_DDSboot","DDSboot_mean_pct")])

write.csv(all_res, "all_triangles_results.csv", row.names = FALSE)
