############################################################
# SENSITIVITY OF THE CL-DDS TO THE DIAGNOSTIC PROPORTION
# Run CL_DDS_corrected.R first (diagnostic_global, bootstrap_cl_dds).
#
# Corrections with respect to sensi_p.R:
# [1] Old diagnostic_global (OLS, wrong projection) replaced.
# [2] M3IR5 converted from incremental to cumulative (incr2cum).
# [3] names(triangles) = names(triangles)[-c(13,14)] removed (name shift).
# [4] The quantity examined is no longer only the bootstrap CV (which was
#     computed on the ultimate with a faulty bootstrap) but, for each
#     proportion: K, rejection counts, DECISION, reserve change vs CL and
#     CV of the reserve (CL-DDS bootstrap with re-selection).
############################################################

library(ChainLadder)
if (!exists("bootstrap_cl_dds")) stop("Run CL_DDS_corrected.R first.")

tri_names <- c("RAA", "GenIns", "MW2014", "ABC", "M3IR5", "MCLpaid",
               "MCLincurred", "MW2008", "Mortgage", "UKMotor",
               "USAApaid", "USAAincurred")
props  <- c(0.3, 0.4, 0.5, 0.6, 0.7)
B_boot <- 500          # bootstrap replications per (triangle, proportion)
do_boot <- TRUE        # set FALSE for a quick run without bootstrap

get_tri <- function(name){
  data(list = name, package = "ChainLadder", envir = environment())
  tri <- unclass(as.triangle(get(name)))
  if (name == "M3IR5") tri <- unclass(incr2cum(as.triangle(tri)))
  tri
}

rows <- list()
for (name in tri_names) {
  tri    <- get_tri(name)
  latest <- apply(tri, 1, function(x) tail(x[!is.na(x)], 1))
  R_cl   <- sum(MackChainLadder(as.triangle(tri))$FullTriangle[, ncol(tri)]) - sum(latest)
  for (p in props) {
    cat(name, p, "\n")
    d <- tryCatch(diagnostic_global(tri, diagnostic_prop = p, verbose = FALSE),
                  error = function(e) NULL)
    if (is.null(d) || d$diagnostic_zone == 0) {
      rows[[length(rows) + 1]] <- data.frame(Triangle = name, prop = p, K = 0,
        R1 = NA, R2 = NA, R3 = NA, R4 = NA, Decision = NA,
        DDS_pct = NA, CV_DDS = NA); next
    }
    st <- d$summary_tests
    cv <- NA
    if (do_boot && R_cl > 0) {
      bs <- tryCatch(bootstrap_cl_dds(tri, B = B_boot, seed = 123, reselect = TRUE,
                                      diagnostic_prop = p), error = function(e) NULL)
      if (!is.null(bs)) cv <- round(100 * bs$se / bs$reserve, 1)
    }
    rows[[length(rows) + 1]] <- data.frame(
      Triangle = name, prop = p, K = d$diagnostic_zone,
      R1 = sum(st$H1_rej), R2 = sum(st$H2_rej),
      R3 = sum(st$H3_rej), R4 = sum(st$H4_rej),
      Decision = d$final_decision,
      DDS_pct = if (R_cl > 0) round(100 * (sum(d$ibnr) / R_cl - 1), 1) else NA,
      CV_DDS = cv)
  }
}
sens <- do.call(rbind, rows)

cat("\n===== Sensitivity: full table =====\n")
print(sens, row.names = FALSE)

cat("\n===== Decision by triangle and proportion =====\n")
print(reshape(sens[, c("Triangle", "prop", "Decision")], idvar = "Triangle",
              timevar = "prop", direction = "wide"), row.names = FALSE)

cat("\n===== CL-DDS reserve vs CL (%) =====\n")
print(reshape(sens[, c("Triangle", "prop", "DDS_pct")], idvar = "Triangle",
              timevar = "prop", direction = "wide"), row.names = FALSE)

cat("\n===== K (number of tested ages) =====\n")
print(reshape(sens[, c("Triangle", "prop", "K")], idvar = "Triangle",
              timevar = "prop", direction = "wide"), row.names = FALSE)

if (do_boot) {
  cat("\n===== CV of the CL-DDS reserve (%) =====\n")
  print(reshape(sens[, c("Triangle", "prop", "CV_DDS")], idvar = "Triangle",
                timevar = "prop", direction = "wide"), row.names = FALSE)
}

write.csv(sens, "sensitivity_results.csv", row.names = FALSE)
