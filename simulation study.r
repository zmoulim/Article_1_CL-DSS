############################################################
# SIMULATION STUDY OF THE CL-DDS DECISION RULE
#
# Purpose
#  1. Size: how often does the scheme declare a triangle "structural" when
#     the Chain-Ladder model is true (delta = 0)?
#  2. Power: how often does it detect a real calendar-year effect (delta > 0)?
#  3. Accuracy: reserve error of the Chain-Ladder, the CL-DDS, the
#     age-period model of clmplus, and of the reserve selected by each
#     decision rule, against the TRUE simulated outstanding claims.
#  Expected results (Python check, 200-300 replications): under delta = 0
#  the per-age rejection rates of H1 and H2 are close to 5%; with onset
#  "mid", n = 14 and delta = 0.05 the scheme is structural in about 90% of
#  the replications, whereas n = 10 (K = 4) is never structural under the
#  rule K >= 6.
#  4. Calibration: the decision rule (K_min, share) is chosen on the
#     simulations (false structural rate <= 5% under delta = 0) and only
#     then applied, unchanged, to the real triangles.
#
# Run CL_DDS_corrected.R first (diagnostic_global, cl_factors).
############################################################

library(ChainLadder)
library(clmplus)

if (!exists("diagnostic_global")) stop("Run CL_DDS_corrected.R first.")

# ---------------- settings ----------------
n_sizes <- c(10, 14)               # K = 4 and K = 6
deltas  <- c(0, 0.02, 0.05, 0.10)  # strength of the calendar-year effect
onsets  <- c("mid", "late")        # when the effect starts
nsim    <- 200                     # replications per scenario (200 if too slow)
phi     <- 3                       # over-dispersion: Var = phi * mean
rules   <- data.frame(K_min = c(6,   4,   4,    4,   6),
                      share = c(0.5, 0.5, 0.75, 1.0, 0.75))
rule_names <- paste0("K>=", rules$K_min, ", share ", rules$share)

# ---------------- data generating process ----------------
# Mack-type recursive model (gamma), so that assumptions (M1)-(M3) hold
# exactly when delta = 0:
#   C_{i,1}   ~ Gamma(mean e_i P_1, CV 10%),  e_i = 1000 * 1.03^(i-1)
#   C_{i,j+1} | C_{i,j} ~ Gamma(mean m, variance phi (f_j - 1) C_{i,j}),
#   m = C_{i,j} * (1 + (f_j - 1) * exp(delta * max(0, t - t0))),
#   t = i + j (calendar period of cell (i, j+1)), f_j from a payment pattern.
# delta > 0: claims inflation on the incremental part, starting at calendar
#   period t0 and continuing in the future. The Chain-Ladder averages past
#   factors and therefore under-estimates the reserve.
# onset "mid":  t0 = n/2 (half of the observed history is affected);
# onset "late": t0 = n - 1 (only the last observed diagonal is affected):
#   a recent change, which no test on the observed data can detect.
# NOTE: an ODP (accident x development) DGP is NOT used for the size study:
# under that model E(C_{j+1} | C_j) is not proportional to C_j, so (M1)
# is false and the RESET test rejects correctly (about 20% per age).
simulate_triangle <- function(n, delta, onset = "mid"){
  t0 <- if (onset == "mid") floor(n / 2) else n - 1
  e  <- 1000 * 1.03^(0:(n - 1))
  p  <- (1:n)^1.5 * exp(-(1:n) / 2); P <- cumsum(p) / sum(p)
  f  <- P[-1] / P[-n]
  C  <- matrix(0, n, n)
  C[, 1] <- rgamma(n, shape = 100, scale = e * P[1] / 100)
  for (i in 1:n) for (j in 1:(n - 1)) {
    t <- i + j
    m <- C[i, j] * (1 + (f[j] - 1) * exp(delta * max(0, t - t0)))
    v <- phi * (f[j] - 1) * C[i, j]
    C[i, j + 1] <- rgamma(1, shape = m^2 / v, scale = v / m)
  }
  obs <- (outer(1:n, 1:n, "+") - 1) <= n
  tri <- C; tri[!obs] <- NA
  latest <- C[cbind(1:n, n:1)]
  list(tri = tri, true_reserve = sum(C[, n]) - sum(latest))
}

cl_reserve <- function(tri){
  f <- cl_factors(tri); n <- nrow(tri); m <- ncol(tri)
  latest <- apply(tri, 1, function(x) tail(x[!is.na(x)], 1))
  lastj  <- apply(tri, 1, function(x) max(which(!is.na(x))))
  ult <- sapply(1:n, function(i) if (lastj[i] < m) latest[i] * prod(f[lastj[i]:(m - 1)]) else latest[i])
  sum(ult) - sum(latest)
}

ap_reserve <- function(tri){
  pp <- AggregateDataPP(cumulative.payments.triangle = tri, eta = 1/2)
  sum(predict(clmplus(pp, hazard.model = "ap"),
              ckj.fc.model = "a", ckj.order = c(0, 1, 0))$reserve)
}

one_run <- function(n, delta, onset){
  s <- simulate_triangle(n, delta, onset)
  d <- tryCatch(diagnostic_global(s$tri, verbose = FALSE), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  st <- d$summary_tests
  data.frame(n = n, delta = delta, onset = onset, true = s$true_reserve,
             CL = cl_reserve(s$tri), DDS = sum(d$ibnr),
             AP = tryCatch(ap_reserve(s$tri), error = function(e) NA),
             K = d$diagnostic_zone,
             R1 = sum(st$H1_rej), R2 = sum(st$H2_rej),
             R3 = sum(st$H3_rej), R4 = sum(st$H4_rej))
}

# ---------------- run ----------------
set.seed(2026)
sims <- list()
for (on in onsets) for (n in n_sizes) for (dl in deltas) {
  if (on == "late" && dl == 0) next          # identical to onset "mid", delta 0
  cat("onset =", on, " n =", n, " delta =", dl, "\n")
  sims[[length(sims) + 1]] <- do.call(rbind, replicate(nsim, one_run(n, dl, on), simplify = FALSE))
}
sims <- do.call(rbind, sims)
write.csv(sims, "simulation_raw.csv", row.names = FALSE)

# ---------------- analysis ----------------
structural <- function(df, r) (df$K >= rules$K_min[r]) &
  (pmax(df$R1, df$R2) >= ceiling(rules$share[r] * df$K))

relerr <- function(est, true) (est - true) / true
summ <- function(e) c(bias = 100 * mean(e, na.rm = TRUE),
                      rmse = 100 * sqrt(mean(e^2, na.rm = TRUE)))

out_dec <- list(); out_err <- list()
for (on in onsets) for (n in n_sizes) for (dl in deltas) {
  df <- sims[sims$n == n & sims$delta == dl & (sims$onset == on | dl == 0), ]
  if (on == "late" && dl == 0) next
  dec <- sapply(seq_len(nrow(rules)), function(r) round(100 * mean(structural(df, r)), 1))
  out_dec[[length(out_dec) + 1]] <- data.frame(onset = on, n = n, delta = dl, t(setNames(dec, rule_names)),
                                               check.names = FALSE)
  e <- rbind(CL = summ(relerr(df$CL, df$true)),
             CL_DDS = summ(relerr(df$DDS, df$true)),
             AP = summ(relerr(df$AP, df$true)))
  for (r in seq_len(nrow(rules))) {
    sel <- ifelse(structural(df, r) & !is.na(df$AP), df$AP, df$DDS)
    e <- rbind(e, summ(relerr(sel, df$true)))
    rownames(e)[nrow(e)] <- paste("Rule:", rule_names[r])
  }
  out_err[[length(out_err) + 1]] <- data.frame(onset = on, n = n, delta = dl, method = rownames(e),
                                               round(e, 2), row.names = NULL)
}
dec_table <- do.call(rbind, out_dec)
err_table <- do.call(rbind, out_err)

cat("\n===== % of triangles declared STRUCTURAL (delta = 0: false structural rate) =====\n")
print(dec_table, row.names = FALSE)
cat("\n===== Reserve error vs TRUE outstanding (%): bias and RMSE =====\n")
print(err_table, row.names = FALSE)
cat("\n===== Average rejection counts =====\n")
print(aggregate(cbind(R1, R2, R3, R4) ~ onset + n + delta, data = sims, FUN = mean))
cat("\nAP model failures:", sum(is.na(sims$AP)), "of", nrow(sims), "\n")

write.csv(dec_table, "simulation_decisions.csv", row.names = FALSE)
write.csv(err_table, "simulation_errors.csv", row.names = FALSE)
