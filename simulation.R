#!/usr/bin/env Rscript
## ============================================================================
## simulation.R  (v5)
## Variable Selection Using Model-X Knockoffs Under Informative Sampling
##
## Key design: ONE shared population (X, e, S_true, sign pattern fixed).
##   Only signal amplitude A varies, so Y_pop(A) = X %*% (A * beta_unit) + e.
##   This makes power curves smooth and monotone in A.
##
## Three informative sampling mechanisms:
##   "none"    : pi_i = expit(gamma0)
##   "error"   : pi_i = expit(gamma0 + gamma1 * e_i)
##   "outcome" : pi_i = expit(gamma0 + gamma1 * Y_i)
##
## Two methods: Oracle (IID subsample) vs Naive (informative sample)
##
## ------------------------------------------------------------
##   Rscript simulation.R --run --B 200 --N 5000 --cores 8
##
## Required packages: MASS, knockoff, glmnet
## ============================================================================

library(MASS)
library(knockoff)
library(glmnet)

## ---------------------------------------------------------------------------
## Utilities
## ---------------------------------------------------------------------------
expit <- function(x) 1 / (1 + exp(-x))

compute_fdp <- function(selected, S_true) {
  if (length(selected) == 0) return(0)
  sum(!(selected %in% S_true)) / length(selected)
}

compute_power <- function(selected, S_true) {
  if (length(S_true) == 0) return(1)
  sum(selected %in% S_true) / length(S_true)
}

## ---------------------------------------------------------------------------
## Calibrate gamma0
## ---------------------------------------------------------------------------
calibrate_gamma0 <- function(v, target_rate = 0.1) {
  ## v = gamma1 * (e or Y), without gamma0
  lo <- -50; hi <- 10
  for (iter in 1:100) {
    mid <- (lo + hi) / 2
    rate <- mean(expit(mid + v))
    if (rate > target_rate) hi <- mid else lo <- mid
  }
  g0 <- (lo + hi) / 2
  ## Warn if calibration hit boundary
  achieved <- mean(expit(g0 + v))
  if (abs(achieved - target_rate) / target_rate > 0.05)
    warning(sprintf("calibrate_gamma0: target=%.3f but achieved=%.3f (gamma0=%.2f)",
                    target_rate, achieved, g0))
  g0
}

## ---------------------------------------------------------------------------
## Generate the shared population skeleton (X, e, S_true, sign pattern)
## Y is computed later for each A
## ---------------------------------------------------------------------------
generate_skeleton <- function(N, p, s0, rho) {
  Sigma  <- rho^(abs(outer(1:p, 1:p, "-")))
  S_true <- sort(sample(1:p, s0))

  ## Unit coefficient vector with random signs (|beta_j| = 1 for j in S)
  beta_unit <- rep(0, p)
  beta_unit[S_true] <- sample(c(-1, 1), s0, replace = TRUE)

  X_pop <- mvrnorm(N, mu = rep(0, p), Sigma = Sigma)
  e_pop <- rnorm(N)

  list(X_pop = X_pop, e_pop = e_pop,
       beta_unit = beta_unit, Sigma = Sigma, S_true = S_true)
}

## ---------------------------------------------------------------------------
## Compute Y for a given amplitude A
## ---------------------------------------------------------------------------
make_Y <- function(skel, A) {
  as.numeric(skel$X_pop %*% (A * skel$beta_unit) + skel$e_pop)
}

## ---------------------------------------------------------------------------
## One Oracle replication
## ---------------------------------------------------------------------------
one_rep_oracle <- function(skel, Y_pop, n_oracle, q = 0.1) {
  N <- nrow(skel$X_pop)
  idx <- sample(N, n_oracle)
  X_or <- skel$X_pop[idx, , drop = FALSE]
  Y_or <- Y_pop[idx]

  res <- tryCatch({
    knockoff.filter(X_or, Y_or,
                    knockoffs = create.second_order,
                    statistic = stat.glmnet_coefdiff,
                    fdr = q, offset = 1)
  }, error = function(e) NULL)

  if (is.null(res))
    return(data.frame(method = "Oracle", fdp = NA_real_,
                      power = NA_real_, n = n_oracle,
                      stringsAsFactors = FALSE))
  sel <- res$selected
  data.frame(method = "Oracle",
             fdp   = compute_fdp(sel, skel$S_true),
             power = compute_power(sel, skel$S_true),
             n     = n_oracle, stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## One Naive replication
## ---------------------------------------------------------------------------
one_rep_naive <- function(skel, Y_pop, pi_vec, q = 0.1) {
  N <- nrow(skel$X_pop)
  p <- ncol(skel$X_pop)

  delta <- rbinom(N, 1, pi_vec)
  n_s   <- sum(delta)

  if (n_s < p + 10)
    return(data.frame(method = "Naive", fdp = NA_real_,
                      power = NA_real_, n = n_s,
                      stringsAsFactors = FALSE))

  S_idx <- which(delta == 1)
  X_s   <- skel$X_pop[S_idx, , drop = FALSE]
  Y_s   <- Y_pop[S_idx]

  res <- tryCatch({
    knockoff.filter(X_s, Y_s,
                    knockoffs = create.second_order,
                    statistic = stat.glmnet_coefdiff,
                    fdr = q, offset = 1)
  }, error = function(e) NULL)

  if (is.null(res))
    return(data.frame(method = "Naive", fdp = NA_real_,
                      power = NA_real_, n = n_s,
                      stringsAsFactors = FALSE))
  sel <- res$selected
  data.frame(method = "Naive",
             fdp   = compute_fdp(sel, skel$S_true),
             power = compute_power(sel, skel$S_true),
             n     = n_s, stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## Parallel dispatch helper
## ---------------------------------------------------------------------------
run_parallel <- function(worker, njobs, ncores, parallel_type,
                         extra_exports = character(0), verbose = FALSE) {
  if (ncores <= 1) {
    res_list <- vector("list", njobs)
    for (i in seq_len(njobs)) {
      res_list[[i]] <- worker(i)
      if (verbose && i %% 200 == 0)
        cat(sprintf("  [%d / %d]\n", i, njobs))
    }
    return(res_list)
  }
  if (parallel_type == "fork" && .Platform$OS.type != "windows") {
    return(parallel::mclapply(seq_len(njobs), worker,
                              mc.cores = ncores, mc.set.seed = FALSE))
  }
  ## PSOCK
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, varlist = extra_exports, envir = parent.frame())
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(MASS))
    suppressPackageStartupMessages(library(knockoff))
    suppressPackageStartupMessages(library(glmnet))
    NULL
  })
  parallel::parLapplyLB(cl, seq_len(njobs), worker)
}

## ---------------------------------------------------------------------------
## Smart core recommendation
## ---------------------------------------------------------------------------
recommend_cores <- function() {
  nc <- parallel::detectCores()
  if (.Platform$OS.type == "windows") max(1L, min(nc - 2L, 6L))
  else max(1L, min(nc - 1L, 96L))
}

## ---------------------------------------------------------------------------
## Summarise results
## ---------------------------------------------------------------------------
summarise_results <- function(results_long) {
  agg <- stats::aggregate(
    results_long[, c("fdp", "power")],
    by = results_long[, c("A", "scenario", "method"), drop = FALSE],
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(agg)[names(agg) == "fdp"]   <- "FDR"
  names(agg)[names(agg) == "power"] <- "Power"
  agg
}

## ---------------------------------------------------------------------------
## Main simulation runner
## ---------------------------------------------------------------------------
run_simulation <- function(
    N        = 5000,
    p        = 200,
    s0       = 50,
    rho      = 0.5,
    A_vec    = c(0.1, 0.2, 0.3, 0.5, 0.75, 1.0),
    scenarios = data.frame(
      pi_type = c("none", "error", "outcome"),
      gamma1  = c(0,      3,       1.0),
      label   = c("Non-informative",
                   "MNAR(e), g1=3",
                   "MNAR(Y), g1=1"),
      stringsAsFactors = FALSE
    ),
    q        = 0.1,
    n_oracle = 300,
    B        = 200,
    ncores   = recommend_cores(),
    parallel_type = c("auto", "psock", "fork"),
    verbose  = TRUE
) {
  parallel_type <- match.arg(parallel_type)
  n_sc <- nrow(scenarios)

  if (verbose) {
    cat("Knockoff simulation under informative sampling (v5)\n")
    cat("  N =", N, ", p =", p, ", s0 =", s0, ", rho =", rho, "\n")
    cat("  A_vec =", paste(A_vec, collapse = ", "), "\n")
    cat("  Scenarios:\n")
    for (i in seq_len(n_sc))
      cat(sprintf("    [%d] %s  (type=%s, gamma1=%.1f)\n",
                  i, scenarios$label[i], scenarios$pi_type[i], scenarios$gamma1[i]))
    cat("  n_oracle =", n_oracle, ", B =", B, ", Cores:", ncores, "\n\n")
  }

  ## --- Generate ONE shared skeleton (X, e, S_true, signs) ---
  if (verbose) cat("Generating shared population skeleton...\n")
  skel <- generate_skeleton(N, p, s0, rho)
  if (verbose) cat(sprintf("  active set = {%s}\n",
                           paste(skel$S_true, collapse = ",")))

  ## --- Pre-compute Y_pop for each A ---
  Y_list <- list()
  for (A in A_vec) {
    Y_list[[as.character(A)]] <- make_Y(skel, A)
  }

  ## --- Pre-compute inclusion probabilities for each (A, scenario) ---
  ## and calibrate gamma0
  pi_table  <- list()  # stores N-vectors of pi_i
  g0_table  <- list()  # stores gamma0 values
  for (A in A_vec) {
    Y_pop <- Y_list[[as.character(A)]]
    for (i in seq_len(n_sc)) {
      key <- paste(A, i, sep = "_")
      pt  <- scenarios$pi_type[i]
      g1  <- scenarios$gamma1[i]

      if (pt == "none") {
        v <- rep(0, N)
      } else if (pt == "error") {
        v <- g1 * skel$e_pop
      } else {
        v <- g1 * Y_pop
      }

      g0 <- calibrate_gamma0(v, target_rate = n_oracle / N)
      g0_table[[key]] <- g0
      pi_table[[key]] <- expit(g0 + v)

      if (verbose) {
        achieved_n <- sum(pi_table[[key]])
        cat(sprintf("  A=%.2f, %-20s => gamma0=%7.3f, E[n]=%.0f\n",
                    A, scenarios$label[i], g0, achieved_n))
      }
    }
  }

  ## Choose backend
  if (parallel_type == "auto")
    parallel_type <- if (.Platform$OS.type != "windows") "fork" else "psock"

  ## ============================================================
  ## Phase 1: Oracle (per A only, independent of scenario)
  ## ============================================================
  tasks_oracle <- expand.grid(rep = seq_len(B), A = A_vec,
                              stringsAsFactors = FALSE)
  if (verbose) cat(sprintf("\nPhase 1: %d Oracle jobs...\n", nrow(tasks_oracle)))

  worker_oracle <- function(i) {
    A     <- tasks_oracle$A[i]
    Y_pop <- Y_list[[as.character(A)]]
    res   <- one_rep_oracle(skel, Y_pop, n_oracle = n_oracle, q = q)
    res$A   <- A
    res$rep <- tasks_oracle$rep[i]
    res
  }

  oracle_list <- run_parallel(
    worker_oracle, nrow(tasks_oracle), ncores, parallel_type,
    extra_exports = c("expit", "compute_fdp", "compute_power",
                      "one_rep_oracle", "skel", "Y_list",
                      "tasks_oracle", "n_oracle", "q"),
    verbose = verbose)
  oracle_long <- do.call(rbind, oracle_list)

  ## Replicate Oracle across all scenarios
  oracle_all <- do.call(rbind, lapply(seq_len(n_sc), function(i) {
    tmp <- oracle_long; tmp$scenario <- scenarios$label[i]; tmp
  }))
  if (verbose) cat("  Oracle done.\n")

  ## ============================================================
  ## Phase 2: Naive (per A × scenario)
  ## ============================================================
  tasks_naive <- expand.grid(rep = seq_len(B), sc_id = seq_len(n_sc),
                             A = A_vec, stringsAsFactors = FALSE)
  if (verbose) cat(sprintf("Phase 2: %d Naive jobs...\n", nrow(tasks_naive)))

  worker_naive <- function(i) {
    A      <- tasks_naive$A[i]
    sc_id  <- tasks_naive$sc_id[i]
    Y_pop  <- Y_list[[as.character(A)]]
    pi_vec <- pi_table[[paste(A, sc_id, sep = "_")]]
    res    <- one_rep_naive(skel, Y_pop, pi_vec, q = q)
    res$A        <- A
    res$scenario <- scenarios$label[sc_id]
    res$rep      <- tasks_naive$rep[i]
    res
  }

  naive_list <- run_parallel(
    worker_naive, nrow(tasks_naive), ncores, parallel_type,
    extra_exports = c("expit", "compute_fdp", "compute_power",
                      "one_rep_naive", "skel", "Y_list", "pi_table",
                      "tasks_naive", "scenarios", "q"),
    verbose = verbose)
  naive_long <- do.call(rbind, naive_list)
  if (verbose) cat("  Naive done.\n")

  ## ============================================================
  ## Combine
  ## ============================================================
  results_long <- rbind(oracle_all, naive_long)
  summary <- summarise_results(results_long)

  list(
    results_long = results_long,
    summary      = summary,
    settings     = list(
      N = N, p = p, s0 = s0, rho = rho,
      A_vec = A_vec, scenarios = scenarios,
      q = q, n_oracle = n_oracle, B = B, ncores = ncores,
      S_true = skel$S_true
    )
  )
}

## ---------------------------------------------------------------------------
## Plotting
## ---------------------------------------------------------------------------
make_plots <- function(sim, filename = "simulation_results.pdf") {
  results <- sim$summary
  scenarios <- sim$settings$scenarios
  A_vals <- sort(unique(results$A))
  sc_labels <- scenarios$label
  n_sc <- length(sc_labels)

  pdf(filename, width = 4 * n_sc, height = 8)
  par(mfrow = c(2, n_sc), mar = c(4.5, 4.5, 3, 1),
      cex.lab = 1.2, cex.main = 1.1)

  methods <- c("Oracle", "Naive")
  cols  <- c(Oracle = "black", Naive = "red")
  pchs  <- c(Oracle = 16, Naive = 17)
  ltys  <- c(Oracle = 1, Naive = 2)

  ## Row 1: FDR
  for (sc in sc_labels) {
    plot(NA, xlim = range(A_vals), ylim = c(0, 0.3),
         xlab = "Signal amplitude A", ylab = "FDR", main = sc)
    abline(h = 0.1, lty = 3, col = "gray50", lwd = 1.5)
    for (m in methods) {
      sub <- results[results$scenario == sc & results$method == m, ]
      sub <- sub[order(sub$A), ]
      if (nrow(sub) > 0 && any(!is.na(sub$FDR))) {
        lines(sub$A, sub$FDR, col = cols[m], lty = ltys[m], lwd = 2)
        points(sub$A, sub$FDR, col = cols[m], pch = pchs[m], cex = 1.2)
      }
    }
    if (sc == sc_labels[1])
      legend("topright", legend = methods, col = cols[methods],
             lty = ltys[methods], pch = pchs[methods], lwd = 2, bty = "n")
  }

  ## Row 2: Power
  for (sc in sc_labels) {
    plot(NA, xlim = range(A_vals), ylim = c(0, 1),
         xlab = "Signal amplitude A", ylab = "Power", main = sc)
    for (m in methods) {
      sub <- results[results$scenario == sc & results$method == m, ]
      sub <- sub[order(sub$A), ]
      if (nrow(sub) > 0 && any(!is.na(sub$Power))) {
        lines(sub$A, sub$Power, col = cols[m], lty = ltys[m], lwd = 2)
        points(sub$A, sub$Power, col = cols[m], pch = pchs[m], cex = 1.2)
      }
    }
    if (sc == sc_labels[1])
      legend("bottomright", legend = methods, col = cols[methods],
             lty = ltys[methods], pch = pchs[methods], lwd = 2, bty = "n")
  }

  dev.off()
  cat("Plot saved to:", filename, "\n")
}

## ============================================================================
## CLI
## ============================================================================
print_usage <- function() {
  cat(paste0(
    "\nUsage:\n",
    "  Rscript simulation.R --run [options]\n\n",
    "Options:\n",
    "  --B <int>         Replications per setting (default 200)\n",
    "  --N <int>         Population size (default 5000)\n",
    "  --cores <int>     CPU cores (default: auto)\n",
    "  --parallel <str>  auto|psock|fork (default auto)\n",
    "  --out <path>      Output prefix (default knockoff_results)\n\n"
  ))
}

get_arg <- function(args, key, default = NULL) {
  pat <- paste0("^--", key, "(=|$)")
  hit <- grep(pat, args)
  if (length(hit) == 0) return(default)
  a <- args[hit[1]]
  if (grepl("=", a, fixed = TRUE)) return(sub(paste0("^--", key, "="), "", a))
  if (hit[1] < length(args)) return(args[hit[1] + 1])
  default
}

to_int <- function(x, default) {
  if (is.null(x) || is.na(x) || x == "") return(default)
  as.integer(x)
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!("--run" %in% args)) { print_usage(); quit(save = "no", status = 0) }

  B     <- to_int(get_arg(args, "B", "200"), 200)
  N     <- to_int(get_arg(args, "N", "5000"), 5000)
  default_cores <- recommend_cores()
  cores <- to_int(get_arg(args, "cores", as.character(default_cores)), default_cores)
  parallel_type <- get_arg(args, "parallel", "auto")
  out_path      <- get_arg(args, "out", "knockoff_results")

  t0 <- proc.time()

  sim <- run_simulation(
    B = B, N = N,
    ncores = cores,
    parallel_type = parallel_type,
    verbose = TRUE
  )

  elapsed <- proc.time() - t0
  cat(sprintf("\nElapsed: %.1f sec (%.1f min)\n",
              elapsed["elapsed"], elapsed["elapsed"] / 60))

  prefix <- sub("\\.rds$|\\.csv$", "", out_path, ignore.case = TRUE)
  utils::write.csv(sim$summary, file = paste0(prefix, "_summary.csv"), row.names = FALSE)
  utils::write.csv(sim$results_long, file = paste0(prefix, "_long.csv"), row.names = FALSE)
  cat("Saved:", paste0(prefix, "_summary.csv"), "\n")
  cat("Saved:", paste0(prefix, "_long.csv"), "\n")

  make_plots(sim, filename = paste0(prefix, "_plots.pdf"))

  cat("\nSummary:\n")
  print(sim$summary, row.names = FALSE)
}
