# R/simulation/sim04_diag.R (2026-11-01): diagnostic for sim04 hanging on N40_K0 (no replicate finished in 20+ min).
# Runs the 16 fits of N40_K0_rep01 one at a time, prints start/end + seconds, kills any INLA call after 120 s.
# Writes nothing to sim04's output. Run:  Rscript R/simulation/sim04_diag.R
suppressPackageStartupMessages({ library(data.table); library(INLA) })
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
ok <- tryCatch({ inla.setOption(inla.timeout = 120); TRUE }, error = function(e) FALSE)
cat("inla.timeout option available:", ok, "| INLA", as.character(packageVersion("INLA")), "\n")
ALPHA_TRUE <- 0.8; FOLLOW_UP <- 50; B0_TRUE <- log(log(2)) - ALPHA_TRUE * log(20)
PRIORS <- list(denis = list(a_param = 5, prec = 0.001, prec_icpt = 0.01), pipeline = list(a_param = 0.1, prec = 1 / 25, prec_icpt = 0))
seed <- 20261102L + 100000L * 3L + 10000L * 1L + 1L          # = JOBS seed of N40_K0_rep01
set.seed(seed); n <- 40
X <- matrix(rnorm(n * 10), n, 10); eta <- B0_TRUE + drop(X[, 1:3] %*% c(0.4, -0.3, 0.2))
T <- (-log(runif(n)) / exp(eta))^(1 / ALPHA_TRUE)
saved <- readRDS("data-output/simulation/sim04_variant_claim/datasets/N40_K0_rep01.rds")
cat("same times as sim04's saved dataset:", isTRUE(all.equal(T, saved$data$T)), "| events <= 50:", sum(T <= FOLLOW_UP), "\n")
d0 <- data.table(id = 1:n, T = T)
observe <- function(d, obs) { d <- copy(d)
  if (obs == "exact") return(d[, `:=`(status = as.integer(T <= FOLLOW_UP), time = pmin(T, FOLLOW_UP), time2 = pmin(T, FOLLOW_UP))][])
  w <- as.numeric(sub("ic", "", obs)); grid <- sort(unique(c(seq(0, FOLLOW_UP, by = w), FOLLOW_UP)))
  ev <- d$T <= FOLLOW_UP; U <- grid[findInterval(d$T, grid, left.open = TRUE) + 1L]; L <- grid[findInterval(d$T, grid, left.open = TRUE)]
  d[, `:=`(status = fifelse(ev, 3L, 0L), time = fifelse(ev, L, FOLLOW_UP), time2 = fifelse(ev, U, FOLLOW_UP))][] }
for (ob in c("exact", "ic1", "ic3", "ic5")) { d <- observe(d0, ob)
  cat(sprintf("\n[%s] status counts: %s | min time among events: %.3f\n", ob, paste(names(table(d$status)), table(d$status), sep = "=", collapse = " "), min(d$time[d$status != 0])))
  f <- if (ob == "exact") inla.surv(time = d$time, event = d$status) ~ 1 else inla.surv(time = d$time, time2 = d$time2, event = d$status) ~ 1
  for (v in 0:1) for (p in names(PRIORS)) { pr <- PRIORS[[p]]
    cat(sprintf("%s  start %s v%d %-8s ... ", format(Sys.time(), "%H:%M:%S"), ob, v, p)); t1 <- proc.time()[[3]]
    m <- tryCatch(inla(f, family = "weibullsurv", data = d,
                       control.family = list(variant = v, hyper = list(alpha = list(prior = "pc.alphaw", param = pr$a_param))),
                       control.fixed = list(prec = pr$prec, mean.intercept = 0, prec.intercept = pr$prec_icpt),
                       control.inla = list(strategy = "adaptive", int.strategy = "grid", control.vb = list(enable = FALSE)), num.threads = "1:1"),
                  error = function(e) list(error = conditionMessage(e)))
    a <- if (!is.null(m$error)) paste("ERROR:", substr(m$error, 1, 100)) else if (isTRUE(m$ok == FALSE)) "ok=FALSE" else sprintf("alpha %.3f (sd %.3f)", m$summary.hyperpar[1, 1], m$summary.hyperpar[1, 2])
    cat(sprintf("%.1f s | %s\n", proc.time()[[3]] - t1, a)) } }
