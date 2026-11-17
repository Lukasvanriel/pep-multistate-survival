# ==============================================================================
# 25_I9b_placebo.R   (no INLA; resumable: re-run until it prints "all placebos done"; run from the project root)
# I9 (b) placebo: each plot's interval climate (z-GDD and z-CMI, jointly) shuffled across its own intervals, scenario A re-predicted with the shuffled
#   climate (22's engine, point parameters); F as in 22. Keeps each plot's climate range, breaks the timing: a real temporal signal -> slope ~ 0.
# 20 permutations; cached per permutation in figures/I/cache_I9.
# ==============================================================================
src <- readLines(here::here("R", "22_hindcast_space_for_time.R")); eval(parse(text = src[1:(grep("^plots <- unique", src) - 1)]))
CI9 <- file.path(FIG_DIR, "I", "cache_I9"); dir.create(CI9, recursive = TRUE, showWarnings = FALSE); N_PERM <- 20L
for (k in seq_len(N_PERM)) { f <- file.path(CI9, sprintf("perm_%02d.rds", k)); if (file.exists(f)) next
  if (as.numeric(Sys.time() - t_start, units = "secs") > 140) { cat("time budget reached; re-run\n"); quit(save = "no") }
  set.seed(900 + k); t0 <- Sys.time()
  ix <- iv[, .(row = .I, id_pe)][, perm := if (.N > 1) sample(row) else row, by = id_pe]$perm
  ps <- copy(iv); ps[, `:=`(zg = iv$zg[ix], zc = iv$zc[ix])]
  out <- rbindlist(lapply(S, function(s) { ii <- which(ps$state_start == s); if (!length(ii)) return(NULL); d <- ps[ii]
    m <- metrics(pred_end(s, d, 1L, d$zg, d$zc), rep(s, nrow(d))); setnames(m, paste0("P_", names(m))); cbind(row = ii, m) }))
  setorder(out, row); saveRDS(cbind(iv[, .(id_pe, interval)], out[, -1]), f)
  cat(sprintf("perm %d (%.0f s)\n", k, as.numeric(Sys.time() - t0, units = "secs"))) }
cat("all placebos done\n")
