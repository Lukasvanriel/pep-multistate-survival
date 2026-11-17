# ==============================================================================
# 25_I4_hysteresis_diagnostic.R   (no INLA; ~1 min; run from the project root). Diagnostic only: no state table is written for the fit.
# I4: how many observed state changes disappear under candidate hysteresis rules (user decides the rule afterwards).
# Rule tested (from the 2nd visit on; 1st visit classified as in 02): (a) dominance: the state changes only if the new dominant species leads the
#   best species of the current state by >= m of plot BA, else the previous state is kept; (b) regen: enter only if BA < lo, leave only if BA > hi
#   (instead of one threshold at 5 m2/ha). States and species mapping as in 02; retained visits only (state_classification.rds).
# ==============================================================================
suppressMessages({library(data.table); library(DBI)}); source(here::here("R", "config.R"))
src <- readLines(here::here("R", "02_classify_states.R")); i0 <- grep("^SPECIES_TO_STATE <- list", src); i1 <- i0 + which(src[(i0 + 1):length(src)] == ")")[1]
eval(parse(text = src[i0:i1])); map <- rbindlist(lapply(names(SPECIES_TO_STATE), function(s) data.table(essence = SPECIES_TO_STATE[[s]], st = s)))
sc <- readRDS(file.path(BUILD_DIR, "state_classification.rds"))
con <- dbConnect(RSQLite::SQLite(), "data/raw/pep/PEP.gpkg", flags = RSQLite::SQLITE_RO)
sp <- as.data.table(dbGetQuery(con, "SELECT id_pe, no_mes, essence, SUM(st_ha) AS ba FROM dendro_arbres WHERE etat IN ('10','12','30','32','40','42','50','52') GROUP BY id_pe, no_mes, essence"))
dbDisconnect(con); sp[, no_mes := as.numeric(no_mes)]; sp <- merge(sp, map, by = "essence", all.x = TRUE)
spl <- split(sp[paste(id_pe, no_mes) %in% sc[, paste(id_pe, no_mes)]], by = c("id_pe", "no_mes"), keep.by = FALSE)
setorder(sc, id_pe, no_mes); sc[, key := paste(id_pe, no_mes, sep = ".")]
grp <- function(s) fifelse(s == "regen", "regen", "established")
run <- function(m, lo, hi) {
  out <- character(nrow(sc)); ids <- sc$id_pe
  for (k in seq_len(nrow(sc))) {
    if (k == 1L || ids[k] != ids[k - 1L]) { out[k] <- sc$state[k]; next }
    prev <- out[k - 1L]; ba <- sc$total_ba[k]; cand <- sc$state[k]
    if (prev == "regen") { if (ba <= hi) { out[k] <- "regen"; next } }                    # stay regen unless BA > hi
    else if (ba < lo) { out[k] <- "regen"; next }                                          # enter regen only below lo
    if (cand == "regen") cand <- { x <- spl[[sc$key[k]]]; if (is.null(x) || !nrow(x)) prev else x[which.max(ba), st] }   # BA in [lo, 5): dominant species
    if (is.na(cand)) cand <- prev
    if (prev == "regen" || cand == prev) { out[k] <- cand; next }
    x <- spl[[sc$key[k]]]; tot <- sum(x$ba); lead <- x[st == cand, max(ba, 0)] - x[st == prev, max(c(ba, 0))]
    out[k] <- if (tot > 0 && lead / tot >= m) cand else prev }
  out }
base_chg <- sc[, .(from = shift(state), to = state), by = id_pe][!is.na(from) & from != to]
cat(sprintf("observed changes now: %d (between established states %d, into regen %d, out of regen %d)\n", nrow(base_chg),
    base_chg[from != "regen" & to != "regen", .N], base_chg[to == "regen", .N], base_chg[from == "regen", .N]))
res <- rbindlist(lapply(list(c(0.05, 4.5, 5.5), c(0.10, 4.5, 5.5), c(0.15, 4.5, 5.5), c(0.05, 4, 6), c(0.10, 4, 6), c(0.15, 4, 6)), function(r) {
  s <- run(r[1], r[2], r[3]); x <- data.table(id_pe = sc$id_pe, st = s)[, .(from = shift(st), to = st), by = id_pe][!is.na(from) & from != to]
  data.table(m = r[1], regen_band = sprintf("%.1f / %.1f", r[2], r[3]), changes = nrow(x), kept_pct = round(100 * nrow(x) / nrow(base_chg)),
             between = x[from != "regen" & to != "regen", .N], into_regen = x[to == "regen", .N], out_of_regen = x[from == "regen", .N],
             reversals = x[, .(a = from, b = to, c = shift(to, -1)), by = id_pe][!is.na(c) & c == a, .N]) }))
base_rev <- base_chg[, .(a = from, b = to, c = shift(to, -1)), by = id_pe][!is.na(c) & c == a, .N]
cat(sprintf("A -> B -> A reversals now: %d\n", base_rev)); print(res); fwrite(res, file.path(FIG_DIR, "I", "I4_hysteresis_diagnostic.csv"))
