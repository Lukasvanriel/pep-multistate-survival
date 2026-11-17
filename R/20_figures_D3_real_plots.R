# ==============================================================================
# 20_figures_D3_real_plots.R   (no INLA; ~1-2 min; run from the project root)
# Exploratory set D3: first-move quantities for the REAL plots (own covariates) instead of a mean stand; averaged effects and maps. NOT manuscript figures.
# Per plot: state = state at its last visit; GDD, CMI, soil = mean over its visit intervals (raw), z-scored with the model scaling; time in state 0;
#   undisturbed baseline. Quantities within 50 y (0.5-y grid, competing exits, as 16): P(leave), expected thermal shift of the next state
#   (sum_j P(next = j) (STI_j - STI_state), STI as in 12).
# Perturbations, applied to every plot: GDD / CMI / soil +1 SD; each disturbance at entry (TSD = time in state; exits without the term at their
#   undisturbed rate). Averaged change = mean over the plots of (perturbed - baseline): over the plots of each state, and over all plots.
# Uncertainty: 50 joint draws (17's draws) on a random subsample of <= 300 plots per state (point estimate uses all plots).
# ==============================================================================
src <- readLines(here::here("R", "17_figures_C.R")); eval(parse(text = src[1:(grep("^# ---- propagation engine", src) - 1)]))
OUTD <- file.path(FIG_DIR, "D"); OUT <- OUTD; set.seed(3); H <- 50; DTG <- 0.5; TG <- seq(0, H, by = DTG); NG <- length(TG) - 1L
sc <- readRDS(file.path(FIT_DIR, "scaling.rds")); zf <- function(x, v) (x - sc[var == v, mean]) / sc[var == v, sd]
STI <- vi[, .(STI = mean(zf(gdd_3, "gdd_3"))), by = .(state = state_end)]; STI <- setNames(STI$STI, STI$state)
con <- DBI::dbConnect(RSQLite::SQLite(), here::here("data", "raw", "pep", "PEP.gpkg"), flags = RSQLite::SQLITE_RO)
xy <- as.data.table(DBI::dbGetQuery(con, "SELECT id_pe, latitude, longitude FROM placette")); DBI::dbDisconnect(con)
setorder(vi, id_pe, interval)
pl <- vi[, .(state = state_end[.N], zg = zf(mean(gdd_3, na.rm = TRUE), "gdd_3"), zc = zf(mean(cmi_gs, na.rm = TRUE), "cmi_gs"), zs = zf(mean(soil_pc1, na.rm = TRUE), "soil_pc1")), by = id_pe]
pl <- merge(pl, xy, by = "id_pe"); pl <- pl[state %in% S & !is.na(zg) & !is.na(zc) & !is.na(zs)]
cat(sprintf("D3: %d plots\n", nrow(pl)))
PERT <- list(baseline = list(), `GDD +1 SD` = list(dz = c(1, 0, 0)), `CMI +1 SD` = list(dz = c(0, 1, 0)), `soil PC1 +1 SD` = list(dz = c(0, 0, 1)),
             fire = list(dist = "fire"), pest = list(dist = "pest"), `light harvest` = list(dist = "harv_light"), `heavy harvest` = list(dist = "harv_heavy"))
# first-move quantities for the plots of one state, one parameter row r; returns matrix [plot, c(pleave, shift)]
first_move <- function(s, P, r, pert) { ex <- EX[[s]]; if (!length(ex$lb)) return(cbind(pleave = rep(0, nrow(P)), shift = 0))
  dz <- if (is.null(pert$dz)) c(0, 0, 0) else pert$dz; um <- (TG[-1] + TG[-length(TG)]) / 2
  dH <- lapply(ex$lb, function(lb) { X <- PAR[[lb]][r, , drop = FALSE]; a <- X[, "alpha"]
    eta <- X[, "b0y"] + gcol(X, "z_gdd_3") * (P$zg + dz[1]) + gcol(X, "z_cmi_gs") * (P$zc + dz[2]) + gcol(X, "z_soil_pc1") * (P$zs + dz[3])
    tsdp <- 0; if (!is.null(pert$dist) && TERMS[[pert$dist]][1] %in% colnames(X)) { tm <- TERMS[[pert$dist]]
      eta <- eta + X[, tm[1]] + gcol(X, tm[2]) * (P$zg + dz[1]); tsdp <- matrix(gcol(X, tm[3]) * f_tsd(um), nrow(P), NG, byrow = TRUE) }
    exp(outer(eta, rep(1, NG)) + tsdp) * matrix(diff(TG^a), nrow(P), NG, byrow = TRUE) })
  tot <- Reduce(`+`, dH); Hc <- t(apply(tot, 1, cumsum)); if (NG == 1) Hc <- t(Hc); Sm <- exp(-(cbind(0, Hc[, -NG, drop = FALSE]) + Hc) / 2)
  cif <- sapply(dH, function(m) rowSums(m * Sm)); if (!is.matrix(cif)) cif <- matrix(cif, nrow = nrow(P))
  pleave <- 1 - exp(-Hc[, NG]); nxt <- cif / pmax(pleave, 1e-12); dest <- S[ex$to]
  cbind(pleave = pleave, shift = drop(nxt %*% (STI[dest] - STI[S[s]]))) }
# point estimates, all plots
PT <- rbindlist(lapply(names(PERT), function(pn) rbindlist(lapply(unique(pl$state), function(st) { P <- pl[state == st]
  q <- first_move(match(st, S), P, 1L, PERT[[pn]]); data.table(id_pe = P$id_pe, state = st, pert = pn, pleave = q[, "pleave"], shift = q[, "shift"]) }))))
B0 <- PT[pert == "baseline", .(id_pe, pleave0 = pleave, shift0 = shift)]; PT <- merge(PT, B0, by = "id_pe")
PT[, `:=`(d_pleave = pleave - pleave0, d_shift = shift - shift0)]
fwrite(merge(PT, pl[, .(id_pe, latitude, longitude, zg, zc, zs)], by = "id_pe"), file.path(OUTD, "D3_plot_level_point.csv"))
# averaged effects with CI (subsample, draws)
sub <- pl[, .SD[sample(.N, min(.N, 300L))], by = state]
AV <- rbindlist(lapply(2:51, function(r) rbindlist(lapply(unique(sub$state), function(st) { P <- sub[state == st]; s <- match(st, S)
  b <- first_move(s, P, r, PERT$baseline)
  rbindlist(lapply(names(PERT)[-1], function(pn) { q <- first_move(s, P, r, PERT[[pn]])
    data.table(draw = r, state = st, pert = pn, n = nrow(P), d_pleave = mean(q[, 1] - b[, 1]), d_shift = mean(q[, 2] - b[, 2])) })) }))))
wst <- pl[, .N, by = state][, w := N / sum(N)]
AVall <- merge(AV, wst[, .(state, w)], by = "state")[, .(state = "all plots", d_pleave = sum(d_pleave * w), d_shift = sum(d_shift * w)), by = .(draw, pert)]
AV <- rbind(AV[, .(draw, state, pert, d_pleave, d_shift)], AVall)
PTm <- rbind(PT[pert != "baseline", .(d_pleave = mean(d_pleave), d_shift = mean(d_shift)), by = .(state, pert)],
             PT[pert != "baseline", .(state = "all plots", d_pleave = mean(d_pleave), d_shift = mean(d_shift)), by = pert])
SUM <- merge(PTm, AV[, .(pl_lo = quantile(d_pleave, .025), pl_hi = quantile(d_pleave, .975), sh_lo = quantile(d_shift, .025), sh_hi = quantile(d_shift, .975)), by = .(state, pert)], by = c("state", "pert"))
fwrite(SUM, file.path(OUTD, "D3_averaged_effects.csv"))
# compare with the mean-stand (z = 0) effect of B6 (same quantities: stay -> pleave, shift)
B6 <- fread(file.path(FIG_DIR, "B", "B6_climate_vs_disturbance.csv"))[quantity %in% c("stay", "shift")]
B6[, `:=`(est = fifelse(quantity == "stay", -est, est))]; B6 <- dcast(B6, origin + contrast ~ quantity, value.var = "est"); setnames(B6, c("origin", "contrast", "stay", "shift"), c("state", "pert", "zero_pleave", "zero_shift"))
SUMc <- merge(SUM, B6, by = c("state", "pert"), all.x = TRUE); fwrite(SUMc, file.path(OUTD, "D3_averaged_vs_mean_stand.csv"))

# ---- figures ------------------------------------------------------------------------------------------------------------------------------------------------
LV <- c("all plots", rev(ORDER)); SUM[, `:=`(state = factor(state, levels = rev(LV)), pert = factor(pert, levels = names(PERT)[-1]))]
PC <- c(`GDD +1 SD` = "#B2182B", `CMI +1 SD` = "#2166AC", `soil PC1 +1 SD` = "#8C510A", fire = "#D95F02", pest = "#7570B3", `light harvest` = "#66A61E", `heavy harvest` = "#1B6E3C")
mk <- function(v, lo, hi, xl) ggplot(SUM, aes(.data[[v]], state, colour = pert)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = .data[[lo]], xmax = .data[[hi]]), width = 0, orientation = "y", position = position_dodge(width = 0.75)) +
  geom_point(size = 1.5, position = position_dodge(width = 0.75)) + scale_colour_manual(values = PC, name = NULL) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) + labs(x = xl, y = NULL)
save_fig(mk("d_pleave", "pl_lo", "pl_hi", "change in P(leaving the state within 50 y)") + mk("d_shift", "sh_lo", "sh_hi", "change in expected ΔSTI of the next state") +
         plot_layout(guides = "collect") + plot_annotation(title = "D3a  Effects averaged over the real plots",
         subtitle = "each plot with its own climate and soil, state at its last visit; disturbance at entry; 95% CI"),
         "D3a_averaged_effects", 15, 8, NULL)
cmp <- melt(SUMc[state != "all plots" & !is.na(zero_pleave)], id.vars = c("state", "pert"), measure.vars = list(real = c("d_pleave", "d_shift"), zero = c("zero_pleave", "zero_shift")), variable.name = "q")
cmp[, q := c("P(leave within 50 y)", "ΔSTI of the next state")[q]]
p <- ggplot(cmp, aes(zero, real, colour = pert)) + geom_abline(linetype = 2, colour = "grey50") + geom_point() + scale_colour_manual(values = PC, name = NULL) +
  facet_wrap(~ q, scales = "free") + labs(x = "effect for a stand at mean climate (B6)", y = "effect averaged over the real plots (D3a)",
  title = "D3b  Mean stand vs real plots", subtitle = "one point per state and perturbation; dashed = 1:1")
save_fig(p, "D3b_mean_stand_vs_real_plots", 12, 5.5, NULL)
mp <- merge(PT[pert %in% c("baseline", "GDD +1 SD")], pl[, .(id_pe, latitude, longitude)], by = "id_pe")
wm <- as.data.table(ggplot2::map_data("world", region = c("Canada", "USA"))); bx <- range(mp$longitude) + c(-1, 1); by_ <- range(mp$latitude) + c(-0.7, 0.7)
mapp <- function(d, v, ttl, lims, pal) ggplot() + geom_polygon(data = wm, aes(long, lat, group = group), fill = "grey97", colour = "grey60", linewidth = 0.2) +
  geom_point(data = d[order(abs(get(v)))], aes(longitude, latitude, colour = .data[[v]]), size = 0.25) + pal + coord_quickmap(xlim = bx, ylim = by_) +
  labs(x = NULL, y = NULL, subtitle = ttl) + theme(panel.background = element_rect(fill = "#EAF2F8"), panel.grid = element_blank())
m1 <- mapp(mp[pert == "baseline"], "pleave", "P(leaving its current state within 50 y)", NULL, scale_colour_viridis_c(option = "magma", direction = -1, limits = c(0, 1), name = NULL))
m2 <- mapp(mp[pert == "baseline"], "shift", "expected ΔSTI of the next state", NULL, scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", limits = c(-1.5, 1.5), oob = scales::squish, name = NULL))
m3 <- mapp(mp[pert == "GDD +1 SD"], "d_pleave", "change in P(leave within 50 y) for +1 SD GDD", NULL, scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", limits = c(-0.3, 0.3), oob = scales::squish, name = NULL))
m4 <- mapp(mp[pert == "GDD +1 SD"], "d_shift", "change in ΔSTI of the next state for +1 SD GDD", NULL, scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", limits = c(-1, 1), oob = scales::squish, name = NULL))
save_fig((m1 + m2) / (m3 + m4) + plot_annotation(title = "D3c  Each plot's predicted first move (undisturbed)",
         subtitle = "own climate and soil, state at last visit; colours capped"), "D3c_maps", 15, 9, NULL)
cat("D3 done\n")
