# ==============================================================================
# 13_fig_helpers.R -- shared helpers + inputs for the exploratory figure scripts 13_figures_A.R and 14_figures_A5.R. Sourced, not run.
# Defines ORDER, PAL, lab(), save_fig(), heat(), fmt_hr(), cap_sym(), hr_breaks(); reads pooled, hr, ar, iy, ts, vi; builds modelled, ev, trow.
# OUT = data-output/<tag>/figures/A/
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(here) })
source(here::here("R", "config.R"))   # RUN_TAG, TABLE_DIR, RESULTS_DIR, FIG_DIR
OUT <- file.path(FIG_DIR, "A"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ORDER <- c("regen", "jack_pine", "white_spruce", "black_spruce", "fir", "white_birch", "trembling_aspen", "other_pioneer",
           "yellow_birch", "red_maple", "sugar_maple", "cedar", "other_temperate")
PAL <- c(regen = "#777777", jack_pine = "#9E9A2B", white_spruce = "#6BAED6", black_spruce = "#08306B", fir = "#238B45",
         white_birch = "#A6761D", trembling_aspen = "#FD8D3C", other_pioneer = "#E6B800", yellow_birch = "#E7298A",
         red_maple = "#CB181D", sugar_maple = "#7A0177", cedar = "#66C2A4", other_temperate = "#9E9AC8")
lab <- function(x) gsub("_", " ", x)
COL_NOTMOD <- "grey92"; COL_NOTERM <- "grey70"
theme_set(theme_bw(base_size = as.numeric(Sys.getenv("FIG_BASE_SIZE", "9"))) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")))

save_fig <- function(p, name, w, h, data = NULL) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  if (!is.null(data)) fwrite(data, file.path(OUT, paste0(name, ".csv")))
  cat(sprintf("  wrote %s (.png/.pdf%s)\n", name, if (is.null(data)) "" else "/.csv"))
}
split_label <- function(dt) dt[, c("from", "to") := tstrsplit(label, "->", fixed = TRUE)][]
full_grid <- function() CJ(from = ORDER, to = ORDER)[from != to]
cap_sym <- function(x, q = 0.98) { m <- quantile(abs(x[is.finite(x)]), q, na.rm = TRUE); max(m, log(1.5)) }   # CONFIRM
fmt_hr <- function(x) ifelse(x >= 100, sprintf("%.0f", x), ifelse(x >= 10, sprintf("%.1f", x), sprintf("%.2f", x)))

# heatmap of origin (rows) x destination (columns); d needs from, to, status in {value, noterm, notmod}, fillv, txt, bold
HR_BREAKS <- c(0.01, 0.03, 0.1, 0.25, 0.5, 1, 2, 4, 10, 30, 100)
hr_breaks <- function(lim) log(HR_BREAKS[abs(log(HR_BREAKS)) <= lim + 1e-9])
heat <- function(d, fill_name, limits = NULL, diverging = TRUE, title = NULL, subtitle = NULL, trans_lab = waiver(), breaks = waiver(), seq_dir = -1) {
  d <- copy(d); d[, from := factor(from, levels = rev(ORDER))]; d[, to := factor(to, levels = ORDER)]
  rng <- if (is.null(limits)) range(d[status == "value", fillv], finite = TRUE) else limits
  d[, rel := if (diverging) abs(fillv) / max(abs(rng)) else (fillv - rng[1]) / diff(rng)]
  if (!diverging && seq_dir == 1) d[, rel := 1 - rel]
  d[, txtcol := fifelse(status == "value" & !is.na(rel) & rel > 0.65, "white", "black")]
  p <- ggplot(d, aes(to, from)) +
    geom_tile(data = d[status == "notmod"], fill = COL_NOTMOD, colour = "white") +
    geom_tile(data = d[status == "noterm"], fill = COL_NOTERM, colour = "white") +
    geom_tile(data = d[status == "bound"], fill = "white", colour = "grey60") +
    geom_tile(data = d[status == "value"], aes(fill = fillv), colour = "white") +
    geom_text(data = d[status %in% c("value", "bound")], aes(label = txt, fontface = ifelse(bold, "bold", "plain"), colour = txtcol), size = 2.1) +
    scale_colour_identity() +
    geom_text(data = d[status == "noterm"], label = "-", size = 2.5, colour = "grey30") +
    scale_x_discrete(labels = lab, drop = FALSE, position = "top") + scale_y_discrete(labels = lab, drop = FALSE) +
    labs(x = "destination", y = "origin", fill = fill_name, title = title, subtitle = subtitle) + coord_equal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0), panel.grid = element_blank())
  if (diverging) p + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = limits,
                                          oob = scales::squish, labels = trans_lab, breaks = breaks)
  else p + scale_fill_viridis_c(option = "magma", direction = seq_dir, limits = limits, oob = scales::squish, labels = trans_lab, breaks = breaks)
}

# ---- inputs -----------------------------------------------------------------------------------------------------------------
pooled <- split_label(fread(file.path(RESULTS_DIR, sprintf("pooled_%s_final.csv", RUN_TAG))))
hr     <- split_label(fread(file.path(RESULTS_DIR, sprintf("hazard_ratios_%s.csv", RUN_TAG)), na.strings = "NA"))
ar     <- split_label(fread(file.path(RESULTS_DIR, sprintf("alpha_report_%s.csv", RUN_TAG))))
iy     <- split_label(fread(file.path(RESULTS_DIR, sprintf("intercept_years_%s.csv", RUN_TAG))))
ts     <- fread(file.path(TABLE_DIR, "transition_set.csv"))
vi     <- readRDS(file.path(TABLE_DIR, "visit_intervals.rds"))
modelled <- unique(pooled$label)
stopifnot(length(modelled) == 90L, all(c(pooled$from, pooled$to) %in% ORDER), nrow(ar) == 90L, nrow(iy) == 90L)
ev <- ts[, .(label = paste0(from, "->", to), from, to, n_events, kept)]
stopifnot(setequal(ev[kept == TRUE, label], modelled))
trow <- merge(ar[, .(label, from, to, class, events)], data.table(label = modelled), by = "label")
trow[, ylab := sprintf("%s (%d)", lab(to), events)]

# scientific-notation helpers (A2, 2026-11-10): cell text "4.9×10⁻³", legend labels 10^k
sup <- function(x) chartr("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079", x)
sci_txt <- function(v) { v <- signif(v, 2); e <- floor(log10(v)); m <- v / 10^e; fifelse(is.na(v), NA_character_, sprintf("%.1f\u00d710%s", m, sup(as.character(e)))) }
pow10_lab <- function(x) parse(text = paste0("10^", x))

gsd <- function(z) ifelse(z == 0, "GDD 0", sprintf("GDD %+d SD", z))   # GDD scenario label; "GDD 0" not "+0 SD" (user 2026-11-10)
