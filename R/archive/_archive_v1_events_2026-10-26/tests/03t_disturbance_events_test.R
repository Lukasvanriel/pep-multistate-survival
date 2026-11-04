# ==============================================================================
# R/tests/03t_disturbance_events_test.R -- tests of the disturbance event logic (R/03a_disturbance_events.R)
# on invented plots, one per rule. No data needed. Run from the project root: Rscript R/tests/03t_disturbance_events_test.R
# Each case lists the exact expected flag and TSD per visit; the script stops at the first mismatch.
# Expected values follow the rules agreed 2026-10-26 (docs/PIPELINE_NOTES.md, "Disturbance event rewrite").
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
source(here::here("R", "config.R"))
source(here::here("R", "03a_disturbance_events.R"))

CODES <- list(fire_o = c("BR", "CRB"), fire_p = "BRP",
              harv_o = c("CT", "CPR", "CBA", "CBT", "CEF", "CRS", "CPH", "CPT", "CS", "CDV", "CPP"),
              harv_p = c("CP", "EPC", "CJ", "CDL", "EC", "CEA", "CAM", "CPI", "CPF", "CJG", "CJP", "CJT", "CPS", "CB",
                         "CTR", "DLD", "CPM", "CPC", "CIP", "CA", "CPX", "ESI", "ECE", "CE", "CPR"),
              pest_o = "ES", pest_p = c("EL", "CE"))
N <- NA_character_; I <- NA_integer_
n_ok <- 0L

check <- function(name, year, ba, o, ao, p, ap, expect) {
  r <- disturbance_plot(as.integer(year), ba, o, as.integer(ao), p, as.integer(ap), CODES, BA_REGEN, BA_LOSS_THRESH)$visits
  for (col in names(expect)) {
    got <- r[[col]]; exp <- expect[[col]]
    same <- length(got) == length(exp) && all((is.na(got) & is.na(exp)) | (!is.na(got) & !is.na(exp) & got == exp))
    if (!same) stop(sprintf("FAIL %s, column %s:\n  expected %s\n  got      %s", name, col,
                            paste(exp, collapse = " "), paste(got, collapse = " ")))
  }
  n_ok <<- n_ok + 1L
  cat(sprintf("ok  %s\n", name))
}

grow <- function(n) seq(10, by = 5, length.out = n)   # basal area that only grows: no ba_event after visit 1

# ---- A: same window, same type -----------------------------------------------------------------------------------
check("A1 same window, later record later year: CT 1973 (v1-2), CT 1974 (v3-4) -> one event, 1974",
      c(1995, 2000, 2005, 2010), grow(4), c("CT", "CT", "CT", "CT"), c(1973, 1973, 1974, 1974), rep(N, 4), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(21, 26, 31, 36)))
check("A2 same window, later record earlier year: CT 1974 (v1-2), CT 1973 (v3-4) -> 1974 kept",
      c(1995, 2000, 2005, 2010), grow(4), c("CT", "CT", "CT", "CT"), c(1974, 1974, 1973, 1973), rep(N, 4), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(21, 26, 31, 36)))
check("A3 same window, levels differ: CP 1983 seen 1996, CT 1984 seen 2003 -> heavy, 1984 (latest year sets level)",
      c(1972, 1979, 1996, 2003, 2018), grow(5), c(N, N, N, "CT", "CPR"), c(I, I, I, 1984, 1984), c(N, N, "CP", N, N), c(I, I, 1983, I, I),
      list(harv = c(0, 0, 2, 2, 2), tsh = c(NA, NA, 12, 19, 34)))
check("A4 same visit, same year, both levels: CT 1990 + CP 1990 -> heavy (tie: more severe)",
      c(1985, 1995, 2005), grow(3), c(N, "CT", N), c(I, 1990, I), c(N, "CP", N), c(I, 1990, I),
      list(harv = c(0, 2, 2), tsh = c(NA, 5, 15)))
check("A5 two stand-origin types same window: BR 1983 + CT 1984 (salvage) -> only the later (harvest) kept",
      c(1975, 1990, 2000), grow(3), c(N, "BR", "CT"), c(I, 1983, 1984), rep(N, 3), rep(I, 3),
      list(fire = c(0, 0, 0), harv = c(0, 2, 2), tsh = c(NA, 6, 16)))
check("A6 partial + stand-origin of different types in one window: EL 1983 + CT 1984 -> both kept",
      c(1975, 1990, 2000), grow(3), c(N, N, "CT"), c(I, I, 1984), c(N, "EL", N), c(I, 1983, I),
      list(pest = c(0, 1, 1), sfp = c(NA, 7, 17), harv = c(0, 2, 2), tsh = c(NA, 6, 16)))
check("S1 same year: CT 1995 (seen 1998) + BR 1995 (seen 2014) (8809416701) -> harvest kept",
      c(1988, 1998, 2014, 2023), grow(4), c(N, "CT", "BR", "BR"), c(I, 1995, 1995, 1995), rep(N, 4), rep(I, 4),
      list(fire = c(0, 0, 0, 0), harv = c(0, 2, 2, 2), tsh = c(NA, 3, 19, 28)))
check("S2 order unknown: undated BR + CT 1950, both before the first visit (7507600602) -> harvest kept",
      c(1975, 1980, 1992), grow(3), c("BR", "BR", "CT"), c(I, I, 1950), rep(N, 3), rep(I, 3),
      list(fire = c(0, 0, 0), harv = c(2, 2, 2), tsh = c(25, 30, 42)))
check("S3 fire + pest, one dated: undated ES + BR 1950 before the first visit -> the dated one (fire) kept",
      c(1975, 1985), grow(2), c("ES", "BR"), c(I, 1950), rep(N, 2), rep(I, 2),
      list(fire = c(1, 1), tsf = c(25, 35), pest = c(0, 0)))
# S4: two undated stand-origin events cannot share a start visit in pee data (one origine code per visit), so the
# selection rule is tested directly: (type, year) -> index kept
stopifnot(choose_stand_origin(1:2, c("fire", "pest"), c(NA, NA)) == 1L,          # both undated fire + pest -> fire
          choose_stand_origin(1:2, c("pest", "fire"), c(NA, NA)) == 2L,
          choose_stand_origin(1:2, c("fire", "pest"), c(NA, 1950L)) == 2L,       # one dated -> the dated one
          choose_stand_origin(1:2, c("fire", "harv"), c(1990L, NA)) == 2L,       # any undated, harvest present -> harvest
          choose_stand_origin(1:2, c("fire", "harv"), c(1995L, 1995L)) == 2L,    # same year -> harvest
          choose_stand_origin(1:2, c("pest", "fire"), c(1995L, 1995L)) == 2L,    # same year, no harvest -> fire
          choose_stand_origin(1:3, c("harv", "fire", "pest"), c(1980L, 1985L, 1983L)) == 2L)  # all dated -> latest year
cat("ok  S4 choose_stand_origin (7 cases)\n"); n_ok <- n_ok + 1L

# ---- 1: dated vs undated, same type --------------------------------------------------------------------------------
check("U1 undated CP seen 1996, later CT dated 1990 (<= 1996) -> same event, dated record sets year and level",
      c(1985, 1996, 2005), grow(3), c(N, N, "CT"), c(I, I, 1990), c(N, "CP", N), c(I, I, I),
      list(harv = c(0, 2, 2), tsh = c(NA, 6, 15)))
check("U2 undated CP seen 1996, later CT dated 2000 (> 1996) -> new event, the most recent sets the values",
      c(1985, 1996, 2005), grow(3), c(N, N, "CT"), c(I, I, 2000), c(N, "CP", N), c(I, I, I),
      list(harv = c(0, 1, 2), tsh = c(NA, NA, 5), tsh_imputed = c(NA, 0, 5), tsh_imputed_prev = c(NA, 11, 5)))

# ---- B: undated record of a type already tracked ---------------------------------------------------------------------
check("B1 CT 1970, CT again undated at v3 with no BA loss -> same event (old year, old level)",
      c(1975, 1985, 1997, 2005), grow(4), c("CT", N, "CT", N), c(1970, I, I, I), rep(N, 4), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(5, 15, 27, 35)))
check("B2 CT 1970, CT again undated at v3 WITH a BA crash -> new event, bounded (1985, 1997]",
      c(1975, 1985, 1997, 2005), c(10, 20, 2, 6), c("CT", N, "CT", N), c(1970, I, I, I), rep(N, 4), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(5, 15, NA, NA), tsh_imputed = c(5, 15, 0, 8), tsh_imputed_prev = c(5, 15, 12, 20),
           tsh_imputed_mid = c(5, 15, 6, 14), harv_dating = c("dated", "dated", "bounded", "bounded")))
check("B3 CT 1970, undated CP at v3 without BA loss -> same event (stays heavy)",
      c(1975, 1985, 1997, 2005), grow(4), c("CT", N, N, N), c(1970, I, I, I), c(N, N, "CP", N), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(5, 15, 27, 35)))
check("B4 CT 1970, undated CP at v3 with BA loss > 25% -> new partial event; most recent -> light",
      c(1975, 1985, 1997, 2005), c(10, 20, 14, 18), c("CT", N, N, N), c(1970, I, I, I), c(N, N, "CP", N), rep(I, 4),
      list(harv = c(2, 2, 1, 1), tsh = c(5, 15, NA, NA), tsh_imputed = c(5, 15, 0, 8)))
check("B5 a code simply repeated with a BA crash also counts (rule applies to repeated codes too)",
      c(1975, 1985, 1997), c(10, 20, 3), c("CT", "CT", "CT"), c(I, I, I), rep(N, 3), rep(I, 3),
      list(harv = c(2, 2, 2), tsh_imputed = c(NA, NA, 0), harv_dating = c("first_visit", "first_visit", "bounded")))

# ---- C(i): heavy then light, different windows ------------------------------------------------------------------------
check("C1 CT 1970, CP 1995 two visits later -> both events; from 1997 the light harvest sets level and TSD",
      c(1975, 1985, 1997, 2005, 2015), grow(5), c("CT", N, N, N, N), c(1970, I, I, I, I), c(N, N, "CP", N, N), c(I, I, 1995, I, I),
      list(harv = c(2, 2, 1, 1, 1), tsh = c(5, 15, 2, 10, 20)))

# ---- D: stand-origin code dated before the previous visit -------------------------------------------------------------
check("D1 EL seen 1989 (not first visit), CT 1955 recorded 1999 and 2014 -> CT back-filled, EL NOT cleared",
      c(1980, 1989, 1999, 2014), grow(4), c(N, N, "CT", "CT"), c(I, I, 1955, 1955), c(N, "EL", N, N), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(25, 34, 44, 59), pest = c(0, 1, 1, 1),
           sfp_imputed = c(NA, 0, 10, 25), sfp_imputed_prev = c(NA, 9, 19, 34)))
check("D2 (7009501001) BR undated 1971/1978, CT undated 1992, BR 1947 recorded 2016 -> fire 1947 back-filled; CT 1992 ends it",
      c(1971, 1978, 1992, 2003, 2016), grow(5), c("BR", "BR", "CT", N, "BR"), c(I, I, I, I, 1947), rep(N, 5), rep(I, 5),
      list(fire = c(1, 1, 0, 0, 0), tsf = c(24, 31, NA, NA, NA), harv = c(0, 0, 2, 2, 2),
           tsh_imputed = c(NA, NA, 0, 11, 24), tsh_imputed_prev = c(NA, NA, 14, 25, 38)))
check("D3 back-fill only as far as chronologically possible: CT 1990 recorded 2005 -> not active at 1985",
      c(1985, 1995, 2005), grow(3), c(N, N, "CT"), c(I, I, 1990), rep(N, 3), rep(I, 3),
      list(harv = c(0, 2, 2), tsh = c(NA, 5, 15)))

# ---- CE: two records -------------------------------------------------------------------------------------------------
check("CE1 undated CE at v2 -> partial harvest AND partial pest, both bounded",
      c(1980, 1990, 2000), grow(3), rep(N, 3), rep(I, 3), c(N, "CE", N), rep(I, 3),
      list(harv = c(0, 1, 1), pest = c(0, 1, 1), tsh_imputed = c(NA, 0, 10), sfp_imputed = c(NA, 0, 10),
           tsh_imputed_prev = c(NA, 10, 20)))

# ---- same-year rule ------------------------------------------------------------------------------------------------
check("Y1 CT dated 1990 recorded at the 1990 visit, no BA drop -> after the measurement: active from the next visit",
      c(1980, 1990, 2000), grow(3), c(N, "CT", N), c(I, 1990, I), rep(N, 3), rep(I, 3),
      list(harv = c(0, 0, 2), tsh = c(NA, NA, 10)))
check("Y2 same, with a BA drop -> before the measurement: active at the 1990 visit",
      c(1980, 1990, 2000), c(20, 2, 6), c(N, "CT", N), c(I, 1990, I), rep(N, 3), rep(I, 3),
      list(harv = c(0, 2, 2), tsh = c(NA, 0, 10)))

# ---- first visit -------------------------------------------------------------------------------------------------------
check("F1 undated CT at the first visit -> active, no TSD, marker first_visit (04 stand-age rule)",
      c(1975, 1985), grow(2), c("CT", N), c(I, I), rep(N, 2), rep(I, 2),
      list(harv = c(2, 2), tsh = c(NA, NA), tsh_imputed = c(NA, NA), harv_dating = c("first_visit", "first_visit")))
check("F2 undated CT at the first visit, later CT dated 1960 (<= 1975) -> dated",
      c(1975, 1985, 1995), grow(3), c("CT", N, "CT"), c(I, I, 1960), rep(N, 3), rep(I, 3),
      list(harv = c(2, 2, 2), tsh = c(15, 25, 35), harv_dating = c("dated", "dated", "dated")))
check("F3 undated EL at the first visit -> never flagged (T3d), also when repeated later without BA loss",
      c(1975, 1985, 1995), grow(3), rep(N, 3), rep(I, 3), c("EL", "EL", N), rep(I, 3),
      list(pest = c(0, 0, 0)))
check("F4 undated ES at the first visit -> active (stand-origin), marker first_visit",
      c(1971, 1978, 1992), grow(3), c("ES", "ES", N), rep(I, 3), rep(N, 3), rep(I, 3),
      list(pest = c(1, 1, 1), sfp_imputed = c(NA, NA, NA), pest_dating = rep("first_visit", 3)))

# ---- resets ----------------------------------------------------------------------------------------------------------
check("R1 EL seen 1990, CT undated seen 2005 -> EL ended at 2005; EL recorded again 2015 (no BA loss) -> new event",
      c(1980, 1990, 1997, 2005, 2015), c(20, 25, 30, 3, 8), c(N, N, N, "CT", N), rep(I, 5), c(N, "EL", N, N, "EL"), rep(I, 5),
      list(pest = c(0, 1, 1, 0, 1), sfp_imputed = c(NA, 0, 7, NA, 0), harv = c(0, 0, 0, 2, 2)))
check("R2 same visit: new CT and an EL already tracked -> EL ended by CT; EL record at that visit = new pest event (tie, not ended)",
      c(1980, 1990, 2000, 2010), c(20, 25, 2, 6), c(N, N, "CT", N), rep(I, 4), c(N, "EL", "EL", N), rep(I, 4),
      list(pest = c(0, 1, 1, 1), sfp_imputed = c(NA, 0, 0, 10), harv = c(0, 0, 2, 2)))

check("R3 two stand-origin types BEFORE the first visit, both dated: BR 1947 + CT 1955 -> the later (CT) kept",
      c(1971, 1981), grow(2), c("BR", "CT"), c(1947, 1955), rep(N, 2), rep(I, 2),
      list(fire = c(0, 0), harv = c(2, 2), tsh = c(16, 26)))

# ---- E0 ------------------------------------------------------------------------------------------------------------------
stopifnot(identical(ba_event_flag(c(3, 10, 7.4, 7.6, 4, 0), 5, 0.25), c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)))
cat("ok  E0 ba_event_flag\n"); n_ok <- n_ok + 1L

cat(sprintf("\nAll %d tests passed.\n", n_ok))
