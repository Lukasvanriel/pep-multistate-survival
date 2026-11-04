# ==============================================================================
# R/tests/03t_disturbance_events_test.R -- tests of R/03a_disturbance_events.R (version 2, 2026-10-26 afternoon)
# on invented and real-pattern plots, one or more per rule. No data needed.
# Run from the project root: Rscript R/tests/03t_disturbance_events_test.R
# Expected values are derived BY HAND from the user's rules (docs/PIPELINE_NOTES.md, "2026-10-26 (afternoon)").
# fv_age used here: harvest 4.5, pest 32 (the medians of the 2026-10-26 check); W1 (stand age at the first visit) per test.
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
source(here::here("R", "config.R"))
source(here::here("R", "03a_disturbance_events.R"))

CODES <- list(fire_o = c("BR", "CRB"), fire_p = "BRP",
              harv_o = c("CT", "CPR", "CBA", "CBT", "CEF", "CRS", "CPH", "CPT", "CS", "CDV", "CPP"),
              harv_p = c("CP", "EPC", "CJ", "CDL", "EC", "CEA", "CAM", "CPI", "CPF", "CJG", "CJP", "CJT", "CPS", "CB",
                         "CTR", "DLD", "CPM", "CPC", "CIP", "CA", "CPX", "ESI", "ECE", "CE", "CPR"),
              pest_o = "ES", pest_p = c("EL", "CE"))
FV <- c(harv = 4.5, pest = 32)
N <- NA_character_; I <- NA_integer_
n_ok <- 0L
check <- function(name, year, ba, o, ao, p, ap, expect, W1 = NA_real_, so = NULL, sp = NULL) {
  r <- disturbance_plot(as.integer(year), ba, o, as.integer(ao), p, as.integer(ap), CODES, BA_REGEN, BA_LOSS_THRESH, W1 = W1, fv_age = FV,
                        s_origine = so, s_perturb = sp)$visits
  for (col in names(expect)) {
    got <- r[[col]]; exp <- expect[[col]]
    same <- length(got) == length(exp) && all((is.na(got) & is.na(exp)) | (!is.na(got) & !is.na(exp) & got == exp))
    if (!same) stop(sprintf("FAIL %s, column %s:\n  expected %s\n  got      %s", name, col, paste(exp, collapse = " "), paste(got, collapse = " ")))
  }
  n_ok <<- n_ok + 1L; cat(sprintf("ok  %s\n", name))
}
grow <- function(n) seq(10, by = 5, length.out = n)

# ---- A1: same type, same visit ---------------------------------------------------------------------------------------
check("A1a CT (no year) + CP 1987 at one visit -> heavy 1987",
      c(1980, 1990, 2000), grow(3), c(N, "CT", N), rep(I, 3), c(N, "CP", N), c(I, 1987, I),
      list(harv = c(0, 2, 2), tsh = c(NA, 3, 13)))
check("A1b CT 1985 + CP 1987 at one visit -> heavy 1985 (higher intensity wins, with its year)",
      c(1980, 1990, 2000), grow(3), c(N, "CT", N), c(I, 1985, I), c(N, "CP", N), c(I, 1987, I),
      list(harv = c(0, 2, 2), tsh = c(NA, 5, 15)))
x <- collapse_same_visit(data.frame(v = c(1L, 1L), type = "harv", level = c(1L, 1L), yr = c(1980L, 1985L)), year = 1990L)
stopifnot(nrow(x) == 1L, x$level == 1L, x$yr == 1985L); cat("ok  A1c equal intensity, same interval -> higher year\n"); n_ok <- n_ok + 1L
check("A1d (9109614802) ES 1987 + EL 2020 at one visit, visits in between -> two events; EL 2020 shown in 2025",
      c(1991, 1999, 2015, 2025), c(6.5, 12.2, 21.4, 18.2), c("CT", "ES", "ES", "ES"), c(I, I, 1987, 1987), c(N, N, N, "EL"), c(I, I, I, 2020),
      list(harv = c(2, 2, 2, 2), tsh_imputed_mid = c(10, 18, 34, 44), pest = c(0, 0, 0, 1), sfp = c(NA, NA, NA, 5)), W1 = 10)
check("A1e CP 1980 + CT 1975 at one visit (same interval) -> heavy 1975",
      c(1970, 1990, 2000), grow(3), c(N, "CT", N), c(I, 1975, I), c(N, "CP", N), c(I, 1980, I),
      list(harv = c(0, 2, 2), tsh = c(NA, 15, 25)))
# ---- pair: undated + dated record of one type at one visit ----------------------------------------------------------
check("P1 CT - + CP 1987 at a later visit, no harvest tracked -> heavy 1987",
      c(1980, 1990, 2000), grow(3), c(N, "CT", N), rep(I, 3), c(N, "CP", N), c(I, 1987, I),
      list(harv = c(0, 2, 2), tsh = c(NA, 3, 13)))
check("P2 CT 1970 tracked, then CT - + CP 1987 -> the CT is a leftover, CP 1987 a new light event",
      c(1975, 1990, 2000), grow(3), c("CT", "CT", N), c(1970, I, I), c(N, "CP", N), c(I, 1987, I),
      list(harv = c(2, 1, 1), tsh = c(5, 3, 13)))
check("P3 as P2 with a BA crash at the pair visit -> still light 1987",
      c(1975, 1990, 2000), c(20, 2, 8), c("CT", "CT", N), c(1970, I, I), c(N, "CP", N), c(I, 1987, I),
      list(harv = c(2, 1, 1), tsh = c(5, 3, 13)))
check("P4 only a partial harvest tracked (CP 1970), then CT - + CP 1987 -> heavy 1987",
      c(1975, 1990, 2000), grow(3), c(N, "CT", N), rep(I, 3), c("CP", "CP", N), c(1970, 1987, I),
      list(harv = c(1, 2, 2), tsh = c(5, 3, 13)))
check("P6 (0119600601) CT - + EPC 1998 (2001), CT 1984 + EPC 1998 (2014): CT - back-filled to 1984 first -> heavy 1984",
      c(2001, 2014, 2019), grow(3), c("CT", "CT", N), c(I, 1984, I), c("EPC", "EPC", "EPC"), c(1998, 1998, 1998),
      list(harv = c(2, 2, 2), tsh = c(17, 30, 35)))
check("P5 CT - + CP 1987 at the first visit -> heavy 1987",
      c(1990, 2000), grow(2), c("CT", N), rep(I, 2), c("CP", N), c(1987, I),
      list(harv = c(2, 2), tsh = c(3, 13)), W1 = 40)

# ---- A2: later visit, same type, same interval as a tracked dated event -> one event, highest intensity with its year; equal -> higher year
check("A2a CT 1974 first, CT 1973 later -> 1974",
      c(1995, 2000, 2005, 2010), grow(4), c("CT", "CT", "CT", "CT"), c(1974, 1974, 1973, 1973), rep(N, 4), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(21, 26, 31, 36)))
check("A2b CT 1973 first, CT 1974 later -> 1974",
      c(1995, 2000, 2005, 2010), grow(4), c("CT", "CT", "CT", "CT"), c(1973, 1973, 1974, 1974), rep(N, 4), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(21, 26, 31, 36)))
check("A2c (7200500902) CP 1983 seen 1996, CT 1984 seen 2003, CPR 1984 -> heavy 1984",
      c(1972, 1979, 1996, 2003, 2018), grow(5), c(N, N, N, "CT", "CPR"), c(I, I, I, 1984, 1984), c(N, N, "CP", N, N), c(I, I, 1983, I, I),
      list(harv = c(0, 0, 2, 2, 2), tsh = c(NA, NA, 12, 19, 34)))
check("A2d (7209900402) CP - , CP - , CP 1962, CP 1962, CT 1962 -> heavy 1962 everywhere",
      c(1972, 1979, 1994, 2003, 2015), c(5.4, 5.5, 15.5, 23.5, 28.4), c(N, N, N, N, "CT"), c(I, I, I, I, 1962),
      c("CP", "CP", "CP", "CP", N), c(I, I, 1962, 1962, I),
      list(harv = c(2, 2, 2, 2, 2), tsh = c(10, 17, 32, 41, 53)))

check("A2e (9109609901) BR 1946 seen 2010, BR 1948 seen 2022 -> 1948 (equal intensity, higher year)",
      c(1991, 1998, 2010, 2022), grow(4), c(N, N, "BR", "BR"), c(I, I, 1946, 1948), rep(N, 4), rep(I, 4),
      list(fire = c(1, 1, 1, 1), tsf = c(43, 50, 62, 74)))
check("A2f CT 1985 seen 1990, CP 1987 seen 2000 (same window) -> heavy 1985 (original)",
      c(1980, 1990, 2000), grow(3), c(N, "CT", N), c(I, 1985, I), c(N, N, "CP"), c(I, I, 1987),
      list(harv = c(0, 2, 2), tsh = c(NA, 5, 15)))

# ---- back-fill -----------------------------------------------------------------------------------------------------
check("BF1 (7009400402) CT -, CT - (BA crash), CT 1969 x3 -> 1969 at every visit; dated, so the BA crash at v5 is ignored",
      c(1971, 1978, 1995, 2007, 2023), c(12.4, 1.1, 14.8, 29.8, 0.2), rep("CT", 5), c(I, I, 1969, 1969, 1969), rep(N, 5), rep(I, 5),
      list(harv = c(2, 2, 2, 2, 2), tsh = c(2, 9, 26, 38, 54)))
check("BF2 CT - (1971), CT 1980 (1995), CT 1965 (2007) -> 1971 gets 1965; 1980 is a new event that resets it",
      c(1971, 1983, 1995, 2007), grow(4), c("CT", N, "CT", "CT"), c(I, I, 1980, 1965), rep(N, 4), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(6, 3, 15, 27)))
check("BF3 CP - (1980), CT 1975 later -> back-filled 1975, then A2: heavy 1975",
      c(1980, 1990, 2000), grow(3), c(N, N, "CT"), c(I, I, 1975), c("CP", N, N), rep(I, 3),
      list(harv = c(2, 2, 2), tsh = c(5, 15, 25)))

# ---- B: undated record of a type already tracked ---------------------------------------------------------------------
check("B1 CT 1970, CT - at v3 without BA event -> same event",
      c(1975, 1985, 1997, 2005), grow(4), c("CT", N, "CT", N), c(1970, I, I, I), rep(N, 4), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(5, 15, 27, 35)))
check("B2 CT 1970, CT - at v3 WITH a BA crash -> new event (1985, 1997]",
      c(1975, 1985, 1997, 2005), c(10, 20, 2, 6), c("CT", N, "CT", N), c(1970, I, I, I), rep(N, 4), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(5, 15, NA, NA), tsh_imputed = c(5, 15, 0, 8), tsh_imputed_prev = c(5, 15, 12, 20),
           tsh_imputed_mid = c(5, 15, 6, 14), harv_dating = c("dated", "dated", "bounded", "bounded")))
check("B3 same as B2 but CT 1970 recorded again at v4 -> back-fill dates the v3 record 1970: one event",
      c(1975, 1985, 1997, 2005), c(10, 20, 2, 6), c("CT", N, "CT", "CT"), c(1970, I, I, 1970), rep(N, 4), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(5, 15, 27, 35)))
check("B4 (b/e) EL seen 1990; CT - and EL - at 2000 with a BA crash -> EL continues (drop explained by CT) and CT resets it",
      c(1980, 1990, 2000), c(20, 25, 2), c(N, N, "CT"), rep(I, 3), c(N, "EL", "EL"), rep(I, 3),
      list(pest = c(0, 1, 0), sfp_imputed = c(NA, 0, NA), harv = c(0, 0, 2)))
check("B5 same without the CT: EL - with a BA crash -> new pest event",
      c(1980, 1990, 2000), c(20, 25, 2), rep(N, 3), rep(I, 3), c(N, "EL", "EL"), rep(I, 3),
      list(pest = c(0, 1, 1), sfp_imputed = c(NA, 0, 0)))

# ---- C(i): two events of one type active -------------------------------------------------------------------------------
check("C1 CT 1970, CP 1995 two visits later -> from 1997 the more recent light harvest",
      c(1975, 1985, 1997, 2005, 2015), grow(5), c("CT", N, N, N, N), c(1970, I, I, I, I), c(N, N, "CP", N, N), c(I, I, 1995, I, I),
      list(harv = c(2, 2, 1, 1, 1), tsh = c(5, 15, 2, 10, 20)))

# ---- resets by known order -----------------------------------------------------------------------------------------
check("R1 user example: EL 1970; CT 1980 + EL 1985 recorded 1990 -> 1970 pest removed, 1985 pest stays",
      c(1975, 1983, 1990), grow(3), c(N, N, "CT"), c(I, I, 1980), c("EL", N, "EL"), c(1970, I, 1985),
      list(pest = c(1, 0, 1), sfp = c(5, NA, 5), harv = c(0, 2, 2), tsh = c(NA, 3, 10)))
check("R2 EL 1983 + CT 1984 in one interval (both dated) -> EL known older, removed",
      c(1975, 1990, 2000), grow(3), c(N, N, "CT"), c(I, I, 1984), c(N, "EL", N), c(I, 1983, I),
      list(pest = c(0, 0, 0), harv = c(0, 2, 2), tsh = c(NA, 6, 16)))
check("R3 EL seen 1990, CT - seen 2005 -> EL ended; EL - again 2015 (no BA event) -> new event (old one was ended)",
      c(1980, 1990, 1997, 2005, 2015), c(20, 25, 30, 3, 8), c(N, N, N, "CT", N), rep(I, 5), c(N, "EL", N, N, "EL"), rep(I, 5),
      list(pest = c(0, 1, 1, 0, 1), sfp_imputed = c(NA, 0, 7, NA, 0), harv = c(0, 0, 0, 2, 2)))
check("R4 first-visit EL (undated) + CT 1955 recorded later -> order unknown, both stay",
      c(1971, 1981), grow(2), c(N, "CT"), c(I, 1955), c("EL", N), rep(I, 2),
      list(pest = c(1, 1), sfp_imputed = c(32, 42), pest_dating = rep("first_visit_age", 2), harv = c(2, 2), tsh = c(16, 26)))

# ---- D: old year recorded later --------------------------------------------------------------------------------------
check("D1 EL seen 1989 (v2); CT 1955 recorded 1999 and 2014 -> CT back to 1955, EL not cleared",
      c(1980, 1989, 1999, 2014), grow(4), c(N, N, "CT", "CT"), c(I, I, 1955, 1955), c(N, "EL", N, N), rep(I, 4),
      list(harv = c(2, 2, 2, 2), tsh = c(25, 34, 44, 59), pest = c(0, 1, 1, 1),
           sfp_imputed = c(NA, 0, 10, 25), sfp_imputed_prev = c(NA, 9, 19, 34)))
check("D2 (7009501001) BR -, BR -, CT -, -, BR 1947 -> fire 1947 (back-fill), ended by the 1992 clear-cut",
      c(1971, 1978, 1992, 2003, 2016), grow(5), c("BR", "BR", "CT", N, "BR"), c(I, I, I, I, 1947), rep(N, 5), rep(I, 5),
      list(fire = c(1, 1, 0, 0, 0), tsf = c(24, 31, NA, NA, NA), harv = c(0, 0, 2, 2, 2),
           tsh_imputed = c(NA, NA, 0, 11, 24), tsh_imputed_prev = c(NA, NA, 14, 25, 38)))
check("D3 CT 1990 recorded 2005 -> filled back only to visits after 1990",
      c(1985, 1995, 2005), grow(3), c(N, N, "CT"), c(I, I, 1990), rep(N, 3), rep(I, 3),
      list(harv = c(0, 2, 2), tsh = c(NA, 5, 15)))

# ---- selection: stand-origin events of different types, unknown order or same year ------------------------------------
check("S1 (8809416701) CT 1995 seen 1998, BR 1995 seen 2014 -> same year: harvest kept",
      c(1988, 1998, 2014, 2023), grow(4), c(N, "CT", "BR", "BR"), c(I, 1995, 1995, 1995), rep(N, 4), rep(I, 4),
      list(fire = c(0, 0, 0, 0), harv = c(0, 2, 2, 2), tsh = c(NA, 3, 19, 28)))
check("S2 (7507600602) BR - (first visit) + CT 1950 -> order unknown, harvest kept",
      c(1975, 1980, 1992), grow(3), c("BR", "BR", "CT"), c(I, I, 1950), rep(N, 3), rep(I, 3),
      list(fire = c(0, 0, 0), harv = c(2, 2, 2), tsh = c(25, 30, 42)), W1 = 40)
check("S3 ES - (first visit) + BR 1950 -> order unknown, no harvest: the dated one (fire) kept",
      c(1975, 1985), grow(2), c("ES", "BR"), c(I, 1950), rep(N, 2), rep(I, 2),
      list(fire = c(1, 1), tsf = c(25, 35), pest = c(0, 0)), W1 = 40)
check("S4 BR 1983 + CT 1984 in one interval (both dated, known order) -> the later (CT) resets the fire",
      c(1975, 1990, 2000), grow(3), c(N, "BR", "CT"), c(I, 1983, 1984), rep(N, 3), rep(I, 3),
      list(fire = c(0, 0, 0), harv = c(0, 2, 2), tsh = c(NA, 6, 16)))
stopifnot(choose_stand_origin(1:2, c("fire", "pest"), c(NA, NA)) == 1L, choose_stand_origin(1:2, c("fire", "pest"), c(NA, 1950L)) == 2L,
          choose_stand_origin(1:2, c("fire", "harv"), c(1990L, NA)) == 2L, choose_stand_origin(1:2, c("fire", "harv"), c(1995L, 1995L)) == 2L,
          choose_stand_origin(1:2, c("pest", "fire"), c(1995L, 1995L)) == 2L)
cat("ok  S5 choose_stand_origin (5 cases)\n"); n_ok <- n_ok + 1L

# ---- CE, same-year rule, year after the visit ---------------------------------------------------------------------------
check("CE1 CE - at v2 -> partial harvest AND partial pest, bounded",
      c(1980, 1990, 2000), grow(3), rep(N, 3), rep(I, 3), c(N, "CE", N), rep(I, 3),
      list(harv = c(0, 1, 1), pest = c(0, 1, 1), tsh_imputed = c(NA, 0, 10), sfp_imputed = c(NA, 0, 10), tsh_imputed_prev = c(NA, 10, 20)))
check("Y1 CT 1990 recorded at the 1990 visit, no BA drop -> active from the next visit",
      c(1980, 1990, 2000), grow(3), c(N, "CT", N), c(I, 1990, I), rep(N, 3), rep(I, 3),
      list(harv = c(0, 0, 2), tsh = c(NA, NA, 10)))
check("Y2 same with a BA drop -> active at the 1990 visit",
      c(1980, 1990, 2000), c(20, 2, 6), c(N, "CT", N), c(I, 1990, I), rep(N, 3), rep(I, 3),
      list(harv = c(0, 2, 2), tsh = c(NA, 0, 10)))
check("H1 (9109609601 part) CEA 2004 recorded at 2003 (after the visit), again at 2017 -> active from 2017",
      c(1991, 1997, 2003, 2017), grow(4), rep(N, 4), rep(I, 4), c(N, N, "CEA", "CEA"), c(I, I, 2004, 2004),
      list(harv = c(0, 0, 0, 1), tsh = c(NA, NA, NA, 13)))

# ---- first visit, undated --------------------------------------------------------------------------------------------
check("F1 CT - at the first visit -> TSD = stand age W1 (+ elapsed)",
      c(1975, 1985), grow(2), c("CT", N), rep(I, 2), rep(N, 2), rep(I, 2),
      list(harv = c(2, 2), tsh = c(NA, NA), tsh_imputed = c(40, 50), tsh_imputed_mid = c(40, 50), tsh_imputed_prev = c(40, 50),
           harv_dating = rep("first_visit_standage", 2)), W1 = 40)
check("F2 CT - at the first visit, CT 1960 later -> dated 1960 (back-fill)",
      c(1975, 1985, 1995), grow(3), c("CT", N, "CT"), c(I, I, 1960), rep(N, 3), rep(I, 3),
      list(harv = c(2, 2, 2), tsh = c(15, 25, 35), harv_dating = rep("dated", 3)), W1 = 40)
check("F3 EL - at the first visit, repeated -> flagged, age 32 at the first visit",
      c(1975, 1985, 1995), grow(3), rep(N, 3), rep(I, 3), c("EL", "EL", N), rep(I, 3),
      list(pest = c(1, 1, 1), sfp_imputed = c(32, 42, 52), pest_dating = rep("first_visit_age", 3)))
check("F4 BRP - at the first visit (partial fire) -> stand age",
      c(1975, 1985), grow(2), rep(N, 2), rep(I, 2), c("BRP", N), rep(I, 2),
      list(fire = c(1, 1), tsf_imputed = c(40, 50), fire_dating = rep("first_visit_standage", 2)), W1 = 40)
check("F5 CP - at the first visit -> age 4.5 at the first visit",
      c(1975, 1985), grow(2), rep(N, 2), rep(I, 2), c("CP", N), rep(I, 2),
      list(harv = c(1, 1), tsh_imputed = c(4.5, 14.5), harv_dating = rep("first_visit_age", 2)))

check("BA1 (9003404902 pattern) fire tracked, BA 1.1 -> 0 with an undated BR: previous BA < 5 -> no BA event -> same fire",
      c(1990, 1999, 2014), c(1.4, 1.1, 0), c("BR", "BR", "BR"), rep(I, 3), rep(N, 3), rep(I, 3),
      list(fire = c(1, 1, 1), tsf_imputed_mid = c(30, 39, 54)), W1 = 30)
# ---- station merge (step 0b): station codes undated; same type -> pee wins; both stand-origin -> pee wins; else both -------
check("ST1 pee CP + station CT at one visit (same type) -> pee wins: light",
      c(1980, 1990, 2000), grow(3), rep(N, 3), rep(I, 3), c(N, "CP", N), c(I, 1987, I),
      list(harv = c(0, 1, 1), tsh = c(NA, 3, 13)), so = c(N, "CT", N), sp = rep(N, 3))
check("ST2 pee BR 1987 + station CT at one visit (both stand-origin) -> pee wins: fire only",
      c(1980, 1990, 2000), grow(3), c(N, "BR", N), c(I, 1987, I), rep(N, 3), rep(I, 3),
      list(fire = c(0, 1, 1), harv = c(0, 0, 0), tsf = c(NA, 3, 13)), so = c(N, "CT", N), sp = rep(N, 3))
check("ST3 pee EL 1987 + station CT (different types) -> both: pest 1987 and heavy harvest dated in the interval (mid 1985)",
      c(1980, 1990, 2000), grow(3), rep(N, 3), rep(I, 3), c(N, "EL", N), c(I, 1987, I),
      list(pest = c(0, 1, 1), sfp = c(NA, 3, 13), harv = c(0, 2, 2), tsh_imputed_mid = c(NA, 5, 15)), so = c(N, "CT", N), sp = rep(N, 3))
check("ST4 pee CT 1987 + station EL (different types) -> both; the CT does not reset the undated EL of the same interval",
      c(1980, 1990, 2000), grow(3), c(N, "CT", N), c(I, 1987, I), rep(N, 3), rep(I, 3),
      list(harv = c(0, 2, 2), tsh = c(NA, 3, 13), pest = c(0, 1, 1), sfp_imputed_mid = c(NA, 5, 15)), so = rep(N, 3), sp = c(N, "EL", N))
check("ST5 (7000700401) station CT at the BA-crash visit 1985, pee CT - in 1994 -> heavy from 1985 (mid 1981.5), continued in 1994",
      c(1970, 1978, 1985, 1994, 2009, 2022), c(25.4, 32.8, 0.3, 0, 11.2, 23.5), c(N, N, N, "CT", N, N), rep(I, 6), c("CP", "CP", N, N, N, N), rep(I, 6),
      list(harv = c(1, 1, 2, 2, 2, 2), tsh_imputed_mid = c(4.5, 12.5, 3.5, 12.5, 27.5, 40.5)), so = c(N, N, "CT", "CT", N, N), sp = rep(N, 6))
check("ST6 station-only CP at the 2nd visit -> light harvest, mid of the interval",
      c(1980, 1990, 2000), grow(3), rep(N, 3), rep(I, 3), rep(N, 3), rep(I, 3),
      list(harv = c(0, 1, 1), tsh_imputed_mid = c(NA, 5, 15)), so = rep(N, 3), sp = c(N, "CP", N))
check("Z1 plot without any disturbance code -> nothing flagged",
      c(1980, 1990, 2000), c(20, 2, 6), rep(N, 3), rep(I, 3), rep(N, 3), rep(I, 3),
      list(fire = c(0, 0, 0), harv = c(0, 0, 0), pest = c(0, 0, 0), tsh_imputed_mid = c(NA, NA, NA)))
stopifnot(identical(ba_event_flag(c(3, 10, 7.4, 7.6, 4, 0), 5, 0.25), c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)))   # 4 -> 0: previous BA already < 5 -> no BA event (2026-10-26 evening)
cat("ok  E0 ba_event_flag\n"); n_ok <- n_ok + 1L
cat(sprintf("\nAll %d tests passed.\n", n_ok))
