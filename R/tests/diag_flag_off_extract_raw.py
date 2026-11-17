# Read-only export of three raw PEP layers to CSV (input of R/tests/diag_flag_off_interval1.R).
# Run from the project root:  python3 R/tests/diag_flag_off_extract_raw.py
import sqlite3, csv, os
os.makedirs("data-output/diagnostics/flag_off_raw", exist_ok=True)
con = sqlite3.connect("file:data/raw/pep/PEP.gpkg?mode=ro", uri=True)
for t in ["pee_ori_sond", "station_pe", "placette_mes"]:
    cur = con.execute(f"select * from {t}")
    keep = [i for i, d in enumerate(cur.description) if d[0] not in ("geom", "geometry")]
    with open(f"data-output/diagnostics/flag_off_raw/{t}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow([cur.description[i][0] for i in keep])
        for r in cur: w.writerow([r[i] for i in keep])
