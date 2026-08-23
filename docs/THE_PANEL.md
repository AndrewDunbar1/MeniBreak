# The breakpoint panel

The panel is **not a fixed hard-coded checklist**. It is a procedure: within each
training fold, the 20 arm-relative bins per autosome arm are ranked by
chi-square association with recurrence and the top *k* (=8) are retained. The
windows below are the **consensus** — those the procedure re-selects across folds
— reported with their selection frequency. A window chosen in 90%+ of folds is a
finding; one at 40% is close to coincidence. (This is why the earlier
"fixed nine-window" description was wrong: the nine windows are the *union* of
what was selected in ≥40% of folds at one or more horizons, not the inputs to
every fold.)

## Consensus windows (GRCh37; arm-relative 0 = telomere, 1 = centromere)

| Arm | Arm-relative | GRCh37 | Sel. 3 yr | Sel. 5 yr | Sel. 7 yr |
|---|---|---|---:|---:|---:|
| 1p | 0.90–0.95 | chr1:112.50–118.75 Mb | **1.00** | — | 0.44 |
| **1p** | **0.95–1.00** | **chr1:118.75–125.00 Mb** | **0.96** | **0.88** | **0.92** |
| 6q | 0.90–0.95 | chr6:66.51–72.01 Mb | **1.00** | 0.64 | — |
| 6q | 0.95–1.00 | chr6:61.00–66.51 Mb | 0.40 | — | 0.68 |
| 19q | 0.50–0.55 | chr19:41.18–42.81 Mb | — | 0.80 | 0.68 |
| 19q | 0.45–0.50 | chr19:42.81–44.45 Mb | — | 0.72 | — |
| 2p | 0.75–0.80 | chr2:69.97–74.64 Mb | — | 0.76 | — |
| 4p | 0.95–1.00 | chr4:47.88–50.40 Mb | 0.44 | — | — |
| 3p | 0.50–0.55 | chr3:45.50–50.05 Mb | — | — | 0.40 |

Source: `results/reference/hotlist_windows.csv` (from `analyses/02_discover_hotlist.py`).

**How to read this table.** The selection frequencies above are the **k = 5**
rows of `hotlist_windows.csv`, and “—” means the window was selected in fewer
than 40% of folds at that horizon. The discovery script sweeps several panel
sizes, so the CSV contains **k ∈ {3, 5, 8, 15}**; filter on `k` before quoting a
number. Note that the *performance* ladder in `results/reference/full_ladder.csv`
("8-window panel") is fitted at **k = 8**, so the stability table and the AUC
table are reported at different k. This is stated explicitly here because the
original report showed the two side by side without saying so, which made the
panel look like a fixed nine-window list. The headline window is more stable, not
less, at the k used for performance: at k = 8 the chr1p 0.95–1.00 window is
selected in **100% / 96% / 96%** of folds at 3 / 5 / 7 years.

## The headline window

The **chr1p pericentromeric window, chr1:118.75–125.00 Mb**, is the single most
stable feature: selected in **96% / 88% / 92%** of folds at 3 / 5 / 7 years
(k = 5; 100% / 96% / 96% at k = 8). It is
recovered **without** the frozen list being an input to selection — an internal
positive control that the procedure rediscovers the known chr1p hotspot rather
than being told about it.

## Positive control

`reference/hotspot_windows_v1.csv` holds the locked model's frozen arm-relative
hotspot windows. These are used **only** as a control to score a rediscovered
window against — never as a training input or label. On the synthetic tutorial
cohort, `apply/menibreak_panel.py` likewise rediscovers a planted chr1p window
(`bp_1p_b18`, chr1:112.5–118.75 Mb) at selection frequency 1.0, which is how you
can confirm the discovery machinery works before trusting it on real data.

## Why breakpoints, not just SCNA

Arm-level SCNA records *that* an arm is lost; the breakpoint records *where the
loss ends*. Two tumours with the same 1p loss burden but different change-point
positions get the same SCNA call and different panel features. The panel is the
only marker whose discrimination **improves** from 3 to 7 years (+0.050 broad,
+0.042 restricted), while WHO grade, integrated grade and SCNA all decline.
