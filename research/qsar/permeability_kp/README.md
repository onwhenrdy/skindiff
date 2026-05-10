# QSAR models for skin permeability coefficient (Kp)

Literature notes on QSARs whose regression *output* is the skin permeability
coefficient Kp (cm/s or cm/h). Sister directories cover models whose output
is D, K, or Jmax separately.

Sorted by year. R^2 is on the training set unless otherwise noted.

| Model              | Year | Descriptors                                  | Training n      | R^2          | Domain                                                | DOI / Source                                            |
|--------------------|-----:|----------------------------------------------|----------------:|-------------:|-------------------------------------------------------|---------------------------------------------------------|
| Potts-Guy          | 1992 | logKow, MW                                   | 93 (Flynn)      | 0.67         | aqueous donor; MW 18-750; logKow -3 to +6             | 10.1023/A:1015810312465                                 |
| Cleek-Bunge        | 1993 | adds aqueous-epidermis correction to any Kp_sc | n/a (correction) | n/a        | lipophilic solutes (logKow > 3) where viable epidermis limits | 10.1023/A:1018981515480                            |
| Abraham (skin LFER)| 1995 | pi2H, alpha2H, beta2H, Vx                    | ~37-50 (Flynn)  | ~0.83-0.86   | neutral organics with Abraham descriptors             | 10.1111/j.2042-7158.1995.tb05725.x                       |
| Wilschut / mod-Robinson | 1995 | logKow, MW (parallel-pathway)            | 123 (99 cmpds)  | best of 5    | aqueous donor; SKINPERM precursor                     | 10.1016/0045-6535(95)00023-2                             |
| Lien-Gao           | 1995 | logKow, (logKow)^2, log MW, Hb               | 53              | ~0.92 R       | human + mouse skin; bilinear in Kow                   | 10.1023/A:1016266316100                                 |
| Barratt            | 1995 | logKow, MV, MPt                              | ~60-65 (Flynn)  | ~0.90        | requires melting point                                | 10.1016/0887-2333(94)00190-6                             |
| Potts-Guy (LFER)   | 1995 | MV, alpha2H, beta2H                          | ~37 (Flynn)     | ~0.94        | neutral organics; H-bond refinement                   | 10.1023/A:1016236932339                                 |
| Abraham (revised)  | 1999 | R2, pi2H, alpha2H, beta2H, Vx                | ~53 -> 119      | ~0.83-0.87   | neutral organics; canonical skin LFER                 | 10.1002/(SICI)1096-9063(199901)55:1<78::AID-PS859>3.0.CO;2-Y |
| Mitragotri (SPT)   | 2002 | Kow, solute radius r                         | ~30-50 fit-anchor | ~0.70 (val) | small hydrophobic solutes; lipid-pathway only         | 10.1002/jps.10071                                       |
| Frasch (random walk) | 2002 | MW, logKow (mechanistic)                   | 94 (Flynn)      | 0.84         | brick-and-mortar SC simulator                         | 10.1111/0272-4332.00024                                 |
| Moss-Cronin        | 2002 | logKow, MW                                   | 116             | 0.82         | Flynn re-fit with corrected steroids                  | 10.1016/S0378-5173(02)00057-1                            |
| Patel-ten Berge-Cronin | 2002 | logKow, MW, ABSQon, SsssCH               | 158 (143 used)  | 0.90         | curated human-skin set; needs E-state + charges       | 10.1016/S0045-6535(02)00114-5                            |
| Magnusson (Jmax)   | 2004 | MW (and Mpt, Ha refinements)                 | 87 - 278        | 0.69 - 0.92  | predicts Jmax not Kp; saturated aqueous donor         | 10.1111/j.0022-202X.2004.22413.x                         |
| ten Berge / SKINPERM | 2009 | logKow, MW (two-pathway split)             | 182             | ~0.78-0.85   | occupational risk; MW 18-584; logKow -3.7 to +5.5     | 10.1080/15459624.2013.831983 (companion)                 |
| Baba (SVR)         | 2015 | many (ML, RBF kernel)                        | 211             | 0.91 (test)  | curated single-protocol aqueous Kp                    | 10.1007/s11095-015-1629-y                                |
| Chen (MLR + SVM)   | 2018 | ALOGP, X3v, Neoplastic-80                    | 274 (139 train) | 0.90 / 0.87 (test, SVM) | extends Patel set; bilinear in ALOGP        | 10.3390/molecules23040911                                |

## How they relate

- **Potts-Guy 1992** is the foundational QSAR. Almost every model below
  benchmarks against it.
- **Cleek-Bunge 1993** is a series-resistance correction layered on top of
  Potts-Guy (or any Kp_sc QSAR). Add it whenever logKow > ~3.
- **Wilschut 1995** is the methodological backbone of **ten Berge 2009 /
  SKINPERM** (the workbook used in EU/US occupational risk assessment).
- **Lien-Gao 1995, Patel 2002, Moss-Cronin 2002, Chen 2018** are
  Potts-Guy-style empirical re-analyses on progressively larger / cleaner
  datasets. Moss-Cronin is the easiest drop-in replacement for Potts-Guy if
  the only inputs are MW and logKow.
- **Abraham 1995 / 1999 (and Abraham-Martins 2004 follow-up)** is the LFER
  family: more accurate but requires five solvation descriptors per compound.
- **Potts-Guy 1995** sits between the empirical and LFER worlds (uses MV +
  H-bond descriptors but no Kow).
- **Barratt 1995** uses melting point in addition to Kow and size --
  adoption-limited because MPt isn't always available.
- **Mitragotri 2002** and **Frasch 2002** are mechanistic models: fewer
  empirical parameters, more physics. Frasch is in NIOSH's calculator;
  Mitragotri is the SPT ancestor of his later 4-pathway model.
- **Magnusson 2004** predicts *Jmax*, not Kp. To convert Kp = Jmax /
  C_v,sat: needs a saturation-solubility QSPR alongside.
- **Baba 2015 / Chen 2018** represent the ML era: highest reported
  test-set R^2 but no closed-form equations (Chen reports an interpretable
  MLR alongside the SVM).

## Datasets the field re-uses

- **Flynn 1990**: ~97 compounds, in-vitro human-skin Kp from aqueous donor.
  The classical benchmark. Inter-laboratory scatter is ~15-fold for some
  compounds.
- **Wilschut 1995 compilation**: 123 Kp values, 99 compounds. Adds
  post-1990 measurements; basis for the modified Robinson model.
- **Patel 2002 / Vecchia-Bunge 2002 / EDETOX**: ~158 to ~186 Kp values for
  ~158 compounds. Most-cited 2000s benchmark.
- **Magnusson 2004 database**: 278 Jmax values; 269 with MPt + Ha.
- **Baba 2015**: 211 Kp under consistent protocol (the cleanest set).
- **HuskinDB (2020)**: open database, ~1900 entries, used by 2020+ ML
  models.

## Mandatory caveats

- All "log Kp" entries should be checked for the unit convention (cm/s vs
  cm/h; the difference is exactly log10(3600) ~= 3.56). Several papers
  publish in cm/s with the implicit -3.56 shift, others in cm/h. The
  per-model files note the convention used.
- Coefficient values quoted here come from the original papers where
  accessible, otherwise from authoritative reviews (Mitragotri 2011 Int. J.
  Pharm. 418:115-129; Patel 2002 Table 2; Kupczewska-Dobecka 2010). Where a
  paper is paywalled and only a review reproduces the equation, that is
  noted in the per-model file's "Notes" or "References" section.
- All Flynn-derived models inherit the Flynn-data caveats: aqueous donor,
  inter-laboratory scatter, steroid outliers (mostly fixed in Moss-Cronin
  2002 and Patel 2002).
- For ionisable compounds at physiological pH, every model in this list
  needs an effective-fraction-neutral correction (or use Zhang 2017
  ionisable extension, DOI: 10.1016/j.ijpharm.2017.02.064).

## Files

- [potts_guy_1992.md](./potts_guy_1992.md) — the foundational equation.
- [cleek_bunge_1993.md](./cleek_bunge_1993.md) — aqueous-epidermis correction.
- [abraham_1995.md](./abraham_1995.md) — first skin LFER.
- [wilschut_1995.md](./wilschut_1995.md) — modified Robinson, SKINPERM precursor.
- [lien_gao_1995.md](./lien_gao_1995.md) — bilinear Kow + H-bond.
- [barratt_1995.md](./barratt_1995.md) — Kow + MV + MPt.
- [potts_guy_1995.md](./potts_guy_1995.md) — H-bond refinement of 1992.
- [abraham_1999.md](./abraham_1999.md) — canonical skin LFER.
- [mitragotri_2002.md](./mitragotri_2002.md) — scaled-particle theory.
- [frasch_2002.md](./frasch_2002.md) — random-walk simulator.
- [moss_cronin_2002.md](./moss_cronin_2002.md) — Flynn re-fit, steroids fixed.
- [patel_2002.md](./patel_2002.md) — 158-compound benchmark.
- [magnusson_2004.md](./magnusson_2004.md) — predicts Jmax, not Kp.
- [ten_berge_2009.md](./ten_berge_2009.md) — SKINPERM / IH-SkinPerm.
- [baba_2015.md](./baba_2015.md) — SVR, curated 211 set.
- [chen_2018.md](./chen_2018.md) — MLR + SVM, 274-compound set.
