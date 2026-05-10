---
name: patel_2002
predicts: Kp
units: cm/s
descriptors: logKow, MW, ABSQon (sum of absolute charges on N, O), SsssCH (E-state index)
domain: neutral organics; outlier-pruned 158-compound human-skin set
status: classical
---

# Patel, ten Berge, Cronin (2002) — Hydrophobicity + Size + H-Bond Refinement

## Citation
Patel, H., ten Berge, W., Cronin, M.T.D. (2002). Quantitative structure-
activity relationships (QSARs) for the prediction of skin permeation of
exogenous chemicals. Chemosphere 48(6):603-613. PubMed: 12143935.
DOI: 10.1016/S0045-6535(02)00114-5.

## Equation
log Kp [cm/s] = -2.47 + 0.681 * log Kow - 0.00653 * MW
                - 0.284 * ABSQon - 0.268 * SsssCH

where:
- ABSQon = sum of absolute partial charges on N and O atoms (proxy for
  H-bond capacity),
- SsssCH = sum of E-state values for >CH< groups (proxy for branching).

(Coefficients from Patel et al. 2002 Eq. 4 / reproduced in Mitragotri 2011.)

## Training set
n = 158 compounds initially; reduced to ~143 after outlier removal
(hydrocortisone derivatives). In-vitro Kp through excised human skin from
aqueous donor. Largest single-curated Kp set as of 2002.

## Reported performance
R^2 = 0.90 on the outlier-pruned set; SE ~ 0.4 log units. The two
electronic descriptors (ABSQon, SsssCH) add ~0.05 to R^2 over a
Potts-Guy-style two-descriptor fit on the same set, and are mechanistically
interpretable: ABSQon captures H-bond donation/acceptance; SsssCH captures
hydrocarbon-skeleton branching effects on diffusion.

## Validity / limitations
- ABSQon and SsssCH are not as widely-tabulated as Kow / MW; need a
  cheminformatics toolkit (E-state via Hall-Kier; partial charges via
  Gasteiger or DFT).
- Steroids are explicitly excluded as outliers (15 compounds dropped).
  Different from Moss-Cronin's strategy of replacing steroid Kp values.
- Same Flynn-style donor/skin caveats.

## Notes
The 158-compound dataset became the de-facto benchmark for the next decade of
Kp QSARs (Khajeh & Modarress, Neely, Lam, etc.). The four-descriptor form is
the best Kp prediction reported for any classical regression as of 2002.
Patel et al. also provide simpler 2-descriptor and 3-descriptor variants in
the same paper for users who can't compute the electronic indices.

## References
- Patel, H., ten Berge, W., Cronin, M.T.D. (2002). Chemosphere 48(6):603-613.
  DOI: 10.1016/S0045-6535(02)00114-5.
- Mitragotri et al. (2011). Mathematical models of skin permeability: an
  overview. Int. J. Pharm. 418:115-129, Table 1.
  DOI: 10.1016/j.ijpharm.2011.02.023.
