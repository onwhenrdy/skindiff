---
name: barratt_1995
predicts: Kp
units: cm/s
descriptors: logKow, MV (molecular volume), MPt (melting point)
domain: neutral organics, Flynn-derived training set
status: classical
---

# Barratt (1995) — logKow + Volume + Melting-Point Model

## Citation
Barratt, M.D. (1995). Quantitative structure-activity relationships for skin
permeability. Toxicology in Vitro 9(1):27-37.
DOI: 10.1016/0887-2333(94)00190-6.

## Equation
log Kp [cm/s] = -5.9163 + 0.82 * log Kow - 0.0093 * MV - 0.039 * MPt

where:
- MV is molecular volume (cm^3/mol),
- MPt is melting point (degrees C).

(Coefficients as reproduced in Patel et al. 2002 Table 2.)

## Training set
Subset of Flynn 1990, restricted to compounds for which a reliable melting
point was available. The original paper reports n in the range 60-65 with
two outlier exclusions.

## Reported performance
R^2 ~ 0.90 reported in the original paper. The melting-point term reflects
the Higuchi-Yalkowsky-style observation that crystal-lattice energy modulates
solubility-driven flux; including it captures additional variance not
explained by Kow + size alone.

## Validity / limitations
- Requires experimental melting point (or a reliable QSPR for MPt). For
  liquids (MPt < 25 degrees C), the convention varies -- some users set MPt
  to room temperature, others use the freezing point.
- MV (not MW) is the size descriptor; MV is well-tabulated for most small
  molecules but conventions differ (McGowan vs van der Waals vs LeBas).
- Does not address the lipophilic-bound aqueous-epidermis limit
  (Cleek-Bunge correction recommended for log Kow > ~3).

## Notes
A useful "second-tier" QSAR: when MPt is available it tightens prediction
versus Potts-Guy. The MPt term is unusual in the Kp literature -- most other
models avoid it because reliable melting points are scarce. Cited mainly as a
reference equation in evaluation studies (Patel 2002 Table 2; Lian-Chen 2008).

## References
- Barratt, M.D. (1995). Quantitative structure-activity relationships for
  skin permeability. Toxicol. in Vitro 9(1):27-37.
  DOI: 10.1016/0887-2333(94)00190-6.
- Patel, H., ten Berge, W., Cronin, M.T.D. (2002). Chemosphere 48:603-613,
  Table 2.
