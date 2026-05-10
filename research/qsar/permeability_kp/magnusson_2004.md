---
name: magnusson_2004
predicts: Jmax (maximum steady-state flux), not Kp directly
units: mol/cm^2/h (log10)
descriptors: MW (primary); Mpt and Ha as refinements
domain: in-vitro human skin from saturated aqueous donor; 87-278 compounds depending on equation form
status: classical
---

# Magnusson, Anissimov, Cross, Roberts (2004) — Molecular Size Drives Jmax

## Citation
Magnusson, B.M., Anissimov, Y.G., Cross, S.E., Roberts, M.S. (2004).
Molecular size as the main determinant of solute maximum flux across the
skin. Journal of Investigative Dermatology 122(4):993-999. PubMed: 15102090.
DOI: 10.1111/j.0022-202X.2004.22413.x.

## Equation
Three regression forms are reported. The simplest (literature dataset, n=87):

log Jmax [mol/cm^2/h] = -3.90 - 0.0190 * MW            (R^2 = 0.847)

Full database, MW only (n=278):

log Jmax [mol/cm^2/h] = -4.52 - 0.0141 * MW            (R^2 = 0.688)

Full database, MW + Mpt (melting point, degrees C) and Ha (H-bond acceptor
count) (n=269):

log Jmax = f(MW, Mpt, Ha)                              (R^2 = 0.917)

The exact 4-coefficient form for the third equation requires the original
paper; abstracts only quote the R^2 progression (R^2 = 0.760 with Mpt; 0.917
with both Mpt and Ha).

## Training set
Two datasets:
- Literature subset: 87 compounds with high-quality saturated-aqueous Jmax
  measurements through human skin in-vitro.
- Full Magnusson database: 278 compounds compiled from 36 references.
- Reduced (Mpt+Ha available): 269 compounds.

## Reported performance
See equations. R^2 from 0.69 (MW only, n=278) to 0.92 (MW+Mpt+Ha, n=269).

## Validity / limitations
- Predicts *Jmax*, not Kp. Conversion: Kp = Jmax / Cv,sat, where Cv,sat is
  the donor saturation concentration -- so Magnusson's Jmax is the
  permeability-relevant quantity for *saturated* aqueous formulations only.
- For risk assessment from ambient (sub-saturated) exposures, Kp is more
  useful and one needs to solve back through Cv,sat (which has its own QSPR
  uncertainty).
- MW dominance is dataset-specific; for organic-vehicle exposures, both
  partition and diffusion change and the simple MW dependence breaks down.

## Notes
Major shift in thinking: for saturated aqueous donors the *flux* is
dominated by MW (R^2 = 0.85 with one descriptor), whereas Kp needs at least
two (Kow + MW). This makes sense mechanistically: in the Kp = D*K/h
decomposition, the diffusion coefficient D scales with MW while the
partition K scales with Kow; for saturated donors, K and Cv,sat enter with
opposite signs and partly cancel, leaving MW (via D) as the dominant driver.
Updated by Zhang 2009 (Pharm. Res. 26:1974-1985) which reports a bilinear,
non-monotone Jmax-Kow relationship for similar-size molecules.

## References
- Magnusson, B.M., Anissimov, Y.G., Cross, S.E., Roberts, M.S. (2004).
  J. Invest. Dermatol. 122(4):993-999. DOI: 10.1111/j.0022-202X.2004.22413.x.
- Magnusson, B.M., Pugh, W.J., Roberts, M.S. (2004). Simple rules defining
  the potential of compounds for transdermal delivery or toxicity. Pharm.
  Res. 21(6):1047-1054. (Companion paper.)
- Zhang, Q. et al. (2009). Skin solubility determines maximum transepidermal
  flux for similar size molecules. Pharm. Res. 26:1974-1985.
- Roberts, M.S. et al. (2022). An updated database of human maximum skin
  fluxes and epidermal permeability coefficients for drugs, xenobiotics, and
  other solutes applied as aqueous solutions. PubMed: 35599823.
