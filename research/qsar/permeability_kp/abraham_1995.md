---
name: abraham_1995
predicts: Kp
units: cm/s
descriptors: pi2H, alpha2H, beta2H, Vx (Abraham solvation parameters)
domain: neutral organics; original LFER fit on a Flynn-derived dataset
status: classical
---

# Abraham, Chadha, Mitchell (1995) — First Skin-LFER

## Citation
Abraham, M.H., Chadha, H.S., Mitchell, R.C. (1995). The factors that influence
skin penetration of solutes. Journal of Pharmacy and Pharmacology
47(1):8-16. DOI: 10.1111/j.2042-7158.1995.tb05725.x.

(Also relevant: Abraham, M.H., Chadha, H.S., Whiting, G.S., Mitchell, R.C.
(1994). Hydrogen bonding. 32. An analysis of water-octanol and water-alkane
partitioning. J. Pharm. Sci. 83:1085-1100.)

## Equation
log Kp [cm/s] = -5.05 - 0.59 * pi2H - 0.63 * sum(alpha2H) - 3.48 * sum(beta2H)
                + 1.79 * Vx

where:
- pi2H is the dipolarity/polarizability,
- sum(alpha2H) is the overall solute hydrogen-bond acidity,
- sum(beta2H) is the overall solute hydrogen-bond basicity,
- Vx is McGowan's characteristic volume in (cm^3 mol^-1)/100.

(Coefficients as reproduced in Mitragotri 2011 / Patel et al. 2002 review
tables; the 1995 original paper gives the same form without the excess-molar
refractivity R2 term that appears in Abraham 1999.)

## Training set
Flynn 1990 dataset, neutral organics. The exact n in the 1995 paper is in the
range 35-50 (paper-specific; the Flynn subset for which all five Abraham
descriptors are tabulated). Fit by ordinary least squares.

## Reported performance
R^2 in the range ~0.83-0.86; SD ~ 0.45-0.50 log units. The hydrogen-bond
basicity (beta2H) and solute volume (Vx) are the dominant terms.

## Validity / limitations
- Requires the five Abraham solvation descriptors. These are tabulated for
  ~5000+ compounds and computable from structure via group-contribution
  methods (Platts, Pharma Algorithms ADME Boxes), but the predicted
  descriptors are noisier than measured ones.
- Steroid sub-class is poorly fit (addressed in Abraham 1997
  J. Pharm. Pharmacol. 49:858-865 and Moss-Cronin 2002).
- Neutral solutes only; ionizable species need an effective-descriptor or
  ion-fraction adjustment (see Zhang 2017 for an ionizable extension).

## Notes
The first attempt at applying Abraham's general LFER framework to skin
permeability. Coefficient signs are mechanistically interpretable: H-bond
basicity hurts permeation strongly (beta2H coefficient is the largest in
magnitude); solute volume helps (lipid-favoured); H-bond acidity has a
modest negative effect. This paper set the template for almost every
LFER-based skin model that followed.

## References
- Abraham, M.H., Chadha, H.S., Mitchell, R.C. (1995). J. Pharm. Pharmacol.
  47(1):8-16. DOI: 10.1111/j.2042-7158.1995.tb05725.x.
- Reproduction of equation in Patel et al. (2002). Chemosphere 48:603-613.
  DOI: 10.1016/S0045-6535(02)00114-5.
- Mitragotri et al. (2011). Int. J. Pharm. 418:115-129.
  DOI: 10.1016/j.ijpharm.2011.02.023.
