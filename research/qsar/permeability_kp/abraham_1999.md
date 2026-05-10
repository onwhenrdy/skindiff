---
name: abraham_1999
predicts: Kp
units: cm/s
descriptors: R2, pi2H, alpha2H, beta2H, Vx (Abraham solvation parameters)
domain: neutral organics; expanded Flynn-derived training set
status: classical
---

# Abraham, Chadha, Martins, Mitchell, Bradbury, Gratton (1999) — Updated Skin LFER

## Citation
Abraham, M.H., Chadha, H.S., Martins, F., Mitchell, R.C., Bradbury, S.P.,
Gratton, J.A. (1999). Hydrogen bonding part 46: A review of the correlation
and prediction of transport properties by an LFER method: physicochemical
properties, brain penetration and skin permeability. Pesticide Science
55(1):78-88. DOI: 10.1002/(SICI)1096-9063(199901)55:1<78::AID-PS859>3.0.CO;2-Y.

Companion: Abraham, M.H., Martins, F., Mitchell, R.C. (1997). Algorithms for
skin permeability using hydrogen bond descriptors: the problem of steroids.
J. Pharm. Pharmacol. 49(9):858-865. PubMed: 9306252.

## Equation
log Kp [cm/s] = -5.13 + 0.44 * R2 - 0.49 * pi2H - 1.48 * sum(alpha2H)
                - 3.44 * sum(beta2H) + 1.94 * Vx

where:
- R2 is the excess molar refraction ((cm^3 mol^-1)/10),
- pi2H, alpha2H, beta2H are dipolarity/polarizability and H-bond
  acidity/basicity,
- Vx is McGowan's characteristic volume ((cm^3 mol^-1)/100).

(Coefficients as reproduced in Patel et al. 2002 Table 2.)

## Training set
n ~= 53 in the 1999 paper for the basic skin equation, expanded to ~119
solutes in the 2004 J. Pharm. Sci. follow-up (Abraham & Martins). Drawn from
Flynn 1990 plus newer steroid/drug data.

## Reported performance
R^2 ~= 0.86-0.87, SD ~ 0.43-0.45 log units. The H-bond acidity coefficient
(alpha2H) becomes much more negative than in the 1995 fit (-1.48 vs -0.63),
which the authors attribute to better-curated data and the addition of R2.
Adding R2 marginally improves the fit and aligns the model with Abraham's
broader transport-property work.

## Validity / limitations
- Five descriptors -- needs all of them (R2, pi2H, sum-alpha2H, sum-beta2H,
  Vx). Abraham descriptor databases (ADME Boxes, UCL group databases) cover
  most small organics, but bespoke compounds need group-contribution
  estimation.
- Neutral-solute regression. Use Zhang 2017 or apply
  effective-fraction-neutral correction for ionisable compounds.
- Steroid problem partially resolved (Abraham 1997) by treating sterol
  scaffolds as a separate fit; the 1999 review aggregates the corrected
  data.

## Notes
The "canonical" skin LFER cited in most reviews. Abraham & Martins 2004
(J. Pharm. Sci. 93:1508-1523, DOI: 10.1002/jps.20094) extends the dataset to
n=119, achieving R^2 = 0.832 / SD = 0.46 log units; the 2004 coefficients are
similar in sign and magnitude to the 1999 set. Used as the comparator in
Kupczewska-Dobecka 2010 (which calls this "ABR" model and finds significantly
different mean log Kp from Potts-Guy at p<0.05).

## References
- Abraham, M.H. et al. (1999). Pesticide Science 55(1):78-88.
- Abraham, M.H., Martins, F., Mitchell, R.C. (1997). J. Pharm. Pharmacol.
  49:858-865. PubMed: 9306252.
- Abraham, M.H., Martins, F. (2004). Human skin permeation and partition:
  general linear free-energy relationship analyses. J. Pharm. Sci.
  93(6):1508-1523. PubMed: 15124209. DOI: 10.1002/jps.20094.
- Patel, H., ten Berge, W., Cronin, M.T.D. (2002). Chemosphere 48:603-613.
  DOI: 10.1016/S0045-6535(02)00114-5 (reproduces 1999 equation in Table 2).
