---
name: frasch_2002
predicts: Kp (via random-walk simulation of biphasic SC)
units: cm/s
descriptors: MW, logKow (used to derive effective D and effective path length)
domain: stratum corneum random-walk through brick-and-mortar; Flynn-derived training set
status: classical
---

# Frasch (2002) — Random Walk Model of Skin Permeation

## Citation
Frasch, H.F. (2002). A random walk model of skin permeation. Risk Analysis
22(2):265-276. PubMed: 12022675. DOI: 10.1111/0272-4332.00024.

## Equation
Frasch's model is *mechanistic*, not a closed-form regression: a 2-D random
walk through a corneocyte-and-lipid stratum-corneum geometry is simulated;
effective diffusivity D_eff(MW, logKow) and effective path length h_eff are
extracted; then

Kp [cm/s] = (D_eff * K_sc/v) / h_eff

with the SC/vehicle partition coefficient K_sc/v fitted from logKow via a
power law, and the effective diffusivity from a free-volume MW dependence.
The actual coefficient values are inside the simulator code (and the paper's
fitted constants for the corneocyte-vs-lipid partition split) rather than in
a single one-line equation.

## Training set
n = 94 compounds from the Flynn 1990 database. Coefficients for the
underlying lipid- and corneocyte-pathway permeabilities are fitted to
minimize residuals against this set.

## Reported performance
R^2 = 0.84, SE = 0.0076 (log Kp units? -- ambiguous in the abstract;
Mitragotri 2011 quotes the same R^2). F-statistic = 154. This is the highest
R^2 achieved on the Flynn set by any model that doesn't use Abraham
descriptors -- a notable improvement over Potts-Guy's R^2 = 0.67 with the
same input descriptors (MW + logKow).

## Validity / limitations
- Simulator-based: requires running the random-walk code rather than
  evaluating a closed-form expression. CDC NIOSH packages it as the "Frasch"
  option in the Skin Permeation Calculator.
- Same Flynn-data caveats (steroids included; inter-laboratory scatter).
- Mechanistic parameters (corneocyte aspect ratio, lipid thickness) baked
  into the geometry; not user-tunable in the calculator.

## Notes
Demonstrates that the Flynn-set residual variance Potts-Guy ascribes to
"data noise" is in fact partially structural -- a more realistic SC geometry
captures it. Frasch's approach is a precursor to the brick-and-mortar
mechanistic models (Wang-Kasting-Nitsche 2006/2007, Naegel et al. 2008/2009)
that compute Kp from first principles.

## References
- Frasch, H.F. (2002). Risk Anal. 22(2):265-276.
  DOI: 10.1111/0272-4332.00024.
- Errata in Frasch & Barbero (2003). J. Pharm. Sci. 92:2196-2207.
- NIOSH Skin Permeation Calculator -- "Frasch" model option:
  https://www.cdc.gov/niosh/skin-exposure/resources/skin-permeation-calculator.html
