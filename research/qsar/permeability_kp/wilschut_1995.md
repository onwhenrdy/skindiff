---
name: wilschut_1995
predicts: Kp
units: cm/h
descriptors: logKow, MW
domain: aqueous vehicle, neutral organics; explicitly accounts for aqueous-epidermis resistance
status: classical
---

# Wilschut, ten Berge, Robinson, McKone (1995) — Modified Robinson Model

## Citation
Wilschut, A., ten Berge, W.F., Robinson, P.J., McKone, T.E. (1995). Estimating
skin permeation. The validation of five mathematical skin permeation models.
Chemosphere 30(7):1275-1296. PubMed: 7749723.
DOI: 10.1016/0045-6535(95)00023-2.

## Equation
The paper validated five published Kp models against an aqueous-vehicle
permeability dataset (123 Kp values, 99 compounds, in-vitro human skin), and
produced a *modified Robinson* form as the best performer. Its core idea is a
parallel-pathway split with MW^0.5 (rather than MW) in the diffusivity term,
giving better behaviour at MW extremes. The modified Robinson model expresses
Kp as a hindered-aqueous + lipoidal-pathway sum; the published coefficients
are part of the SKINPERM family and were re-fit by ten Berge (2009).

The paper does not give a closed-form one-line equation suitable for
copy-paste; the practical implementation is the SKINPERM workbook (ten Berge
2009), which uses the same parallel-pathway philosophy.

For comparison, Wilschut's group also re-evaluated the Guy-Potts (1992) form
by re-fitting on their 123-point dataset, obtaining (Patel-style restatement,
Kupczewska-Dobecka 2010 Eq. 5):

log Kp [cm/h] = 0.781 * log Kow - 0.01115 * MW - 2.19   (n=123)

## Training set
n = 123 measured Kp values for 99 compounds applied in vitro to human skin in
aqueous solution. Compiled from the literature; covers monoaromatic
hydrocarbons, halogenated aliphatics, phenols, steroids, and small drugs.

## Reported performance
Coefficients of all five models (Brown-Rossi, Fiserova-Bergerova, McKone-Howd,
Guy-Potts-style, modified Robinson) were re-estimated by non-linear multiple
regression. The modified Robinson model gave the *smallest residual variance*
across the 123-point set; it predicts highly hydrophilic and highly lipophilic
compounds more accurately than the Guy-Potts form. Specific R^2 / RSE values
were reported per-model in the original Chemosphere tables (paywalled).

## Validity / limitations
- Aqueous donor only.
- Inherits the inter-laboratory scatter of the Wilschut compilation
  (~15-fold range of reported Kp for the same compound across labs).
- The "modified" part is the MW^0.5 substitution -- this captures
  hindered-pore behaviour for water-pathway-dominant compounds but is
  empirical, not derived from first principles.
- Outliers similar to Flynn dataset (steroids, certain lipophilics).

## Notes
Wilschut 1995 is the methodological backbone of ten Berge's IH-SkinPerm tool
(2009 onwards) and the EPA/NIOSH "Modified Robinson" calculator option. It
matters historically because it was the first study to explicitly *validate*
multiple competing Kp models on a common dataset rather than fit a single new
one. The aqueous-pathway emphasis distinguishes it from the lipid-only
Potts-Guy formulation.

## References
- Wilschut, A., ten Berge, W.F., Robinson, P.J., McKone, T.E. (1995).
  Estimating skin permeation: The validation of five mathematical skin
  permeation models. Chemosphere 30(7):1275-1296.
  DOI: 10.1016/0045-6535(95)00023-2.
- Robinson, P.J. (1993). A composite model for predicting dermal penetration
  in vivo. Procter & Gamble Human and Environmental Safety Division
  technical report; cited in Wilschut 1995.
- Kupczewska-Dobecka et al. (2010). Calculating the dermal flux of chemicals
  with OELs based on their molecular structure. Environ. Toxicol. Pharmacol.
  30:95-102, see Eq. 5. DOI: 10.1016/j.etap.2010.06.005.
- NIOSH Skin Permeation Calculator (Modified Robinson option):
  https://www.cdc.gov/niosh/skin-exposure/resources/skin-permeation-calculator.html
