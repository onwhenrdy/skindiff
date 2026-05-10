---
name: cleek_bunge_1993
predicts: Kp (with viable-epidermis correction)
units: cm/s (or cm/h, dimensionally consistent with input Kp)
descriptors: Kp_sc (from a base QSAR like Potts-Guy), MW
domain: lipophilic solutes (log Kow > ~3) where the viable epidermis becomes rate-limiting
status: classical
---

# Cleek and Bunge (1993) — Aqueous-Boundary-Layer / Viable-Epidermis Correction

## Citation
Cleek, R.L., Bunge, A.L. (1993). A new method for estimating dermal absorption
from chemical exposure. 1. General approach. Pharmaceutical Research
10(4):497-506. PubMed: 8483831. DOI: 10.1023/A:1018981515480.

Companion: Bunge, A.L., Cleek, R.L. (1995). A new method for estimating dermal
absorption from chemical exposure. 2. Effect of molecular weight and
octanol-water partitioning. Pharm. Res. 12(1):88-95. PubMed: 7724493.

## Equation
Kp_adjusted = Kp_sc / (1 + B)

with B = (Kp_sc * sqrt(MW)) / 2.6   (Kp in cm/h, MW in Da)

Equivalently, if Kp_sc is the Potts-Guy 1992 prediction in cm/s:

Kp_adj [cm/s] = Kp_sc / (1 + 1400 * Kp_sc * sqrt(MW))

(Mitragotri 2011 Eq. 14 gives the cm/s form with the 1400 prefactor; the
cm/h form uses the 2.6 g/(cm^2 h) effective aqueous-epidermis permeability.)

## Training set
This is a *correction*, not a regression QSAR -- it has no compounds of its
own. It is derived from a series-resistance argument: total skin permeability
is the harmonic sum of SC permeability (from Potts-Guy or any other Kp_sc
QSAR) and an aqueous-epidermis permeability that scales as ~2.6/sqrt(MW)
cm/(g/h). The 2.6 factor comes from a literature value for an aqueous-pore
permeability through the viable epidermis.

## Reported performance
Not a stand-alone fit -- evaluated by demonstrating the correction tightens
the experiment-vs-prediction ratio for ~14 highly lipophilic compounds
(steroids, hexylnicotinate, etc.) to within an order of magnitude (Mitragotri
2011, Fig. 2; Guy 2010).

## Validity / limitations
- Only matters when B is non-negligible -- i.e. for highly lipophilic
  compounds (log Kow > ~3) with smallish MW. For polar drugs B << 1 and the
  correction reduces to Kp_adj ~= Kp_sc.
- Inherits all limitations of the underlying Kp_sc model (usually Potts-Guy).
- Assumes a fixed "effective viable-epidermis permeability" of ~2.6/sqrt(MW)
  cm h^-1 (g/mol)^-1/2; in reality this scaling has noticeable error and was
  derived from a small set of solutes diffusing through hydrated dermis.
- A series-resistance model -- it cannot describe transient effects, only
  shifts the steady-state Kp downward.

## Notes
The standard "next step after Potts-Guy" in occupational risk assessment.
US EPA's 2007 dermal-exposure guidance and NIOSH's Skin Permeation Calculator
both apply this correction on top of Potts-Guy. ten Berge's IH-SkinPerm
applies the same idea via a separate viable-epidermis permeability term.
Trivial to implement as a wrapper around any Kp_sc model.

## References
- Cleek, R.L., Bunge, A.L. (1993). Pharm. Res. 10(4):497-506.
  DOI: 10.1023/A:1018981515480.
- Bunge, A.L., Cleek, R.L. (1995). Pharm. Res. 12(1):88-95. PubMed: 7724493.
- Mitragotri et al. (2011). Mathematical models of skin permeability: an
  overview. Int. J. Pharm. 418:115-129, see Eq. 14.
  DOI: 10.1016/j.ijpharm.2011.02.023.
- Guy, R.H. (2010). Predicting the rate and extent of fragrance chemical
  absorption into and through the skin. Chem. Res. Toxicol. 23:864-870.
