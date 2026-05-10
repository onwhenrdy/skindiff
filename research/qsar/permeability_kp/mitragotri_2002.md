---
name: mitragotri_2002
predicts: Kp (mechanistic via scaled-particle theory of lipid-bilayer diffusion)
units: cm/s
descriptors: Kow, r (solute radius in Angstroms)
domain: small hydrophobic solutes diffusing through SC lipid; assumes lipid pathway only
status: classical
---

# Mitragotri (2002) — Scaled-Particle-Theory Mechanistic Model

## Citation
Mitragotri, S. (2002). A theoretical analysis of permeation of small
hydrophobic solutes across the stratum corneum based on Scaled Particle
Theory. Journal of Pharmaceutical Sciences 91(3):744-752. PubMed: 11920759.
DOI: 10.1002/jps.10071.

## Equation
Kp [cm/s] ~ 5.6e-6 * Kow^0.7 * exp(-0.46 * r^2)

(Form as reproduced in Patel et al. 2002 Table 2; r is the solute radius in
Angstroms calculated from molecular volume.)

The diffusion coefficient through the lipid bilayer in the SPT framework is

D_lip [cm^2/s] = 2e-5 * exp(-0.46 * r^2)

(Mitragotri 2011 Eq. 32.) This is combined with a Kow^0.7 partition law into
the bilayer to give the Kp expression above.

## Training set
Coefficients (the 5.6e-6 prefactor, the 0.7 partition exponent, the 0.46 r^2
constant) are derived theoretically from scaled-particle theory of lipid
bilayers, then anchored to a small subset of the Flynn data (~ 30-50
compounds) for the prefactor. Not a regression in the OLS sense.

## Reported performance
On the Flynn-style validation set: R^2 ~ 0.70, MAE ~ 0.09 log units
(comparison in Lian-Chen 2008). Performs better than Potts-Guy on small
hydrophobic solutes; worse on hydrophilic solutes because the model assumes
the lipid pathway is the only one.

## Validity / limitations
- Lipid-pathway-only: under-predicts Kp for hydrophilic solutes that go
  through the corneocyte / aqueous pore pathway (water, urea, glycerol).
- Requires an estimate of solute radius r in Angstroms; if computed from
  molecular volume V via r = (3V / 4*pi)^(1/3), small estimation errors
  amplify in the exp(-0.46 r^2) factor.
- Single-membrane assumption: no aqueous-epidermis correction.
- Highly accurate in its domain of validity (small, neutral, lipophilic).

## Notes
Different philosophy from the empirical regressions: parameters from
statistical mechanics, not OLS. Notable for being one of the few models that
*can* explain the Kp pattern of small hydrophobics with so few free
parameters. Mitragotri's later "four-pathway" 2003 extension
(J. Control. Release 86:69-92) adds porous, shunt, and corneocyte routes,
covering hydrophilic solutes too -- but at the cost of being more like a
simulator than a QSAR.

## References
- Mitragotri, S. (2002). J. Pharm. Sci. 91(3):744-752.
  DOI: 10.1002/jps.10071.
- Mitragotri, S. (2003). Modeling skin permeability to hydrophilic and
  hydrophobic solutes based on four permeation pathways. J. Control. Release
  86(1):69-92. (Extension to four pathways.)
- Mitragotri et al. (2011). Mathematical models of skin permeability: an
  overview. Int. J. Pharm. 418:115-129, Eqs. 32 and surrounding.
  DOI: 10.1016/j.ijpharm.2011.02.023.
