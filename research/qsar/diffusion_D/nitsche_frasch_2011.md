---
name: nitsche_frasch_2011
predicts: D
layer: corneocyte
units: cm^2/s (effective)
descriptors: D_aq, partition coefficient K_kw, binding rate constants
domain: SC corneocyte interior with reversible solute-keratin binding
status: modern
---

# Nitsche–Frasch 2011 — Diffusion + reversible binding in heterogeneous SC

## Citation
Nitsche JM, Frasch HF. "Dynamics of diffusion with reversible binding in
microscopically heterogeneous membranes: General theory and applications
to dermal penetration." Chem Eng Sci. 2011;66(10):2019–2041.
DOI: 10.1016/j.ces.2011.01.024.
(Extended in Nitsche & Kasting 2022, J Pharm Sci 111(7):1923–1935.
DOI 10.1016/j.xphs.2022.02.014.)

## Equation
The 2011 paper is a *theoretical framework* deriving the effective
membrane diffusion coefficient when solutes undergo fast reversible
binding to immobile sites (here keratin in corneocytes). Two limits:

- **Fast binding limit (local equilibrium):**
```
D_eff = D_free / (1 + R_B)
```
where R_B = (free + bound)/free = 1 + K_b * c_sites is a partition-derived
retardation factor. Recovers the Crank "diffusion with linear isotherm"
result.

- **Slow binding limit:**
```
D_eff -> D_free  (binding stalls behind the diffusion front; bound mass
                  is left behind as a tail)
```

The intermediate-rate regime requires the full coupled PDE solution
(memory kernel form, also derivable as an Onsager-de Groot effective
medium).

For SC corneocytes, the practical formula (combining hindered diffusion
+ binding) becomes:

```
D_cor_eff = D_aq * H(lambda) * (1 - phi_f) / (1 + R_B)
```

with H(lambda) the Bungay-Brenner polynomial and R_B from keratin-binding
data (Hansen 2011 keratin-binding database, ~150 compounds).

## Training set
Theoretical, no fitted training set; the 2022 follow-up regresses against
the Hansen extended-keratin-binding database (~150 compounds, MW 30–400)
to derive forward/reverse binding rate constants for use in transient
SC simulations.

## Reported performance
2022 follow-up: predicted vs. measured macroscopic forward binding
rate constants within factor 2–3 across the database; lipophilicity
maximum near log K_ow ~ 2 (more than Hansen's data suggested).

## Validity / limitations
- The framework assumes linear (concentration-independent) binding
  isotherms; non-linear binding requires extension.
- Practical D_cor_eff predictions are coupled to keratin binding QSARs
  (e.g. Yamaguchi 2009, Hansen 2011) for R_B.
- For transient (finite-dose) skin transport the slow-binding regime
  matters; standard "lumped D" models can be biased.

## Notes
Underpins the modern view that "D in corneocyte" is best computed as
**aqueous D × hindrance × (1 + binding retardation)^-1**, with binding
the new ingredient versus older multiphase models. Important when
computing transient release (skin reservoir effect) — at steady state
the binding retardation drops out of P but lengthens the lag time.

## References
- DOI 10.1016/j.ces.2011.01.024 (Nitsche & Frasch 2011).
- DOI 10.1016/j.xphs.2022.02.014 (Nitsche & Kasting 2022 transient
  binding framework).
- Hansen S et al. J Pharm Sci. 2011;100(5):1655–1668.
  DOI 10.1002/jps.22397 (extended keratin binding dataset).
