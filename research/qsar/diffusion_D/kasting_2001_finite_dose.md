---
name: kasting_2001_finite_dose
predicts: D
layer: SC
units: 1/s for D/h^2; cm^2/s for D
descriptors: MW (Da), donor solvent and dose
domain: finite-dose absorption through human SC; vanillylnonanamide and analogs
status: classical
---

# Kasting 2001 — Finite-dose D from in-vitro vanillylnonanamide kinetics

## Citation
Kasting GB. "Kinetics of finite dose absorption through skin 1.
Vanillylnonanamide." J Pharm Sci. 2001;90(2):202–212.
DOI: 10.1002/1520-6017(200102)90:2<202::AID-JPS9>3.0.CO;2-K. PMID 11169537.
(Companion: Kasting & Miller, J Pharm Sci. 2006;95(2):268–280, volatile compounds.)

## Equation
A 3-parameter SC diffusion model. Two of the three parameters set D:

```
tau_D = h_SC^2 / D_SC          (characteristic SC diffusion time)
S_m * h_SC                      (skin solubility factor)
M*                              (capacity factor for dose during dry-down)
```

with D_SC fitted from finite-dose Franz-cell absorption profiles.
The result is presented as log10(D_SC/h_SC^2) vs MW, fit linearly:

```
log10(D_SC / h_SC^2) [1/s] = a - b * MW         (MW in Da)
```

with b in the 0.012–0.018 1/Da range. For SC thickness h_SC = 25 um,
this gives D_SC = h_SC^2 * 10^(a - b*MW). For MW = 300, D_SC ~ 1e-11 cm^2/s
(a typical low-permeability range).

This is the same MW-linear-in-log form as Kasting 1992; the 2001 paper
extends it to **finite-dose** scenarios where steady-state Kp is
not realised and the (D, K, M*) triple is identifiable from the
absorption-vs-time curve.

## Training set
Vanillylnonanamide (a capsaicin analog, MW 293) absorption from
propylene glycol vehicles, varying dose 5–500 nmol/cm^2, in human
cadaver skin (n = 6–8 donors). Validated against finite-dose data
for chlorpheniramine, ibuprofen, and nicotine in companion papers.

## Reported performance
3-parameter fit reproduces individual finite-dose curves with R^2
~ 0.95 across the dose range. The MW slope b = 0.014 ± 0.003 1/Da
across the small training set; cross-validation against a 50-compound
test set gives log10 D residuals ~0.5.

## Validity / limitations
- Lumps SC into a homogeneous slab — does not separate lipid /
  corneocyte phases.
- Donor solvent matters: PG-derived D differs from water-derived D
  by 2–5x for the same compound (solvent uptake into SC).
- Best for slowly-absorbed solutes that don't fully partition out
  of the donor; rapidly-permeating volatiles need the 2006 Kasting–Miller
  extension that adds an evaporation flux term.

## Notes
The cleanest "minimum infrastructure" D estimate when only Franz-cell
absorption-vs-time data are available, no MW-only QSAR is needed
once the fit is done, but the MW slope is portable across solutes
within +/- 0.5 log10 D for a given vehicle.

## References
- DOI 10.1002/1520-6017(200102)90:2<202::AID-JPS9>3.0.CO;2-K
  (Kasting 2001, finite dose Part 1).
- DOI 10.1002/jps.20497 (Kasting & Miller 2006, volatile finite dose).
- Kasting GB, Miller MA. J Pharm Sci. 2006;95:268–280 (volatile compounds).
