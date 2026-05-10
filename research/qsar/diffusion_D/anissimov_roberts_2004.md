---
name: anissimov_roberts_2004
predicts: D
layer: SC
units: cm^2/s (relative position-dependent)
descriptors: depth z within SC; relative D(z)/D(0)
domain: variable diffusion coefficient through SC depth; not a structural QSAR
status: classical
---

# Anissimov–Roberts 2004 — Variable diffusion-coefficient SC profiles

## Citation
Anissimov YG, Roberts MS. "Diffusion modeling of percutaneous absorption
kinetics: 3. Variable diffusion and partition coefficients, consequences
for stratum corneum depth profiles and desorption kinetics." J Pharm Sci.
2004;93(2):470–487. DOI: 10.1002/jps.10567. PMID 14705201.
(Earlier installments: Part 1, J Pharm Sci 1999;88(11):1201–1209;
Part 2, J Pharm Sci 2001;90(4):504–520.)

## Equation
Anissimov–Roberts is *not* a structural QSAR for D in terms of MW; it
is a class of **depth-dependent D(z) functional forms** for SC, fit
post-hoc to desorption / depth-profile data:

```
Linear:        D(z) = D_0 * (1 + a * z / L)
Exponential:   D(z) = D_0 * exp(-b * z / L)
Two-slab:      D(z) = D_top  for z < L/2;  D_bot  for z > L/2
```

with L = SC thickness, z depth from skin surface, D_0 surface diffusivity.
Parameters (a, b) or (D_top / D_bot ratio) are fit to experimental
concentration-depth profiles (CDP) and desorption curves.

A key result: **steady-state permeation P is insensitive to the *direction*
of D and K asymmetry** — a result of solving the steady transport
equation. Under finite-dose / desorption conditions, however, the
asymmetry strongly shapes the depth profile and the lag time.

## Training set
Examples in the paper: published D(z) and K(z) profiles for nicotine,
salicylic acid, and steroids in human SC (Anigbogu, Roberts groups);
self-consistency tests against synthetic data.

## Reported performance
Demonstrates that conventional homogeneous-D analysis can misestimate
SC reservoir capacity and lag time by 50–200% if D(z) varies by more
than 10x with depth — quantitative result depends on solute.

## Validity / limitations
- Provides *functional forms*, not *parameter QSARs*. Each new compound
  needs its own fit; no MW or log K_ow scaling is given for D_0, a, b.
- Only useful when CDP or desorption data exist. For pure forward
  prediction from chemical structure, fall back to Mitragotri or
  Wang–Kasting.
- The two-slab and linear forms are first-pass approximations to the
  underlying microheterogeneity captured better by Wang–Kasting
  bricks-and-mortar.

## Notes
This work is foundational to the modern view that "D in SC" is not a
single number but a depth-dependent property reflecting SC delipidation
gradient and the corneocyte size/orientation gradient. For PBPK-style
average D it can be averaged out; for finite-dose toxicology it cannot.

## References
- DOI 10.1002/jps.10567 (Anissimov & Roberts 2004, Part 3).
- DOI 10.1002/(SICI)1520-6017(199911)88:11<1201::AID-JPS17>3.0.CO;2-7
  (Part 1, 1999).
- DOI 10.1002/jps.1011 (Part 2, 2001).
