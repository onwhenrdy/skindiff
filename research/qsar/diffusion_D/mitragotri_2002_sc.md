---
name: mitragotri_2002_sc
predicts: D
layer: SC_lipid
units: cm^2/s
descriptors: molecular_radius_r (Angstrom), MW
domain: small (MW < 500 Da) hydrophobic solutes; lipid-bilayer pathway
status: classical
---

# Mitragotri 2002 — Scaled Particle Theory for SC lipid bilayers

## Citation
Mitragotri S. "A theoretical analysis of permeation of small hydrophobic solutes
across the stratum corneum based on Scaled Particle Theory." J Pharm Sci.
2002 Mar;91(3):744–752. DOI: 10.1002/jps.10031. PMID 11920759.

## Equation
Mitragotri factors steady-state P into K and D contributions. The lipid-bilayer
diffusion coefficient is the exponential-in-area term in the lumped permeability:

```
P [cm/s] = 5.6e-6 * K_ow^0.7 * exp(-0.46 * r^2)         (r in Angstrom)
```

Decomposing P = K_lip * D_lip / h_lip with K_lip ~ K_ow^0.7 gives the
diffusion-coefficient scaling

```
D_lip(r) = D_0 * exp(-0.46 * r^2)         (r in Angstrom, A^2)
```

i.e. an **exponential dependence on solute cross-sectional area**, derived
from scaled particle theory (work of cavity formation in a fluid of lipid
chain segments). For typical small solutes (r ~ 2–4 A) D_lip lies in the
range ~1e-9 to ~1e-7 cm^2/s. Translating r^2 to MW via spherical solute
gives an effective scaling D ~ exp(-c * MW^(2/3)) — distinct from sqrt(MW)
or 1/3-power forms used elsewhere.

## Training set
Tested against ~60 lipophilic permeants from Flynn's compiled human-skin
permeability database (the same dataset Potts–Guy used). MW range
~30–500 Da, log K_ow range −2 to 6.

## Reported performance
P predictions correlate with measured P over four log decades; the model
later won Lian et al. 2008's seven-model bake-off as the best mechanistic
QSAR (lowest mean squared error vs. 124-compound database).

## Validity / limitations
Only the lipid-bilayer route. Free-volume limit (D > 0) is not embedded —
the exponential keeps decaying without bound at large r, so the model
underestimates D for MW > ~500 Da where transport collapses anyway.
Hydrophilic / charged solutes are out-of-domain; the four-pathway 2003
extension covers them via aqueous-pore terms.

## Notes
The 0.46 prefactor in the exponent is fit; the exponential **form** is
predicted from cavity-creation thermodynamics with Hamiltonian parameters
of skin lipid chains (chain length, packing density). Unlike sqrt(MW)
hopping forms, this is an *area*-scaling, which becomes critical in
the 200–500 Da window where models diverge by orders of magnitude.

## References
- DOI 10.1002/jps.10031 (Mitragotri 2002, PubMed 11920759).
- Lian, Chen, Han. J Pharm Sci 2008, 97:584–598. DOI 10.1002/jps.21074
  (model evaluation, 124 compounds).
