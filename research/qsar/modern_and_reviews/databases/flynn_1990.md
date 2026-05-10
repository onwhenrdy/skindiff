---
name: flynn_1990
type: database
url: not directly hosted; reproduced as a table in many subsequent papers (Potts-Guy 1992, Magnusson 2004, Tsakovska 2017)
license: free (publication-derived data)
---
# Flynn 1990 dataset — the foundational Kp compendium

## Description
Gerald L. Flynn's 1990 book chapter compiled in-vitro human skin permeability coefficient (Kp) data for ~94 compounds from the prior 30 years of literature, mostly from Franz-cell experiments dosed from aqueous vehicle. This was the dataset on which the Potts-Guy equation was fit (1992) and which has, in slightly modified form, served as the de-facto skin-Kp benchmark ever since.

## Coverage
- ~94 compounds (later papers report 90-97 depending on inclusion criteria).
- Predominantly small organic molecules: alcohols, phenols, steroids, simple aromatic compounds, parabens.
- Kp values from in-vitro permeation through human cadaver skin (variable preparation: full-thickness, dermatomed, epidermal).
- Aqueous donor vehicle (mostly).
- Kp range spans roughly 4-5 log units (1e-7 to 1e-2 cm/h).
- No vehicle effect data, no finite-dose data.

## Access
Not hosted as a single file. The compound list and Kp values appear:
- Original: Flynn GL, in *Principles of Route-to-Route Extrapolation for Risk Assessment* (Gerrity TR, Henry CJ, eds.), Elsevier, 1990, pp. 93-127.
- Reproduced in: Potts RO, Guy RH. Pharm Res 1992; 9: 663 (table); Magnusson BM et al. J Invest Dermatol 2004; 122: 993 (with additions); Vecchia BE, Bunge AL. In *Transdermal Drug Delivery* 2nd ed. (Guy RH, Hadgraft J, eds.) 2003 (curated table). Also as supplementary material in Tsakovska et al. Toxicology 2017.

## Citation
Flynn GL. "Physicochemical determinants of skin absorption." In: Gerrity TR, Henry CJ, eds. *Principles of Route-to-Route Extrapolation for Risk Assessment*. Elsevier; 1990: 93-127.

## Notes
- Heavy heterogeneity in skin source (different cadaver donors, preparations, ages).
- Some Kp values have since been re-measured with very different results — the dataset has known noise floor of ~0.5 log units.
- Vecchia & Bunge curated a "cleaned" subset; many ML papers use this rather than Flynn's full list.
- Despite known limitations, this remains the dataset every new skin-Kp QSAR is benchmarked against.
