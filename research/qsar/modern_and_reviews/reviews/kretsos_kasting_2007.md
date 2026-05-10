---
name: kretsos_kasting_2007
year: 2007
scope: review of D and K (diffusivity and partition) in skin compartments — vehicle, stratum corneum, viable epidermis, dermis
---
# Kretsos & Kasting 2007 — Diffusivity and partition coefficient review

## Citation
Kretsos K, Kasting GB. "A geometrical model of dermal capillary clearance." Mathematical Biosciences 2007; 208: 430-453. (The "review of D and K in skin" reference is a separate Kretsos-Kasting paper from the same era — likely Kretsos K, Miller MA, Zamora-Estrada G, Kasting GB. "Partitioning, diffusivity and clearance of skin permeants in mammalian dermis." Int J Pharm 2008; 346: 64-79.) Verify the precise citation; the Kretsos / Kasting collaboration produced several papers in 2005-2010 reviewing different aspects of dermal D and K.

## Coverage
The Kretsos-Kasting 2008 paper specifically reviews:
- Measured and modelled D and K values in stratum corneum, viable epidermis, and dermis.
- Capillary clearance from the dermis as a "sink" boundary condition (a key piece for in-vivo IVIVE).
- Compartment-specific partition coefficients (K_SC/water, K_VE/water, K_dermis/water) and how they relate to logP.
- Effective diffusivity in each layer; SC has D ~ 1e-10 to 1e-12 cm^2/s for typical drugs, four to six orders of magnitude lower than aqueous diffusivity.

## Key takeaways
1. Partition coefficient K_SC/water roughly tracks K_octanol/water for moderately lipophilic compounds but deviates strongly at extremes; QSARs based purely on logP lose accuracy outside the middle range.
2. Viable epidermis and dermis have D values much closer to aqueous (~ 1e-7 to 1e-8 cm^2/s) than SC; for most compounds these layers are not rate-limiting and a single-compartment SC model often suffices.
3. Dermal capillary clearance is the right physical boundary in vivo; treating the systemic compartment as a perfect sink is an approximation that overestimates flux for slowly-cleared compounds.
4. For the kind of multi-layer simulator `skindiff` implements, sensible D and K defaults are: D_SC ~ 1e-10 cm^2/s, K_SC/water tracking logP weakly; D_VE ~ 1e-7 cm^2/s with K ~1; D_dermis ~ 1e-7 cm^2/s with K ~1.

## References
- Kretsos K, Kasting GB. Math Biosci 2007; 208: 430.
- Kretsos K et al. Int J Pharm 2008; 346: 64.
- Anissimov YG, Watkinson A. "Modelling skin penetration." Skin Pharmacol Physiol 2013; 26: 314 — complementary review of D and K.
- Frasch HF, Barbero AM. "The transient dermal exposure: theory and experimental examples." Skin Pharmacol Physiol 2008; 21: 251.
