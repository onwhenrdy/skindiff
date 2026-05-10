---
name: tsakovska_2017
year: 2017
scope: comprehensive review and benchmarking of skin permeability QSARs from 1989 to 2017, with focus on cosmetic ingredients
---
# Tsakovska et al. 2017 — QSAR for skin permeability under the EU cosmetics testing ban

## Citation
Tsakovska I, Pajeva I, Al Sharif M, Alov P, Fioravanzo E, Kovarich S, Worth AP, Richarz AN, Yang C, Mostrag-Szlichtyng A, Cronin MTD. "Quantitative structure-skin permeability relationships." Toxicology 2017; 387: 27-42. DOI: 10.1016/j.tox.2017.06.008.

A companion / earlier publication in ATLA (Alternatives to Laboratory Animals) covers similar ground from a regulatory perspective.

## Coverage
The Toxicology 2017 review is the most comprehensive single reference for the skin Kp QSAR field. It covers:
- All major Kp regression equations from Potts-Guy 1992 onward (Magnusson 2004, Wilschut 1995, Patel 2002, Moss-Cronin 2002, ten Berge 2009, Chen-Lian 2007, Lian 2008, Kupczewska-Dobecka 2010, etc.).
- The two foundational datasets (Flynn 1990, Magnusson 2004) and their derivatives.
- Statistical, neural-network, fragment-based, and mechanistic approaches.
- Domain-of-applicability considerations.
- The European regulatory context: the 2013 cosmetic-animal-testing ban, REACH, and OECD TG 428 in-vitro permeation testing.

Each model is summarized with its descriptors, training set, performance metrics, and a critical assessment.

## Key takeaways
1. The Potts-Guy equation, despite being thirty years old, remains a respectable baseline; many "improved" QSARs gain only marginal accuracy on the original Flynn data.
2. Differences between published models are often smaller than the experimental noise in the underlying Kp data (skin source, vehicle, lab-to-lab variation).
3. No single QSAR is universally best; consensus modelling outperforms any single regression.
4. The field's data ceiling is the key limit: any model trained on Flynn-derived data inherits Flynn's heterogeneity. Better data (controlled experimental conditions, single skin source) would let a much simpler model do much better.
5. Extrapolation to cosmetic ingredients (often outside the training-data physicochemical space) is unreliable — high logP, surface-active, or polymeric ingredients lie outside the QSAR domain of applicability.
6. Combining QSAR with mechanistic diffusion simulation (i.e. the kind of pipeline `skindiff` enables) is more useful for risk assessment than either alone.

## References
- Potts RO, Guy RH. Pharm Res 1992; 9: 663.
- Magnusson BM et al. J Invest Dermatol 2004; 122: 993.
- Wilschut A et al. Chemosphere 1995; 30: 1275.
- Cronin MTD, Dearden JC, Moss GP, Murray-Dickson G. Eur J Pharm Sci 1999; 7: 325.
- Many further refs in Tsakovska 2017 itself (~150 cited papers).
