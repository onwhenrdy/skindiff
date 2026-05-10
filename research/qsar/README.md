# QSAR / structure-property models for skin permeation

Survey of published QSAR / structure-property models that predict the inputs
the `skindiff` simulator consumes. This directory is excluded from the
package build (`^research$` in `.Rbuildignore`) — it is reference material
for choosing which QSARs (if any) to wrap as helper functions later.

## What's here

```
research/qsar/
├── README.md             (this file — start here)
├── permeability_kp/      Kp regressions (lumped permeability)
├── diffusion_D/          D (diffusion coefficient) per skin layer
├── partition_K/          K (partition coefficient) per skin layer
└── modern_and_reviews/   ML models, review articles, databases, toolboxes
    ├── ml_models/
    ├── reviews/
    └── databases/
```

`skindiff` takes **D and K** as direct per-layer inputs, so `diffusion_D/`
and `partition_K/` are the folders most relevant to wrapping. Kp models
are useful as sanity-check predictions of the lumped quantity
`Kp ≈ D_SC · K_SC / h_SC` and as the largest single source of skin
permeation literature; `permeability_kp/` collects them.

Every per-model file follows the same template (YAML frontmatter +
Citation / Equation / Training set / Performance / Validity / Notes /
References). Each sub-folder has its own `README.md` with a comparison
table.

## File counts

| Folder | Models | What |
|---|---|---|
| [permeability_kp/](permeability_kp/) | 16 | Potts-Guy lineage; Flynn 1990 / Magnusson 2004 re-analyses |
| [diffusion_D/](diffusion_D/) | 17 | D in aqueous, SC, SC-lipid, corneocyte, VE, dermis |
| [partition_K/](partition_K/) | 13 | K vs water for each layer + logKow predictors overview |
| [modern_and_reviews/ml_models/](modern_and_reviews/ml_models/) | 11 | Post-2015 ML/DL approaches (RF, SVM, GNN, transformer) |
| [modern_and_reviews/reviews/](modern_and_reviews/reviews/) | 8 | Major review articles + OECD/WHO guidance |
| [modern_and_reviews/databases/](modern_and_reviews/databases/) | 16 | Datasets + toolboxes + web servers |
| **Total** | **81** | |

## Headline findings (read these first)

### The MW exponent in `D` disagrees by 5+ log decades

`D` in the stratum corneum is the most uncertain input to any skin
diffusion model. Published QSARs disagree by 5+ orders of magnitude on
`D_SC` for a typical drug (MW=200, logKow=2), driven entirely by what
exponent they put on MW:

| Form | Source | Effective MW scaling |
|---|---|---|
| `D ∝ MW^(-1/3)` | Stokes-Einstein | very weak |
| `D ∝ MW^(-0.6)` | Wilke-Chang, Hayduk-Laudie (aqueous) | mild |
| `D ∝ exp(-c·MW^(2/3))` | Mitragotri 2002 SPT (lipid) | steep, area-exponential |
| `D ∝ MW^(-2.43)` + floor | Wang-Kasting 2007 lateral lipid | moderate, with kinetic floor |
| `D ∝ 10^(-0.79·MW^(1/3))` | Wang-Kasting 2007 trans-bilayer | exponential of cube-root |
| `D ∝ 10^(-b·MW)`, b~0.014 | Kasting 1992/2001 | exponential in raw MW (very steep) |

See [diffusion_D/README.md](diffusion_D/README.md) for the cheat-sheet
table with prefactors and validity ranges.

### The Kow exponent in `K` is split by phase

K-vs-Kow exponent in `K_layer/water = (prefactor) · K_ow^β` is much
better-behaved, but separates cleanly by phase:

| Phase | β consensus | Sources |
|---|---|---|
| SC lipid / water | 0.67–0.81 (~0.7–0.75) | Anderson-Raykar, Frasch-Barbero, Mitragotri, Wang-Kasting, Hansen, COSMOmic |
| Corneocyte / keratin / water | 0.31–0.32 (very tight) | Anderson-Raykar, Hansen 2011, Wang-Kasting |
| Whole-tissue SC / water | 0.71–0.81 | Potts-Guy implied (0.71), Cleek-Bunge (0.74) |
| Viable epidermis / water | ~0 baseline + small lipid bump | Nitsche-Kasting 2013 |
| Dermis / water | ~0 baseline + small lipid bump | Kretsos-Kasting 2008 |

See [partition_K/README.md](partition_K/README.md) for the per-model
table and recommended defaults for `skindiff`.

### Performance is data-limited, not model-limited

Modern ML approaches (Baba SVR, Sun RF, Chen GNN, Wu transformer)
plateau around R² ≈ 0.75–0.91 on the same Flynn 1990 / Magnusson 2004
datasets that classical regressions used. The bottleneck is dataset
size (~150–300 reliable in-vitro measurements, with ~0.5 log unit
inter-lab variance on the same compound), not the regression algorithm.
Consensus modelling (averaging over multiple QSARs) is the most
reliable improvement; vehicle / mixture handling is the open frontier.
See [modern_and_reviews/README.md](modern_and_reviews/README.md) for
the synthesis.

## Suggested reading order

If you have an hour:

1. [permeability_kp/potts_guy_1992.md](permeability_kp/potts_guy_1992.md) — the foundational regression every later paper benchmarks against.
2. [diffusion_D/README.md](diffusion_D/README.md) — the MW-exponent cheat sheet.
3. [partition_K/README.md](partition_K/README.md) — the Kow-exponent table + recommended defaults.
4. [modern_and_reviews/reviews/mitragotri_2011.md](modern_and_reviews/reviews/mitragotri_2011.md) — the community-consensus review.
5. [modern_and_reviews/reviews/tsakovska_2017_atla.md](modern_and_reviews/reviews/tsakovska_2017_atla.md) — most cited single QSAR review.

If you have ten minutes: read the three sub-folder READMEs and
this file.

## Caveats / TODO

- Equations and coefficients were copied from primary papers where the
  PDFs were accessible, and from authoritative secondary reviews
  (Mitragotri 2011, Kupczewska-Dobecka 2010) where the originals were
  paywalled. Files explicitly flag which is which. Before using any
  coefficient in production, verify against the primary source.
- A few specific items have known low-confidence flags
  (Lien-Gao 1995 intercept, ten Berge 2009 SKINPERM coefficients,
  Mansoor lab identity, Bouwman/Geraets 2024 — see the relevant files
  for details).
- This is a static literature snapshot from the survey date. New
  papers appear regularly; OPERA and ADMETlab in particular are
  actively developed.
