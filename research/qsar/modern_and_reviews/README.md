# Modern ML, Reviews, and Databases for Skin Permeation QSAR

This subdirectory complements the classical Kp / D / K regression notes
maintained by sister agents. It covers three areas that the classical
notes do not:

1. **Modern statistical / ML approaches** to skin Kp prediction
   (post-2015), including ensembles, deep learning, GNNs, and
   transformers.
2. **Major review articles** that survey the skin QSAR field from
   multiple angles — mechanistic, statistical, regulatory.
3. **Public databases, toolboxes, and web servers** that contain skin
   permeation data or pre-built Kp / K predictors the user could call.

## Layout

```
modern_and_reviews/
├── README.md                          (this file)
├── ml_models/                         (one .md per model)
│   ├── ashrafi_gp_2010s.md
│   ├── baba_ann_2015.md
│   ├── baba_svm_2017.md
│   ├── chen_gnn_2022.md
│   ├── lim_deeppk_2022.md
│   ├── mansoor_qsar_lab.md
│   ├── modi_lstm_2017.md
│   ├── pawar_riviere_2018.md
│   ├── sun_rf_gbm_2018.md
│   ├── tsakovska_consensus_2017.md
│   └── wu_transformer_2023.md
├── reviews/
│   ├── anissimov_2018.md
│   ├── bouwman_geraets_2024.md
│   ├── chen_2013_2015.md
│   ├── kretsos_kasting_2007.md
│   ├── mathes_brandner_lehr_2014.md
│   ├── mitragotri_2011.md
│   ├── oecd_who_guidance.md
│   └── tsakovska_2017_atla.md
└── databases/
    ├── admetlab.md
    ├── biovia_dermal.md
    ├── chemaxon.md
    ├── comptox_dashboard.md
    ├── cosmetics_europe_lrss.md
    ├── cosmotherm_skin.md
    ├── edetox.md
    ├── epi_suite.md
    ├── flynn_1990.md
    ├── heeds_dermal.md
    ├── magnusson_2004.md
    ├── oecd_qsar_toolbox.md
    ├── opera.md
    ├── qikprop.md
    ├── r_python_packages.md
    └── simulations_plus_dermsim.md
```

## Top-level synthesis

### Modern ML for skin Kp — what's the state of play?

Three observations that emerge from the post-2015 ML literature:

1. **Performance is data-limited, not model-limited.** Across ANN
   (Baba 2015), SVM/RF/GBM (Baba 2017, Sun-Moss series, GP variants),
   GNN (Chen-Lian), LSTM (Modi-style), and transformer (Wu-style)
   approaches, the cross-validated R^2 on the canonical Flynn /
   Magnusson sets clusters in 0.78-0.86 for almost all well-tuned
   methods. The variance between architectures is smaller than the
   variance between data-curation choices on the same architecture.
   The bottleneck is the noise floor of the underlying in-vitro Kp
   measurements (~0.3-0.5 log units across labs), not the
   sophistication of the regressor.
2. **Consensus modelling is the most reliable improvement.**
   Tsakovska's 2017 cosmetics consensus and the Sun-Moss group's
   nested-CV ensembling consistently outperform any single-model
   approach. This is not surprising — averaging across heterogeneous
   noisy regressors is the textbook way to improve generalization on
   small noisy datasets.
3. **Vehicle / mixture handling is the open frontier.** Most published
   QSARs assume infinite-dose aqueous donor — wildly unrepresentative
   of cosmetic / topical drug exposures. The Riviere lab's LSER /
   Abraham-descriptor approach and Baba's vehicle-aware SVM are early
   attempts. The Cosmetics Europe LRSS data generation and recent
   industry-academic collaborations are pushing this forward, but
   there's no widely accepted vehicle-aware QSAR yet.

### Reviews — where to start

For a single-source orientation:
- **Tsakovska et al. Toxicology 2017** — best comprehensive review of
  skin Kp QSARs; covers all major models, datasets, and the EU
  regulatory context. Start here.
- **Mitragotri et al. Int J Pharm 2011** — best community-consensus
  review of the underlying mathematical models (steady-state QSAR,
  finite-element diffusion, brick-and-mortar microscopic). Useful
  context for what a tool like `skindiff` does.
- **Anissimov et al. Adv Drug Deliv Rev 2013** — best review of the
  analytical and inverse-fitting theory; relevant for understanding
  what `skindiff`'s `skin_fit()` is doing under the hood.

For regulatory context:
- **OECD GD 156 (2011)** + **WHO EHC 235 (2006)** — define what
  acceptable dermal absorption evidence looks like for risk assessment.

### Databases / toolboxes — what to actually use

For a practical open-source pipeline today:

1. **OPERA** + **CompTox Dashboard** — free, EPA-maintained, validated
   property predictor with skin Kp endpoint and domain-of-applicability
   flag. Best free starting point.
2. **Flynn 1990** + **Magnusson 2004** datasets — the standard QSAR
   training / validation reference data. Available as supplementary
   material in many published papers.
3. **`skindiff`** (this package) — open-source 1D Crank-Nicolson
   simulator, parameter fitting, and unit-safe R interface. The
   diffusion-simulator slot in an open pipeline.
4. **OECD QSAR Toolbox** — for read-across workflows and access to
   curated public dermal data; less useful as a Kp calculator.
5. **ADMETlab 3.0** — free web server for quick Kp estimates among 90+
   ADMET endpoints.

Commercial alternatives (BIOVIA / COSMOlogic, Schrödinger QikProp,
Simulations Plus DermSim, ChemAxon) have specific strengths but lock
in a closed-source pipeline.

### What's still missing

- **A modern, large, well-controlled, open-access skin Kp dataset.**
  The field has been working on the same Flynn / Magnusson compounds
  for 30 years. Cosmetics Europe's LRSS work and EPA / Hartung CITER
  efforts may eventually fill this gap.
- **A vehicle-aware QSAR with broad domain of applicability.** Riviere
  lab and Baba have made progress; nothing has reached widespread
  adoption.
- **Open-source mechanistic K prediction.** Mechanistic K (COSMO-RS
  family) is currently the province of commercial tools. An
  open-source alternative (ML on COSMO-RS-generated training data, or
  end-to-end DFT-light methods) is an unmet need.
- **Reproducible, version-controlled QSAR distribution.** Most
  published QSARs are bare equations in PDFs. Re-implementing each
  against curated descriptors is a per-paper effort.

## File-naming conventions used here

- ML models: `<lead-author>_<algorithm-or-name>_<year>.md`.
- Reviews: `<lead-author>_<year>.md`.
- Databases: `<short-name>.md`.

Each file has YAML front-matter for machine-readable metadata
(name, year, predicts/scope/type, performance/coverage, license/access)
followed by structured prose.

## Caveats on coverage

This research note was compiled from the agent's training knowledge
plus citation patterns common in the skin Kp QSAR literature. Several
specific citations should be verified before use in formal work:

- **Bouwman / Geraets 2024** — exact citation not unambiguously
  identified.
- **Lim et al. DeepPK** — likely refers to one of several recent
  Korean-group multi-task ADMET deep-learning papers; verify which.
- **Modi LSTM** — specific Modi reference for LSTM-on-SMILES skin
  permeation not pinned down; broader Modi-Unilever papers exist.
- **Mansoor lab** — lab identity not confirmed; multiple "Mansoor"
  authors in the dermal literature.
- **HEEDS dermal absorption** — most likely a workflow involving
  Dassault / SIMULIA tools driven by HEEDS, not a dermal-specific
  product.
- **Baba 2017 SVM** — DOI / citation pattern inferred from the series;
  verify exact issue.

Where uncertainty exists, the individual files flag it explicitly.
