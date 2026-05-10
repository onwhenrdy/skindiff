---
name: opera
type: software package + web server (open source)
url: https://github.com/NIEHS/OPERA ; https://comptox.epa.gov/dashboard/ (web access via CompTox)
license: open source (MIT)
---
# OPERA — OPEn structure-activity Relationship App

## Description
OPERA is an open-source suite of QSAR models maintained by the U.S. EPA (formerly NIEHS NICEATM); it predicts a wide range of physicochemical and ADMET properties from chemical structure. The current version includes 17+ endpoints including some dermal-relevant ones (logP, water solubility, vapor pressure, melting point, BCF). Skin permeability is a more recent addition in v2.x — uses a kNN regression on Flynn-derived data with curated descriptors.

Also delivers domain-of-applicability flags (one of OPERA's signature features) so users see whether a query compound is in or out of the training distribution.

## Coverage
- ~17 endpoints in current release, including a skin Kp model.
- Training data drawn from public sources (Flynn / Magnusson for skin Kp, EPI Suite databases for general properties).
- ~1500-12000 compounds per endpoint depending on which.
- Skin Kp model: smaller training set (~300 compounds, Flynn + Magnusson curated) than the general property models.

## Access
- Open source on GitHub: https://github.com/NIEHS/OPERA
- Pre-built Windows / Linux binaries.
- Integrated into the EPA CompTox Chemicals Dashboard (https://comptox.epa.gov/dashboard) — type any chemical, get OPERA predictions including Kp where applicable.
- Command-line and KNIME node interfaces for batch processing.

## Citation
Mansouri K, Grulke CM, Judson RS, Williams AJ. "OPERA models for predicting physicochemical properties and environmental fate endpoints." Journal of Cheminformatics 2018; 10: 10. DOI: 10.1186/s13321-018-0263-1.

Updated periodically; recent versions documented in subsequent J Cheminform / EHP papers from the EPA CompTox group.

## Notes
- Best free, openly licensed, validated multi-endpoint property prediction tool. Strongly recommended as the first stop for any property prediction including skin Kp screening.
- Skin Kp model is one of the smaller / less-validated endpoints in the suite — for serious Kp work, supplement with dedicated tools.
- CompTox Dashboard integration means a quick web lookup for any chemical (by CAS / SMILES / name) returns the OPERA Kp prediction in seconds.
- Domain-of-applicability flag is unusually rigorous - OPERA explicitly tells you when the query is outside the training distribution.
