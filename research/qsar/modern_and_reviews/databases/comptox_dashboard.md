---
name: comptox_dashboard
type: web server + database
url: https://comptox.epa.gov/dashboard/
license: free (U.S. government public data)
---
# EPA CompTox Chemicals Dashboard

## Description
The CompTox Chemicals Dashboard is the U.S. EPA's public web portal aggregating chemical property and toxicity data for ~1.2 million substances. Each chemical's page shows experimental and predicted property values from many sources, including OPERA model predictions (which include a skin Kp endpoint), in-vivo and in-vitro toxicity data from EPA's ToxCast / Tox21 programs, exposure model outputs, and links to safety data sheets.

For dermal absorption work specifically: the Dashboard's value is as a one-stop lookup of (a) any measured Kp data EPA has aggregated, (b) the OPERA Kp prediction with domain-of-applicability flag, and (c) related properties (logP, solubility, MW) that feed downstream Kp QSARs.

## Coverage
- ~1.2 million chemical substances.
- Skin Kp data: only for compounds where EPA has aggregated published values (small subset — overlaps with Flynn/Magnusson and a few hundred others).
- Predicted skin Kp from OPERA: available for any chemical with parsable structure (SMILES).
- Connected to PubChem and other public chemistry databases via cross-references.

## Access
- Free public web at https://comptox.epa.gov/dashboard/
- Search by name, CAS, SMILES, InChI.
- Batch search and bulk download for hundreds-to-thousands of chemicals via the Dashboard's batch interface.
- API: https://api-ccte.epa.gov/docs/ — REST API for programmatic access.
- All underlying data downloadable as CSV / SDF.

## Citation
Williams AJ, Grulke CM, Edwards J, McEachran AD, Mansouri K, Baker NC, Patlewicz G, Shah I, Wambaugh JF, Judson RS, Richard AM. "The CompTox Chemistry Dashboard: a community data resource for environmental chemistry." Journal of Cheminformatics 2017; 9: 61. DOI: 10.1186/s13321-017-0247-6.

## Notes
- Best free, government-maintained source for an integrated view of any chemical's properties including skin Kp.
- OPERA model integration means quick predicted Kp values with calibrated uncertainty / domain-of-applicability indicator for any new compound.
- Data quality is generally high (curated, version-controlled, linked to original sources). The Dashboard is the modern successor to EPA's older DSSTox portal.
- A natural front-end / property source for any skin permeation pipeline that needs to handle arbitrary user-supplied chemicals.
