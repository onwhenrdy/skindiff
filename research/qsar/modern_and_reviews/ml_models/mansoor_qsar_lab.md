---
name: mansoor_lab
year: ~2015-2022
predicts: dermal absorption / Kp; some published work on transdermal delivery enhancers
algorithm: classical regression, MLR, SVM, depending on study
descriptors: physicochemical + topological
training_set: varies; small focused sets
performance: study-dependent
---
# Mansoor lab and other smaller groups

## Citation
The user-provided list mentions "Mansoor lab" — this is most likely a reference to one of several research groups working on transdermal delivery / penetration enhancement and Kp QSARs. Without more context (institution, country) it's hard to identify the specific group. Candidates include:
- Mansoor TA et al. (transdermal patches and chemical permeation enhancers, multiple papers).
- A. Mansoor in penetration enhancement studies.

The honest answer is that this entry could not be unambiguously verified. The user should clarify the lab identity if a specific reference is needed.

## Approach
Generic pattern shared by smaller QSAR groups working in this space: small in-house Kp dataset (often combined with literature data), classical descriptors from MOE / Dragon / RDKit, MLR or SVM, internal cross-validation. Few of these models become widely used because the larger benchmarking studies (Tsakovska 2017, Moss / Sun series) demonstrate that the underlying data quality limits any model's accuracy.

## Performance
Study-specific. Generally R^2 in the 0.6-0.8 range on small test sets, with the usual caveat that small-test-set R^2 is optimistic.

## Availability
Mostly papers only; little public code or data release.

## Limitations
The small-group QSAR pattern in this field has produced many redundant papers on the same Flynn / Magnusson data with marginal improvements. The Tsakovska 2017 review is the right entry point for a balanced view of the field.

## References
- See Tsakovska et al. Toxicology 2017; 387: 27-42 for a survey of which QSAR groups have produced skin Kp models.
- For penetration-enhancer QSARs (a different but related question): Karande P, Mitragotri S. Adv Drug Deliv Rev 2009; 61: 1147.
