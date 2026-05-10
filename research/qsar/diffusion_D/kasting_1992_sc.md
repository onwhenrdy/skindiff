---
name: kasting_1992_sc
predicts: D
layer: SC
units: cm^2/s (D/h^2 reported as 1/s)
descriptors: MW (Da)
domain: small molecules through human stratum corneum; finite/infinite dose
status: classical
---

# Kasting–Smith–Anderson 1992 — MV-scaled D in stratum corneum

## Citation
Kasting GB, Smith RL, Anderson BD. "Prodrugs for dermal delivery: solubility,
molecular size, and functional group effects." In: Sloan KB, ed. *Prodrugs:
Topical and Ocular Drug Delivery*. New York: Marcel Dekker; 1992:117–161.
A widely-cited follow-on extension is Kasting GB. J Pharm Sci. 2001;
90(2):202–212 (vanillylnonanamide finite-dose) — DOI 10.1002/1520-6017.

## Equation
Kasting fit log P_SC and log(D_SC/h^2) of human SC against molecular weight
in linear-in-MW form:

```
log10(D_SC / h_SC^2) [1/s] = a - b * MW         (MW in Da)
```

with a, b > 0, slope b ~ 0.012–0.018 1/Da depending on dataset cut. Holding
h_SC fixed (e.g. 25 um) gives D_SC ~ 10^(-b*MW) cm^2/s, i.e.
**exponential in MW**. For MW = 100, D_SC ~ 1e-9; for MW = 500, D_SC ~ 1e-15
cm^2/s — the same blow-up the lipid-bilayer SPT exponential gives, in a
different parameterisation.

The Kasting 1992 / 2001 line was later refined by Wang–Kasting–Nitsche
(see `wang_kasting_2007_dlip`), where the Kasting MW dependence was
mapped onto the trans-bilayer hopping coefficient log10(k_trans) =
A − B * MW^(1/3) — converting the MW-linear log form into a power-law
in MW^(1/3).

## Training set
Subset of Flynn's compiled human stratum corneum permeability database
plus Kasting's own in-house in-vitro data. ~50–80 compounds, MW 18–400,
log K_ow −2 to 5. The 2001 paper revisits the regression and adds
finite-dose vanillylnonanamide kinetics.

## Reported performance
log10 D root-mean-squared residuals ~0.4–0.6 (i.e. ~3x in D); explains
~70% variance in the training fit. Standard deviation in slope b across
sub-datasets is the dominant uncertainty.

## Validity / limitations
Lumps SC into a single homogeneous layer with effective D — no separation
into lipid vs. corneocyte phases. Best suited as a starting estimate for
D_SC when only MW is known. Breaks down for very small (MW < 30) or very
large (MW > 600) molecules and for ionised species.

## Notes
Functionally interchangeable with Potts–Guy (which is a Kp QSAR, not a D
QSAR) — Kp = D K / h_SC, so given Kp from Potts–Guy and an estimate of K,
one obtains the same D ballpark this model gives directly.

## References
- Kasting GB, Smith RL, Anderson BD, in *Prodrugs* ed. Sloan KB. Marcel
  Dekker, 1992 (book chapter).
- Kasting GB. J Pharm Sci. 2001;90(2):202–212. DOI
  10.1002/1520-6017(200102)90:2<202::AID-JPS9>3.0.CO;2-K
- Kasting GB, Miller MA. J Pharm Sci. 2006;95(2):268–280. DOI
  10.1002/jps.20497 (volatile finite-dose extension).
