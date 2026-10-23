# Test Ω × excavated depth — proposed rule

**Status when written: design awaiting approval; no job submitted.**
This document precedes any new sensitivity calculation. The values of the Ω test
are already known references, used only to lock
the instrument. The result and verdict fields of the new test remain empty.

## 1. Question

Does the effect of the burial velocity ω on the annual total of remineralised
carbon depend on the excavated thickness p_d, and by how much, with a fixed
hopper balance and then with a propagated hopper balance?

## 2. Reading rule written beforehand

### Design and quantities

Grid kept: ω = {0.1; 0.42; 1; 5} cm yr⁻¹ × p_d = {0.3; 0.5; 1} m,
for smooth and flag separately. The velocities
are the values of the YAML in service. The 0.3 m case is
motivated by the passage of Porz et al. (2026) read at the source; 0.5 m is an intermediate
design point, not an observation. No scenario is declared more realistic.

Two tests, analysed separately, each use the whole grid:

- **A — `frac_fixe` (main test)**: in-place fraction of `_r6` constant
  in each cell. SAR stays fixed; SVR becomes SAR × p_d; the depth
  shares and W_eff are recomputed. A isolates the depth effect from the
  change in hopper balance. This counterfactual implies that V_retire
  would vary in proportion to the swept volume if it were given
  a volumetric interpretation; it is not a new observed fleet.
- **B — `frac_recalculee` (propagation test)**: same changes, with
  V_retire and the cell area fixed to the table in service,
  V_balaye = SAR × p_d × cell_km2 × 10⁶,
  frac = max(V_balaye − V_retire, 0) / V_balaye. This test includes the change
  in the weights of the cells in the total when some are brought to zero.
  It does not prove that the same hopper removal is achievable at any depth.

Cells missing from the hopper table keep frac = 1, as in the merge;
their number is reported. A missing value in a row
present in the table is, however, a validation failure. grid_id are
never deleted: a cell that became inactive stays with a zero contribution.

For each test and label, T(ω,p) = Σ C_ri and
**Rω(p) = T(ω,p) / T(0.42,p)**. The denominator always has the same thickness,
the same hopper treatment and the same label. Totals are reported,
but their variation with p is not enough to diagnose an interaction.

Main contrast: ω = 5 and p = 0.3 m against p = 1 m.
All the planned contrasts of the grid are also printed:

- Iω(p) = Rω(p) / Rω(1) − 1: relative change of the multiplicative factor;
- ΔRω(p) = Rω(p) − Rω(1): difference of the ratios (×100 for percentage points);
- Eω(p) = Rω(p) − 1: effect relative to the reference ω;
- Qω(p) = Eω(p) / Eω(1): share of the excess retained, distinct from the ratio of the R;
- Bω(p) = Rω,B(p) / Rω,A(p) − 1: effect of hopper propagation on the ratio.

Q is undefined when Eω(1) is zero within the absolute numerical tolerance.
If T(0.42,p) = 0, the ratio concerned is **NULL**, never 0 or 1 by default.
Valid results at the other depths are kept, but the overall verdict
for this test and this label is incomplete. No confidence
interval or statistical test: this is a deterministic contrast.

### Validity lock — before any calculation at p_d other than 1 m

The job first checks the full fingerprints, the schema of the tables, the
unique and non-missing grid_id, identical fi/cri universes, finite values
and their domain. It checks that p_d is indeed 1 m in the active
cells, SVR = SAR × 1 m, C0i = C0i_mean_carbon_stock and d_i = 1.
The inactive exception of the merge (SAR=0, p_d and SVR missing)
is kept with its NA mask, before the NA→0 of step6; a missing SVR
in an active cell makes the lock fail.

The reconstruction of frac is checked on the whole table; SAR and SVR are
compared with `_r6` on the common cells. An old, geometrically
incompatible table **stops the whole test**: its V_balaye is not silently
corrected to make the lock pass. Cells of the table outside `_r6` do not
take part in the total; their identifiers and values stay in the frozen
input. `_r6` cells without a table follow the explicit fallback described above.

The instrument must recover, in both tests A/B and for both labels:

1. The W_eff, f_i before frac, f_i after capping and C_ri cell by cell
   of `_r6`, then the W_eff, f_i and C_ri of the four existing Ω reruns at p_d = 1 m.
   The covariates, stocks and sets of cells of the reruns are compared with
   `_r6`, to prevent agreement of totals from hiding a change of input.
2. The unrounded totals and ratios of the Ω Parquet files, to **10⁻⁶ relative**.
3. The numbers printed below to **0.000051 absolute**, in the unit of the
   number compared (dimensionless ratio, or tonnes of C for the totals).

| ω (cm a⁻¹) | smooth | flag |
|---|---:|---:|
| 0,1 | 0,8958 | 0,8958 |
| 0,42 | 1 | 1 |
| 1 | 1,1872 | 1,1870 |
| 5 | 1,9340 | 1,9335 |

Expected central totals: smooth 4 038.1542 t C; flag 5 654.3864 t C.
Source: reference outputs of the uncertainty tests (`read.R`); rounded values, never
substituted for the unrounded totals to compute the new ratios.

Per-cell tolerance: |calculation − reference| ≤ 10⁻¹² + 10⁻⁹ × |reference|,
with identical NA masks. **Design choice, no source**:
the absolute floor avoids a relative requirement undefined at zero; the relative term
is close to V0 (lock of the uncertainty tests), while checking the cells individually.
The 10⁻⁶ threshold on totals takes up V0; it is much stricter than the
interpretation threshold. The 0.000051 tolerance is a **design choice, no
source**: half of the last displayed digit is 0.00005,
with a margin of 0.000001 for numerical rounding. It does not exempt from the
check against the unrounded Parquet files.

**A single failure ⇒ invalid/incomplete instrument: stop, no calculation at p≠1 m,
no conclusion on the interaction, no automatic relaxation of the threshold.**
The diagnostics of the lock are saved, including on failure. A later
error also prevents the COMPLETE marker and any substantive interpretation.

### Substantive verdict, only after all locks

Hypothesis tested: the Ω factor obtained at 1 m is transportable to the depths
tested with less than 5 % relative difference.

- **Main contrast maintained** if |I5(0.3)| < 0.05; **failure of
  transportability** otherwise. This failure is a result of the model, not a failure
  of the instrument: the factor is conditional on p_d. Report the sign,
  R5(0.3), R5(1), ΔR and Q; do not confuse a reduction of R with a reduction of E.
- **Transportability maintained over the grid** if max |Iω(p)| < 0.05 for
  ω other than 0.42 and p other than 1 m; **conditional** otherwise.
  Conclusion common to both labels only if both satisfy
  the criterion; otherwise present both. This is a finite grid, not a
  demonstration for all intermediate depths/velocities.
- The effect at ω = 5 becomes **secondary at 0.3 m** if |E5(0.3)| < 0.05.
  This is the operational definition here of a "collapse" to a small residual
  effect. If E ≤ −0.05, the sign is reversed; if E ≥ 0.05, the effect remains
  positive and not secondary, even if Q indicates a strong attenuation.
- The influence of hopper propagation is **material** if |Bω(p)| ≥ 0.05.
  An effect present only in B cannot be attributed to the profile alone.
  The A/B contrast is deterministic and is not a causal additive
  decomposition when caps exist.

**Provenance of the 5 %**: threshold of the earlier uncertainty tests (test G, then test Ω).
Its extension to I and B is a **design choice, no source**,
justified by a convention common to Table 2; it is neither
a universal threshold nor a value taken from Atwood, Middelburg or Porz.
Equality at the threshold belongs to the material/non-transportable case.

## 3. Raw result

<!-- To be filled in after approval and execution. -->

## 4. Verdict against the rule

<!-- To be filled in after approval and execution. -->

## 5. Provenance and instrument retained

**Recalculation from the `_r6` outputs, run on GRIT by sbatch.**
This choice avoids the raster reprojection and the heavy GeoTIFF outputs of step6,
which depend on neither ω nor p_d. Above all, it makes it possible to observe the values
before each cap. It is not a simple rescaling of C_ri.

The code parses the merge whose fingerprint is guarded and evaluates only the
definitions of `weight_depth()` and `weight_eff()`, in an environment
containing the parameters of the YAML. It never sources the body of the pipeline.
It reconstructs, in the order of the merge: SVR, w1/w2/w3/H3, W_eff, p_l_corr,
annual kinetics, f_i before frac, multiplication by frac then cap [0,1].
Then it applies the NA→0 of step6 and min(C0i × f_i × d_i, C0i).

Held fixed: SAR, p_l, k_used, fresh_fact, D_b, the fractions and constants
of the pools, p_r, d_i, coordinates, spatial universe and C0i. These parameters do not
depend on ω/p_d in the files checked. The source columns are
checked against the existing reruns; no division by a capped
or zero f_i is used. No new AIS classification or hopper attribution.

**A full-rerun alternative is possible without modifying the pipeline**: write
copies of the tiles to a new output, replace only sum_dw_pd with
p_d × sum_dw, keep sum_dw (so the owner/max sum_dw selection is unchanged), then
pass these tiles through MERGE_IN_DIR and a frac table specific to the test through FI_FRAC_FILE.
The launcher logs and the temporary files of step6 would also have to be
redirected, since they otherwise write outside the authorised folder. This route is
**not implemented**, because the locked recalculation is enough for the question asked.
Modifying only W_eff, without SVR and without the frac options, would be insufficient.

Fingerprints of the versions executed (before the translation of their comments):

| Object | Full md5 fingerprint | Use |
|---|---|---|
| pipeline/step5_merge_tiles_corr.R | 85e94d4758aeda3e76f74e8e290aa2c6 | depth functions; SAR/SVR/depth; kinetics/frac/cap |
| pipeline/step6_calculate_cri_corr.R | d258c7f009038dda1410ff549d61ce0b | NA, raster, C0i, d_i and C_ri |
| fi_parameters_with_freshness.yaml | f298996c85c16bc56e47e93e7bf27c03 | central scenario and ω |
| pipeline/merge_rerun.sbatch | f37758f4de423486d582f2b6895dbd9f | file actually executed |
| pipeline/step6_pr.sbatch | e09e9e4c0c3d90d63b2b75bc15a0e7dc | input and script executed |

`inputs.tsv` gives the full path, the bytes and the md5 of **each** source
file, fi_grid, cri_final and frac table. A Parquet file has no
"text line": observations are identified by grid_id and the
columns are described in `compute.R`. The Ω paths are frozen after requiring
a single Parquet file per case/label/type; never "take the most recent".
`bundle.md5` freezes the rule, the scripts and the manifest; `provenance.tsv`, the
fingerprints of the results and the Slurm job accompany each new number.

### Explicit position on the five pitfalls

1. **Double entry of p_d**: SVR and the normalised weights are changed together.
   R uses a denominator of the same depth. The linear volume effect
   cancels in this ratio only without redistribution or capping;
   it is therefore not assumed that normalisation removes all effects of p_d.
2. **Hopper**: A and B are distinct. B recomputes V_balaye, V_sur_place and frac
   in memory; V_retire stays fixed. n_frac_zero, n_new_frac_zero,
   the fallbacks to frac=1 and the C0i stock carried by cells at frac=0 are reported.
3. **Caps**: count of raw f_i ≥1, of f_i strictly >1, of f_i ≤0,
   of NA replaced by zero, then of C_ri at the stock and strictly clipped;
   shares of ΣC_ri attributable to capped cells. Total count and active
   count for f_i. Cells with C0i=0 are excluded from the carbon cap so as not
   to count the trivial equality 0=0 as saturation. No arbitrary threshold
   "close to the cap". With d_i=1 and f_i≤1, the second cap is
   algebraically redundant, but its diagnostic is printed anyway.
   A diagnostic total without caps and its ratio are also saved.
4. **C0i stock**: Atwood et al. 2020, *Materials and Methods / C Mapping*, specifies
   a stock of the first metre; Fig. 6 describes the same horizon. Step6 extracts
   `mean_carbon_stock` and renames it C0i, without p_d. This
   reference stock is kept in A and B: the perturbed depth already
   enters through SVR and W_eff. Multiplying it again by p_d would change the
   convention of the model and double the geometric factor under the assumption
   of a uniform stock. **Assumption**, not a demonstration of a real profile:
   the link between this integrated stock and the reactive normalisation over 0–5 cm
   remains that of the pipeline. The raster does not provide the observed stock over
   the top 30 cm alone; the test does not estimate one. The download chain
   from the local raster to the published file was not checked here.
5. **Limits**: the permanent profile does not represent the succession of
   deposits and maintenance dredging; no value of ω/p_d is validated
   as the right one. No depletion correction, no OCIM, no airborne
   fraction, no probabilistic uncertainty. The test demonstrates neither a
   real underestimate nor a lower bound of the emissions.

A precision from the code: ω also enters the two roots of the
0–10 cm profile and its 0–5 cm normalisation (`weight_depth`).
The result must not be attributed exclusively to the layer below 10 cm.
The test does not decompose the contributions by layer and by pool.

Primary sources consulted:
[Atwood 2020, C Mapping and C Model Error and Uncertainty](https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2020.00165/full)
and [Porz 2026, Methods passage on the initialisation of the model](https://www.nature.com/articles/s41467-025-68105-5).
Porz averages the reactivity classes over the top 30 cm, in relation to
the dredging furrows; this motivates a scenario, not a universal depth
for the ten vessels. The primary references cited by Porz at this point
were not reread. Middelburg was not reread here: values attributed
via the YAML, distinct from their verified presence.

## 6. What each outcome authorises in the paper

| Outcome | Reading | Proposed authorisation, always subject to approval of the text |
|---|---|---|
| A lock fails / incomplete calculation | invalid or incomplete instrument | Nothing new in Table 2; keep and explain the failure in the project log |
| Zero denominator at one depth | NULL contrast at that depth | Qualitative limit of the hopper balance; no invented value of R; other valid contrasts explicitly bounded |
| A transportable, B transportable | factor stable under the convention retained over the grid | Numbered line giving the observed envelope and its domain, without claiming universal independence |
| A not transportable | dependence on p_d even with fixed frac | Numbered line conditional on p_d, giving at least R5(0.3) and R5(1), and a limit of the model |
| A transportable, B not transportable | change due to hopper propagation and its weights | Present A and B separately; do not attribute to the depth of the profile alone |
| Secondary effect at 0.3 m | the excess at high ω becomes <5 % | Quantify the attenuation and withdraw any generalisation of the ×1.93 to this depth |
| Two labels or tests disagree | narrower domain of validity | Separate values; common qualitative wording only if justified |

A negative transportability result remains publishable in the Assessment.
A "maintained" verdict can also be reported with its maximum gap; it
is not an empirical validation. No automatic change to the manuscript.

## Execution after approval only

`prepare.py` freezes the paths and fingerprints, without R and without touching the calculation
outputs. `launch.sh` is DRY=1 by default. For a real run, it requires
DRY=0 **and** APPROVED_RULE_MD5 equal to the md5 of this version of RULE.md after
approval. It atomically reserves a submission, refuses any non-empty
output, checks the fingerprints before submission and at job start, and
never stays connected to wait. Any R execution goes through sbatch,
partition grit_nodes and a complete exclusion list.

The job writes only under `$DFP_SCRATCH/output_V6_incert/omega_pd/`
(temporary files included), the submission logs under
`$DFP_REPO/assessment/omega_pd/`. The commands displayed are filtered by
`grep -viE 'token|secret|key|password|client_id|auth'`.

`read.sh` reads the saved results, checks their md5 and refuses
interpretation without LOCK_PASSED and COMPLETE. The job produces
`LOG_ENTRY.md` with the six fields (or `FAILURE.md` on failure):
**copy it into the project log as soon as it is read,
before any interpretation or other test**. The cluster cannot
write into the log itself, which is on the local computer.
