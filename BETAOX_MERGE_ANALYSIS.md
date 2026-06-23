# Analysis — Merging β-Oxidation into the Collins et al. Model

*Scoping study: what it would take to port the fatty-acid β-oxidation module
from `mitochondria-metabolism` (Meyer et al. / van Eunen et al. 2013) into
`mitochondrial-carb-tca-oxphos` (Collins et al.).*

> This is an **analysis only** — no code has been changed. It identifies the
> coupling seams, the convention mismatches, the exact wiring points in
> `dXdT.m`, the decisions someone must make, and an effort estimate.

---

## 0. Verdict (TL;DR)

**Feasible, but it is a genuine model-integration project, not a copy-paste.**
The β-oxidation pathway is a *self-contained kinetic island* with only **four
narrow chemical seams** to the rest of metabolism (acetyl-CoA, NADH, CoA, and
FADH₂). That is the good news. The hard news is that the donor module and
Collins use **different unit systems, volume bases, moiety-pool conventions, and
thermodynamic philosophies**, so all four seams must be re-derived rather than
reused.

| | |
|---|---|
| **Core mechanism transplant** | Low risk — the 12 enzyme rate laws are pure functions of concentrations |
| **Unit / volume reconciliation** | **High risk** — donor is *per-mg-protein*; Collins is *per-L-water / per-L-mito* |
| **Moiety-pool wiring (CoA/NAD/FAD)** | Medium risk — Collins uses explicit state vars where Meyer uses conserved-pool algebra |
| **Thermodynamic consistency** | Decision required — accept a "kinetic island" or reformulate |
| **Rough effort** | ~1–2 focused weeks to wire + integrate; **several more weeks to re-tune & validate** |

---

## 1. What the β-ox module actually is

Source: `Cell_dxdT.m` lines ~1378–1790 (donor model). Lineage: van Eunen et al.,
*PLoS Comput Biol* 2013 (Bakker lab), exported from SBML.

- **47 state variables** (donor indices 63–109): for 7 chain lengths
  (C16→C14→C12→C10→C8→C6→C4), each tracks acyl-carnitine (cytosol + matrix),
  acyl-CoA, enoyl-CoA, 3-OH-acyl-CoA, 3-keto-acyl-CoA — plus a single matrix
  **FADH₂** pool.
- **12 enzymes** as pure MATLAB functions of concentration:
  `CPT1, CACT, CPT2, VLCAD, LCAD, MCAD, SCAD, CROT, MSCHAD, MCKAT(A/B), MTP`,
  plus a phenomenological `RES` FADH₂ sink.
- **Parameters:** `boxVmax(1:12)`; a global scale `SF = 1e6` (µM→M) dividing
  every Km; chain-length **specificity factors** `sf*` (in `s_F.mat`/`sf_out.mat`);
  held-constant pools `const_species_*` (cytosolic carnitine, CoA, **malonyl-CoA**
  = CPT1 inhibitor, matrix carnitine, total FAD).
- **Native conventions (critical):**
  - Concentrations in **M**, time in **s**.
  - Volume basis is **per mg mitochondrial protein**:
    `compartment_VMAT = 1.8e-6 L mito/mg`, `compartment_VCYT = 1.8e-5 L cyto/mg`.
  - Rate laws are **plain reversible Michaelis-Menten** — *no* H⁺/Mg²⁺/K⁺
    binding, *no* membrane-potential dependence, fixed `Keq`.

---

## 2. The four coupling seams (donor side)

In the donor's in-vivo model the β-ox block touches the shared core in exactly
four places (`Cell_dxdT.m`):

| # | Species | Donor wiring | Notes |
|---|---------|--------------|-------|
| 1 | **Acetyl-CoA** | `J_boxaccoa` added to `f(iACCOA_x)` (matrix AcCoA, idx 13); consumed by citrate synthase | Σ of all thiolase (MCKAT) + MTP fluxes; C4 thiolase ×2 |
| 2 | **NADH** | `J_boxnadh` added to `f(iNADH_x)` (matrix NADH, idx 10); consumed by Complex I | Σ of MSCHAD + MTP dehydrogenase fluxes; MSCHAD reads **free** NADH (`NADH_m`, line 627) |
| 3 | **CoA** | Conserved moiety: `CoAMAT = param(45) − Σ(all matrix acyl-CoA)` — *algebraic, not an ODE* | β-ox consumes/releases CoA via the acyl-CoA pools |
| 4 | **FADH₂** | **Not coupled to the ETC.** Drained by a 1st-order sink `vfadhsink = Ks·(FADH_m − K1)` (line 1737); `f(iQH2_x)` has **no** β-ox term | Acyl-CoA dehydrogenase electrons are **dropped** — they do *not* drive O₂/proton pumping |

**Implication of seam 4:** in the donor, fat's only energetic contribution to
OxPhos is through **acetyl-CoA → TCA → NADH (Complex I)** and the **MSCHAD
NADH**. The ~½ of fat-derived reducing equivalents that physiologically travel
via ETF→ubiquinone are *not* modeled — they vanish into the sink. Porting this
verbatim into Collins (which models the Q pool explicitly) would systematically
**under-predict fat-driven JO₂**. See Decision A.

The volume/unit bridge the donor uses for seams 1–2:
```
J_box•  =  ( (1/compartment_VMAT) · Σ reaction_v• ) / W_x
```
i.e. convert the module's per-mg-protein rate to per-L-mito, then to
per-L-matrix-water (`W_x`). **This factor is host-specific and must be rebuilt
for Collins.**

---

## 3. Convention mismatch: donor vs Collins

| Aspect | Donor (Meyer / van Eunen β-ox) | Collins (`dXdT.m`) | Consequence |
|--------|-------------------------------|--------------------|-------------|
| Matrix water | `W_x = 0.9·W_m ≈ 0.651 L/L mito` | `VWater_matrix = 0.4705 L/L mito` | Different normalization → conversion factor changes |
| β-ox flux basis | **per mg mito protein** | per L matrix water / per L mito | **Need mg-protein-per-L-mito constant** — the crux unit problem |
| Acetyl-CoA | shared matrix pool, idx 13 | explicit state `x(5)`, `f(5)` | Add source to `f(5)` |
| CoA | **conserved algebra** (`param(45) − Σ`) | **explicit state** `x(2)`, `f(2)` | Must wire CoA flux into `f(2)`, not just close a moiety |
| NAD/NADH | total `NADtot`; free NADH for MSCHAD | `NAD=x(3)`, `NADH=NADtot−NAD`; also vestigial `f(6)` | Add `−J_boxnadh` to `f(3)` |
| FADH₂ | own state + sink | no β-ox FAD pool; explicit Q pool (`coQH2=x(21)`) | Decision A: sink vs ETF→Q |
| Kinetics | plain MM, no pH/ΔΨ | full binding polynomials `P(i)`, `Keq(H,Mg,K)`, ΔΨ | Decision D: island vs reformulate |
| Proton balance | β-ox ignores H⁺ stoichiometry | careful matrix H⁺/charge balance (`f(54)`, Φ sums) | β-ox reactions perturb a balanced pH model |
| Compartments | cytosol + matrix | matrix + IM + **external buffer (cuvette)** — no cytosol | Map acyl-carnitine "cytosol" → buffer compartment |
| Substrate driver | malonyl-CoA/pyruvate/fat ICs | substrate combos from Excel | Add palmitoyl-carnitine as a buffer substrate |
| Extra Collins biology | — | ROS, NADPH, transhydrogenase, GDH, ammonia | Orthogonal to β-ox (no conflict; indirect redox coupling only) |

---

## 4. Concrete wiring points in Collins `dXdT.m`

Exact lines that would gain a β-ox source/sink term (current code shown):

```matlab
% line 2190  f(2)  COAS_matrix :  ... + J_SCS_matrix ) / VWater_matrix
%   → + (β-ox net CoA release/consumption)               ... / VWater_matrix
% line 2191  f(3)  NAD_matrix  :  ... + J_MPTT - J_MALIC ) / VWater_matrix
%   → − J_box_nadh_collins                                ... / VWater_matrix
% line 2193  f(5)  acetylcoA   :  + J_PDH - J_CTS )       / VWater_matrix
%   → + J_box_accoa_collins                               / VWater_matrix
% line 2212  f(21) coQH2 (mol/L mito, NOTE: no /VWater)   ← only if Decision A = couple
%   → + J_box_etfqo   (units: mol / L mito / s)
```

Note `f(21)` is normalized **per L mito**, not per L matrix water, and
oxaloacetate `f(12)` carries an extra `/beta` — unit care is non-negotiable.

**State-vector extension:** append the 47 β-ox states as **`x(70:116)`** (leave
Collins' 1–69 untouched — this avoids re-indexing every binding polynomial and
dissociation constant). FADH₂ becomes `x(116)` *or* is removed in favor of the
Q pool (Decision A). Update:
- the 47 new `f(70:116)` derivative lines (port from donor 1746–1790, re-based),
- `maincode.m` `x0(70:116)` initial conditions,
- the `J(...)` flux-output vector,
- `odeset('NonNegative', 1:62 …)` → extend to the new acyl-CoA states,
- the matrix proton/charge balance sums (Decision E).

---

## 5. Decisions someone must make

**A. FADH₂: keep the sink, or couple ETF→Q? (scientific)**
- *Sink (port as-is):* trivial, but fat JO₂ is under-counted and Collins' whole
  selling point is accurate respirometry. Defeats much of the purpose.
- *Couple:* add an ETF-QO-like reaction reducing `coQ→coQH2` (source on `f(21)`),
  remove the sink. More correct, but needs a new rate law + re-tuning. **This is
  the scientifically right but more expensive path.**

**B. CoA pool model.** Collins integrates CoA explicitly (`x(2)`); the donor
closes it by conservation. Porting requires writing the **exact CoA stoichiometry**
of CPT2 (consumes CoA) and the thiolases/MTP (release CoA) into `f(2)`, and
verifying total CoA (free + all acyl-CoA) stays conserved. Easy to get a slow
drift wrong here.

**C. Unit reconciliation (the crux).** Donor `boxVmax` is **per mg mito protein**;
Collins is volume-based. You must either (i) introduce a `mg protein / L mito`
density constant and convert all 12 `boxVmax` + the coupling fluxes, or
(ii) re-fit `boxVmax` directly in Collins units against fat-substrate data.
Because matrix water also differs (0.651 vs 0.4705), the naive transfer is wrong
by a constant — expect to **re-tune `boxVmax`** regardless.

**D. Thermodynamic treatment.** Cheapest: keep β-ox as a *kinetic island* with
its own fixed `Keq` and no ion binding (what Meyer did). "Correct": reformulate
all 12 reactions with Collins' binding-polynomial / ΔG°f machinery — large,
data-hungry effort. Recommend island first, reformulate only if needed.

**E. Proton / charge balance.** β-ox dehydrogenations release H⁺ that Collins'
matrix pH model would normally account for. The donor ignores this. Decide
whether to add proton stoichiometry to the coupling (consistent but laborious)
or accept a small pH/charge approximation (pragmatic).

**F. Compartment mapping.** No cytosol exists in Collins. Map the donor's
cytosolic acyl-carnitine + CPT1 to the **external buffer** compartment (CPT1 sits
on the outer membrane facing the medium in an isolated-mito prep), and supply
carnitine/CoA/malonyl-CoA as buffer constants or experimental inputs.

**G. Free vs total NADH.** Donor MSCHAD uses a *free* NADH (binding-corrected,
line 627). Collins' dehydrogenases use total `NADtot−NAD`. Pick one convention
for the ported MSCHAD/MTP to avoid a subtle bias.

---

## 6. Suggested phased plan

1. **Driver first.** Add palmitoyl-carnitine as a Collins buffer substrate +
   ICs, with no β-ox yet (sanity baseline). 
2. **Transplant the island.** Copy the 12 rate-law functions + parameter block;
   append `x(70:116)`; port the 47 derivative lines using the *donor's* internal
   `compartment_VMAT/VCYT` basis unchanged (self-consistent island).
3. **Wire the 4 seams** into `f(2)/f(3)/f(5)` with a single explicit
   per-mg-protein→Collins-units conversion constant (Decision C). FADH₂: start
   with the sink (Decision A = defer).
4. **Conserve & verify.** Check total CoA, total carnitine, NAD, FAD conservation
   numerically; confirm steady state with fat = 0 reproduces baseline Collins
   exactly (regression guard).
5. **Re-tune `boxVmax`** against a fat-substrate respirometry/NAD(P)H dataset.
6. **(Optional, higher fidelity)** Replace FADH₂ sink with ETF→Q coupling
   (Decision A) and re-validate JO₂.

---

## 7. Validation & expected differences

- **Conservation tests** (CoA, carnitine, NAD, Q, FAD) are the first-line bug net.
- **Regression:** with no fatty-acid substrate, the merged model **must** reduce
  numerically to stock Collins.
- **New observable:** palmitoyl-carnitine-driven JO₂ and NAD(P)H — the natural
  validation target (and where the FADH₂ decision shows up most).
- Expect β-ox-derived NADH to **shift the NADH/NAD redox poise**, which in
  Collins feeds back (indirectly) into **ROS and transhydrogenase/NADPH**
  predictions — a feature, but it means fat substrate changes more outputs than
  in the donor model.

---

## 8. Alternative worth weighing

Rather than porting β-ox *into* Collins, consider the reverse: the donor model
already integrates β-ox + whole-cell + the QSP layer. If the goal is "fat
metabolism **with** Collins' ROS/NADPH/redox detail," it may be cheaper to port
Collins' **ROS + transhydrogenase + explicit-CoA** features into the donor model
(which already has the β-ox seams wired and the cellular context) than to rebuild
β-ox, a cytosol, and a CoA pool inside the isolated-mitochondria Collins
framework. Pick the direction by which **end use** matters more: isolated-mito
respirometry (→ merge into Collins) vs. whole-heart energetics/pharmacology
(→ merge into the donor).
