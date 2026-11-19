# Volumes of canonical knot complements: code and data

Code and data for the computations in

> M. Banarse, T. Cremaschi and J. McDonnell, *Volumes of canonical knot complements in tangent bundles* (2026).

For a closed geodesic γ on a hyperbolic surface S, the canonical lift γ̂ is the knot in the unit tangent bundle UT(S) made of the unit tangent vectors to γ. The repository computes the hyperbolic volume of its complement M_γ = UT(S) ∖ γ̂:

- **Modular surface.** All 2,807,194 closed geodesics on the modular surface whose code word has length at most 25 (Sections 3 and 6 of the paper).
- **Once-punctured torus.** Every filling closed curve of word length at most 11, computed through a finite cover of the modular surface.

**Interactive explorer:** <https://jackmcdonnell011-cmd.github.io/canonical-lift-volumes/explorer/volume_explorer.html>

## Data

| File | Contents |
|---|---|
| [`modular_geodesics_L25.csv.gz`](https://github.com/jackmcdonnell011-cmd/canonical-lift-volumes/releases/download/v1.0/modular_geodesics_L25.csv.gz) (release v1.0, 52 MB) | One row per closed geodesic on the modular surface with code word length 2 to 25 (2,807,194 rows).<br>Columns: `word, word_length, trace, hyperbolic_length, volume`.<br>SHA-256 `70a64ed3661f2182e0cb48f75e3eadf2aa06c15d6af9f2a5a05a00b9aebe44e4` |
| [`data/twist_parents.csv`](data/twist_parents.csv), [`data/twist_parents_summary.txt`](data/twist_parents_summary.txt) | The parent manifolds N of the twist families X^kY, X^kY², X^kY³ and the fillings of N that give the family members (the paper's table of parents). |
| [`data/torus11_volumes_N11.csv.gz`](data/torus11_volumes_N11.csv.gz) | Once-punctured torus: one row per filling closed curve of word length ≤ 11 in the free basis a, b of π₁(S) (12,047 curves).<br>Gives the curve's surface word, its modular code word w and degree m (the curve is a lift of the geodesic of wᵐ), the trace and length on S, the volume and the checks. |
| [`data/torus11_curves_N11.csv.gz`](data/torus11_curves_N11.csv.gz) | Once-punctured torus: every closed curve of word length ≤ 11, with whether it fills S (exact test). |
| [`data/torus11_twist_family_B2Ak.csv`](data/torus11_twist_family_B2Ak.csv) | The twist family B²Aᵏ on the once-punctured torus, k = 1, …, 24. |

**Conventions.**
- A code word is a cyclic word in X = [[1,1],[0,1]] and Y = [[1,0],[1,1]], listed in its lexicographically smallest rotation, one per closed geodesic.
- Hyperbolic length is 2 arccosh(|tr M_w| / 2).
- Volumes were computed with SnapPy.
- The DT codes are not stored, as they are regenerated in seconds by `modular/dt_codes.py`.

## Code

Python 3 with the packages in `requirements.txt`: `pip install -r requirements.txt`.

### `modular/`: the modular surface

By Ghys, UT(Mod) is the complement of the trefoil, and the canonical lift of a closed geodesic is a Lorenz knot. M_w is therefore the complement of a two-component link in S³, which SnapPy builds from a DT code.

| Script | What it does |
|---|---|
| `dt_codes.py` | The code words up to a given length (primitive words in their least rotation; 2,807,194 up to length 25), and `get_dt`, which turns a code word into the DT code of the Lorenz knot together with the trefoil (Williams' algorithm). |
| `volumes.py` | Volume (`calculate_knot_volume`) and hyperbolic length (`calculate_length`) for every code word up to a given length, in parallel. |
| `worked_example.py` | The paper's worked example w = XXXYXYY: ordering, DT code, volume 9.0467, trace 20, length 5.986. |
| `twist_parents.py` | The parent manifold N of each twist family is found by drilling the shortest geodesic out of M_w. The script checks that every member k = 1, …, 20 is a 1/n filling of N (`is_isometric_to`). |
| `make_figures.py` | The three figures of Section 3, from the release data file. |

```
cd modular
python worked_example.py
python volumes.py 12 volumes_L12.csv          # every geodesic up to code length 12
python twist_parents.py ../data
python make_figures.py modular_geodesics_L25.csv.gz figures
```

The original notebooks that produced the data set, `Modular_DT_Codes.ipynb` and `Modular Links Sage.ipynb`, are in `notebooks/`. The scripts above contain the same functions.

### `torus/`: the once-punctured torus

The once-punctured torus is S = H²/Γ′, where Γ′ is the commutator subgroup of PSL₂(ℤ), of index 6. UT(S) is then a 6-fold cover of UT(Mod), and M_γ is obtained from the 6-sheeted cover of M_w by Dehn filling every lift of the Lorenz knot other than γ̂ along its meridian.

`run_torus.sh <N> [processes]` runs the whole computation up to surface word length N:
1. `surface_census.py` lists every closed curve on S and tests whether it fills.
2. `run_surface.py` and `pipeline.py` build each manifold.
3. `merge_parts.py` merges the results and tests equal volumes for isometry.
4. `filling_exact.py` compares the results with the exact filling test.

Every manifold is kept only if it passes five checks:
1. The Ghys map ρ: π₁(M_w) → PSL₂(ℤ) is a representation.
2. The cover has the genus and cusps of S.
3. The core of the lift maps to a conjugate of M_wᵐ.
4. The lifts match the cycles of M_w on the cosets of Γ′.
5. Filling every lift gives back UT(S); this is compared through cusps, H₁ and H₁ of the double covers.

The curves that fill S are exactly those whose complement SnapPy finds hyperbolic, on every row.

For two curves SnapPy found no Ghys map, so their manifolds were taken from the reversed curve (`reverse_rho.py`). The involution v ↦ −v of UT(S) gives an isometry between M_γ and M_γ⁻¹.

`verify_sage.py` (run inside SageMath) recomputes each hyperbolic structure and volume with interval arithmetic. `twist_family.py` computes the family B²Aᵏ, and `explorer_json.py` writes the explorer's data from the results of `run_torus.sh`.

```
cd torus
bash run_torus.sh 7 4      # a small run: 108 code words, about a minute
```

### `explorer/`

`volume_explorer.html` is a self-contained page: open it in a browser, or use the link above. It is rebuilt with

```
python make_explorer.py modular_geodesics_L25.csv.gz volume_explorer.html torus11.json torus11_twists.json
```

## Licence

Code: MIT licence (see `LICENSE`). Data: CC BY 4.0. Please cite the paper above when using either.
