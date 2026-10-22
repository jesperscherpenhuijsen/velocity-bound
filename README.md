# velocity-bound

Code, results and manuscript for **"A velocity bound for positive-energy warp shells"** (J. E. Scherpenhuijsen, v4, October 2026).

- `papers/submission/velocity_bound.tex`: the manuscript.
- `code/`: the analysis pipeline. `wfcheck.py` is the core library (metric generators, Einstein tensor, energy conditions). The other scripts are listed in `code/README.md`.
- `results/`: raw JSON output and run logs that the paper's numbers come from.

## Reproducing the numbers

Requires Python 3 with numpy, scipy and sympy.

```bash
cd code
python3 vbound.py        # first-order thresholds, weak-field formula, bound, bisection (8 shells)
python3 exactJ.py        # exact 1D momentum density and thresholds
python3 mask_check.py pub 1e-6
```

Outputs go to `results/` by default. Set `FTL_RESULTS=/some/dir` to write elsewhere.

The Warp Factory toolkit used for comparison is at https://github.com/NerdsWithAttitudes/WarpFactory and is not included here.
