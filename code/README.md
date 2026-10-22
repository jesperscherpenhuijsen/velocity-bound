# Code

Python 3 with numpy, scipy, sympy. All scripts expect to sit in one folder (they import wfcheck).

| Script | What it does | Writes |
| --- | --- | --- |
| wfcheck.py | Core pipeline: metric generators (ported from Warp Factory, incl. Warp Shell TOV construction), 4th-order finite-difference Einstein tensor, Eulerian tetrad, sampled and Hawking-Ellis energy conditions | (library) |
| formula_test.py, formula_test2.py | Symbolic check of the Ricci-tensor formula used in wfcheck.py | stdout |
| validate.py | Validation: exact sympy Einstein tensor for Alcubierre, closed-form energy density, Schwarzschild vacuum, Minkowski | stdout |
| run_checks.py, run_mt.py | Phase 0: every metric shipped with Warp Factory, energy conditions | results.json, results_mt.json |
| shell_diag.py, shell_scan.py | Phase 0: Warp Shell velocity threshold diagnosis and scan | stdout |
| appendixA_check.py | Phase 0 Appendix A identities | stdout |
| vbound.py | Velocity-bound paper: first-order threshold, weak-field formula, bound, bisection, for 8 shells | vbound_results.json |
| vbound_rest.py, vbound_last.py, vbound_fix.py | Helper runners used to finish the 8-shell scan within time limits | vbound_results.json |
| exactJ.py | Exact linearized momentum density on a static spherically symmetric background (sympy) and 1D thresholds | exactJ_results.json, exactJ_expressions.txt |
| mask_check.py | Matter-mask sensitivity check (usage: python3 mask_check.py pub 1e-6, or c07 1e-6) | mask_*.json |
| lapse_identity.py | Paper 3: numerical check of the lapse-weighted null-average identity | stdout |
| dbg.py | Debug script from the validation phase | stdout |

Warp Factory reference: https://github.com/NerdsWithAttitudes/WarpFactory (main branch, cloned Oct 20, 2026; MIT licence). Not included here.
