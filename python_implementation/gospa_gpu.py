"""GPU (NVIDIA cuOpt) implementation of the graph GOSPA metric.

Drop-in replacement for graphGOSPA.LP_graph_GOSPA and
graphGOSPA.LP_graph_GOSPA_directed: same arguments, same return values.

Usage
-----
Put this file next to graphGOSPA.py and call it like the CPU version::

    import gospa_gpu as M
    dxy, loc, fa, miss, edge = M.LP_graph_GOSPA(
        X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon)

Requirements: numpy, scipy, cuopt (pip install cuopt-cu13), an NVIDIA GPU
and Linux. The LP is assembled by graphGOSPA._build_lp, so the model is
identical to graphGOSPA.py; only the solver changes.

Solver
------
method="pdlp" (default) runs cuOpt's GPU PDLP solver only. cuOpt's
"concurrent" mode also starts dual simplex and barrier on the CPU and waits
for them to stop, which makes small solves several times slower; on these
LPs dual simplex does not finish in reasonable time and barrier runs out of
memory, so they never win.

Accuracy
--------
presolve=True reproduces the CPU result to solver tolerance (tolerance,
default 1e-8). presolve=False is 4-8x faster on large graphs (n >= 300) but
the objective is accurate to ~1e-9 relative, not bit-identical.
presolve="auto" (default) turns presolve off above PRESOLVE_VAR_THRESHOLD.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sps

import graphGOSPA as _cpu

try:  # cuOpt is Linux-only; keep the import failure for a clear message.
    from cuopt.linear_programming import DataModel
    from cuopt.linear_programming.solver import Solve
    from cuopt.linear_programming.solver_settings import (SolverMethod,
                                                          SolverSettings)

    _CUOPT_IMPORT_ERROR = None
except Exception as _exc:  # pragma: no cover - depends on the machine
    DataModel = Solve = SolverSettings = SolverMethod = None
    _CUOPT_IMPORT_ERROR = _exc


PRESOLVE_VAR_THRESHOLD = 10_000
_INT32_MAX = np.iinfo(np.int32).max
_METHODS = {"pdlp": "PDLP", "concurrent": "Concurrent",
            "barrier": "Barrier", "dual_simplex": "DualSimplex"}


def gpu_available() -> bool:
    """True when cuOpt could be imported."""
    return _CUOPT_IMPORT_ERROR is None


def gpu_status() -> dict:
    """Diagnostic information about the GPU backend."""
    return {
        "available": gpu_available(),
        "import_error": None if gpu_available() else repr(_CUOPT_IMPORT_ERROR),
        "presolve_var_threshold": PRESOLVE_VAR_THRESHOLD,
    }


def _require_cuopt() -> None:
    if not gpu_available():
        raise ImportError(
            "cuOpt is not available. Install it with "
            "`pip install --extra-index-url https://pypi.nvidia.com cuopt-cu13` "
            "on Linux with an NVIDIA GPU."
        ) from _CUOPT_IMPORT_ERROR


def solve_lp(c, A_ub, b_ub, A_eq, b_eq, *, presolve="auto",
             tolerance=1e-8, time_limit=1800, method="pdlp"):
    """Solve min c'x s.t. A_ub x <= b_ub, A_eq x = b_eq, x >= 0.

    A_ub and b_ub may be None. method is one of "pdlp", "concurrent",
    "barrier" or "dual_simplex". Returns (x, objective). objective is
    evaluated as c @ x so it is always consistent with the returned vector.
    """
    _require_cuopt()
    if method not in _METHODS:
        raise ValueError(f"method must be one of {sorted(_METHODS)}")

    # cuOpt wants contiguous float64 vectors.
    c = np.ascontiguousarray(c, dtype=np.float64).ravel()
    n = c.size
    # The graph GOSPA LP built by graphGOSPA._build_lp has equality
    # constraints only; an empty inequality block keeps the code below generic.
    if A_ub is None:
        A_ub = sps.csr_matrix((0, n))
        b_ub = np.zeros(0)
    b_ub = np.ascontiguousarray(b_ub, dtype=np.float64).ravel()
    b_eq = np.ascontiguousarray(b_eq, dtype=np.float64).ravel()

    # Presolve gives tighter agreement with the CPU result but is slow on
    # large models, so "auto" enables it only for small ones.
    if presolve == "auto":
        presolve_on = n <= PRESOLVE_VAR_THRESHOLD
    else:
        presolve_on = bool(presolve)

    # cuOpt uses ranged rows, lower <= G x <= upper, so both blocks go into
    # one matrix:
    #   inequality rows: -inf <= A_ub x <= b_ub
    #   equality rows:   b_eq <= A_eq x <= b_eq
    G = sps.vstack([A_ub, A_eq], format="csr")
    # cuOpt takes int32 CSR indices; casting larger values would wrap silently.
    if G.nnz > _INT32_MAX or n > _INT32_MAX:
        raise ValueError(
            f"LP has {G.nnz} nonzeros, more than cuOpt's int32 limit "
            f"({_INT32_MAX})")
    lower = np.concatenate([np.full(A_ub.shape[0], -np.inf), b_eq])
    upper = np.concatenate([b_ub, b_eq])

    # Build the cuOpt model from the CSR arrays and the bounds.
    dm = DataModel()
    dm.set_csr_constraint_matrix(
        np.ascontiguousarray(G.data, dtype=np.float64),
        np.ascontiguousarray(G.indices, dtype=np.int32),
        np.ascontiguousarray(G.indptr, dtype=np.int32),
    )
    # The scipy matrix is no longer needed; free it before the solve so that
    # host memory peaks only once.
    del G
    dm.set_constraint_lower_bounds(lower)
    dm.set_constraint_upper_bounds(upper)
    # x >= 0, no upper bound (W <= 1 already follows from the row sums).
    dm.set_variable_lower_bounds(np.zeros(n))
    dm.set_variable_upper_bounds(np.full(n, np.inf))
    dm.set_objective_coefficients(c)
    dm.set_maximize(False)

    settings = SolverSettings()
    # Relative optimality tolerance of the solver.
    settings.set_optimality_tolerance(float(tolerance))
    settings.set_parameter("time_limit", int(time_limit))
    # PDLP only by default; see "Solver" in the module docstring.
    settings.set_parameter("method", int(getattr(SolverMethod, _METHODS[method])))
    if not presolve_on:
        settings.set_parameter("presolve", 0)

    sol = Solve(dm, settings)
    reason = str(sol.get_termination_reason())
    # A merely feasible point would overestimate the metric, so only accept
    # an optimal solution.
    if reason != "Optimal":
        raise RuntimeError("cuOpt did not solve the LP to optimality: " + reason)

    x = np.asarray(sol.get_primal_solution(), dtype=np.float64).ravel()
    return x, float(c @ x)


def _gpu_solver(presolve, tolerance, time_limit, method):
    # Adapt solve_lp to the solve(f, Aeq, beq) -> (x, objective) interface
    # expected by graphGOSPA._graph_GOSPA. Input checks, LP assembly and the
    # split into cost components are shared with the CPU backend, so the two
    # backends differ only in the solver.
    def solve(f, Aeq, beq):
        return solve_lp(f, None, None, Aeq, beq, presolve=presolve,
                        tolerance=tolerance, time_limit=time_limit,
                        method=method)
    return solve


def LP_graph_GOSPA(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, *,
                   presolve="auto", tolerance=1e-8, time_limit=1800,
                   method="pdlp"):
    """GPU version of graphGOSPA.LP_graph_GOSPA (undirected graphs)."""
    return _cpu._graph_GOSPA(
        X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, False,
        _gpu_solver(presolve, tolerance, time_limit, method))


def LP_graph_GOSPA_directed(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, *,
                            presolve="auto", tolerance=1e-8, time_limit=1800,
                            method="pdlp"):
    """GPU version of graphGOSPA.LP_graph_GOSPA_directed."""
    return _cpu._graph_GOSPA(
        X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, True,
        _gpu_solver(presolve, tolerance, time_limit, method))
