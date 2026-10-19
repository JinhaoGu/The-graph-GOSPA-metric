"""GPU (NVIDIA cuOpt) implementation of the graph GOSPA metric.

Drop-in replacement for graphGOSPA.LP_graph_GOSPA and
graphGOSPA.LP_graph_GOSPA_directed: same arguments, same return values.

Usage
-----
Put this file next to graphGOSPA.py and call it like the CPU version::

    import graphGOSPA_gpu as M
    dxy, loc, fa, miss, edge = M.LP_graph_GOSPA(
        X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon)

Requirements: numpy, scipy, cuopt (pip install cuopt-cu13), an NVIDIA GPU
and Linux. The LP is assembled by the CPU module's own helpers, so the model
is identical to graphGOSPA.py; only the solver changes.

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
    from cuopt.linear_programming.solver_settings import SolverSettings

    _CUOPT_IMPORT_ERROR = None
except Exception as _exc:  # pragma: no cover - depends on the machine
    DataModel = Solve = SolverSettings = None
    _CUOPT_IMPORT_ERROR = _exc


PRESOLVE_VAR_THRESHOLD = 10_000


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
             tolerance=1e-8, time_limit=1800):
    """Solve min c'x s.t. A_ub x <= b_ub, A_eq x = b_eq, x >= 0.

    Returns (x, objective). objective is evaluated as c @ x so it is always
    consistent with the returned vector.
    """
    _require_cuopt()

    c = np.ascontiguousarray(c, dtype=np.float64).ravel()
    b_ub = np.ascontiguousarray(b_ub, dtype=np.float64).ravel()
    b_eq = np.ascontiguousarray(b_eq, dtype=np.float64).ravel()
    A_ub = A_ub.tocsr()
    A_eq = A_eq.tocsr()

    n = c.size
    if presolve == "auto":
        presolve_on = n <= PRESOLVE_VAR_THRESHOLD
    else:
        presolve_on = bool(presolve)

    G = sps.vstack([A_ub, A_eq], format="csr")
    lower = np.concatenate([np.full(A_ub.shape[0], -np.inf), b_eq])
    upper = np.concatenate([b_ub, b_eq])

    dm = DataModel()
    dm.set_csr_constraint_matrix(
        np.ascontiguousarray(G.data, dtype=np.float64),
        np.ascontiguousarray(G.indices, dtype=np.int32),
        np.ascontiguousarray(G.indptr, dtype=np.int32),
    )
    dm.set_constraint_lower_bounds(lower)
    dm.set_constraint_upper_bounds(upper)
    dm.set_variable_lower_bounds(np.zeros(n))
    dm.set_variable_upper_bounds(np.full(n, np.inf))
    dm.set_objective_coefficients(c)
    dm.set_maximize(False)

    settings = SolverSettings()
    settings.set_optimality_tolerance(float(tolerance))
    settings.set_parameter("time_limit", int(time_limit))
    if not presolve_on:
        settings.set_parameter("presolve", 0)

    sol = Solve(dm, settings)
    reason = str(sol.get_termination_reason())
    if reason not in ("Optimal", "PrimalFeasible"):
        raise RuntimeError("cuOpt did not solve the LP: " + reason)

    x = np.asarray(sol.get_primal_solution(), dtype=np.float64).ravel()
    return x, float(c @ x)


def _build_undirected(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon):
    """Mirror of graphGOSPA.LP_graph_GOSPA's model assembly."""
    n_x = X_attr.shape[0]
    n_y = Y_attr.shape[0]
    DAB = _cpu.locCostComp(X_attr, Y_attr, c, p)

    nxny = n_x * n_y
    WLen = (n_x + 1) * (n_y + 1)
    nParam = WLen + 1 + nxny

    f = np.zeros(nParam)
    f[:WLen] = DAB.reshape(WLen, order="F")
    f[WLen] = epsilon ** p / 2

    Aeq, beq = _cpu._assignment_constraints(n_x, n_y, nParam)

    A1 = sps.hstack((
        sps.csr_matrix((1, WLen)),
        sps.csr_matrix(np.array([[-1.0]])),
        np.ones((1, nxny)),
    ), format="csr")

    G = sps.hstack((
        _cpu._edge_difference_matrix(n_x, n_y, X_adj, Y_adj),
        sps.csr_matrix((nxny, n_x + 1)),
    ), format="csr")
    neg_I = -sps.eye(nxny, format="csr")
    A2 = sps.hstack((G, sps.csr_matrix((nxny, 1)), neg_I), format="csr")
    A3 = sps.hstack((-G, sps.csr_matrix((nxny, 1)), neg_I), format="csr")

    A = sps.vstack((A1, A2, A3), format="csr")
    b = np.zeros(A.shape[0])
    return f, A, b, Aeq, beq, DAB, n_x, n_y, WLen


def _build_directed(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon):
    """Mirror of graphGOSPA.LP_graph_GOSPA_directed's model assembly."""
    n_x = X_attr.shape[0]
    n_y = Y_attr.shape[0]
    DAB = _cpu.locCostComp(X_attr, Y_attr, c, p)

    nxny = n_x * n_y
    WLen = (n_x + 1) * (n_y + 1)
    nParam = WLen + 2 + 2 * nxny

    f = np.zeros(nParam)
    f[:WLen] = DAB.reshape(WLen, order="F")
    f[WLen] = epsilon ** p / 4
    f[WLen + 1] = epsilon ** p / 4

    Aeq, beq = _cpu._assignment_constraints(n_x, n_y, nParam)

    minus_one = sps.csr_matrix(np.array([[-1.0]]))
    A1_1 = sps.hstack((
        sps.csr_matrix((1, WLen)), minus_one, sps.csr_matrix((1, 1)),
        np.ones((1, nxny)), sps.csr_matrix((1, nxny)),
    ), format="csr")
    A1_2 = sps.hstack((
        sps.csr_matrix((1, WLen)), sps.csr_matrix((1, 1)), minus_one,
        sps.csr_matrix((1, nxny)), np.ones((1, nxny)),
    ), format="csr")
    A1 = sps.vstack((A1_1, A1_2), format="csr")

    G1 = sps.hstack((
        _cpu._edge_difference_matrix(n_x, n_y, X_adj, Y_adj),
        sps.csr_matrix((nxny, n_x + 1)),
    ), format="csr")
    G2 = sps.hstack((
        _cpu._reverse_edge_difference_matrix(n_x, n_y, X_adj, Y_adj),
        sps.csr_matrix((nxny, n_x + 1)),
    ), format="csr")

    zeros_e = sps.csr_matrix((nxny, 2))
    neg_I_h1 = sps.hstack((-sps.eye(nxny, format="csr"),
                           sps.csr_matrix((nxny, nxny))), format="csr")
    neg_I_h2 = sps.hstack((sps.csr_matrix((nxny, nxny)),
                           -sps.eye(nxny, format="csr")), format="csr")

    A2 = sps.hstack((G1, zeros_e, neg_I_h1), format="csr")
    A3 = sps.hstack((-G1, zeros_e, neg_I_h1), format="csr")
    A4 = sps.hstack((G2, zeros_e, neg_I_h2), format="csr")
    A5 = sps.hstack((-G2, zeros_e, neg_I_h2), format="csr")

    A = sps.vstack((A1, A2, A3, A4, A5), format="csr")
    b = np.zeros(A.shape[0])
    return f, A, b, Aeq, beq, DAB, n_x, n_y, WLen


def LP_graph_GOSPA(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, *,
                   presolve="auto", tolerance=1e-8, time_limit=1800):
    """GPU version of graphGOSPA.LP_graph_GOSPA (undirected graphs)."""
    X_attr, Y_attr, X_adj, Y_adj = _cpu._check_inputs(
        X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon)
    n_x = X_attr.shape[0]
    n_y = Y_attr.shape[0]
    if n_x == 0 or n_y == 0:
        return _cpu._empty_result(n_x, n_y, c, p)

    f, A, b, Aeq, beq, DAB, n_x, n_y, WLen = _build_undirected(
        X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon)
    x, obj = solve_lp(f, A, b, Aeq, beq, presolve=presolve,
                      tolerance=tolerance, time_limit=time_limit)

    Wx = x[:WLen].reshape((n_x + 1, n_y + 1), order="F")
    loc_cost = np.sum(DAB[:n_x, :n_y] * Wx[:n_x, :n_y])
    false_cost = np.sum(DAB[n_x, :n_y] * Wx[n_x, :n_y])
    miss_cost = np.sum(DAB[:n_x, n_y] * Wx[:n_x, n_y])
    edge_cost = epsilon ** p / 2 * x[WLen]

    return _cpu._metric_components(obj, loc_cost, false_cost, miss_cost,
                                   edge_cost, p)


def LP_graph_GOSPA_directed(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, *,
                            presolve="auto", tolerance=1e-8, time_limit=1800):
    """GPU version of graphGOSPA.LP_graph_GOSPA_directed."""
    X_attr, Y_attr, X_adj, Y_adj = _cpu._check_inputs(
        X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon)
    n_x = X_attr.shape[0]
    n_y = Y_attr.shape[0]
    if n_x == 0 or n_y == 0:
        return _cpu._empty_result(n_x, n_y, c, p)

    f, A, b, Aeq, beq, DAB, n_x, n_y, WLen = _build_directed(
        X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon)
    x, obj = solve_lp(f, A, b, Aeq, beq, presolve=presolve,
                      tolerance=tolerance, time_limit=time_limit)

    Wx = x[:WLen].reshape((n_x + 1, n_y + 1), order="F")
    loc_cost = np.sum(DAB[:n_x, :n_y] * Wx[:n_x, :n_y])
    false_cost = np.sum(DAB[n_x, :n_y] * Wx[n_x, :n_y])
    miss_cost = np.sum(DAB[:n_x, n_y] * Wx[:n_x, n_y])
    edge_cost = epsilon ** p / 4 * (x[WLen] + x[WLen + 1])

    return _cpu._metric_components(obj, loc_cost, false_cost, miss_cost,
                                   edge_cost, p)
