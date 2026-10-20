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
and Linux. CuPy is optional (it is installed together with cuopt-cu13).

LP assembly
-----------
build="gpu" assembles the LP on the GPU with CuPy (_build_lp_gpu). The
result is the same matrix as graphGOSPA._build_lp, entry for entry, but
it avoids SciPy's intermediate copies, which need several times the final
matrix in host memory. build="cpu" uses graphGOSPA._build_lp.
build="auto" (default) uses the GPU when CuPy is available and the LP has
at least GPU_BUILD_MIN_NONZEROS entries; smaller LPs build faster on the CPU.

Solver
------
method="pdlp" (default) runs cuOpt's GPU PDLP solver only. cuOpt's
"concurrent" mode also starts dual simplex and barrier on the CPU and waits
for them to stop, which makes small solves several times slower; on these
LPs dual simplex does not finish in reasonable time and barrier runs out of
memory, so they never win.

Accuracy
--------
Results agree with the CPU solver to the optimality tolerance (tolerance,
default 1e-8); the largest relative difference in dxy measured on random
graphs was 2.5e-9. Presolve is slow on large models, so presolve="auto"
(default) turns it off above PRESOLVE_VAR_THRESHOLD variables.
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

try:  # CuPy is only needed for build="gpu".
    import cupy as _cp
    import cupyx.scipy.sparse as _csp

    _CUPY_IMPORT_ERROR = None
except Exception as _exc:  # pragma: no cover - depends on the machine
    _cp = _csp = None
    _CUPY_IMPORT_ERROR = _exc


PRESOLVE_VAR_THRESHOLD = 10_000
_INT32_MAX = np.iinfo(np.int32).max
_METHODS = {"pdlp": "PDLP", "concurrent": "Concurrent",
            "barrier": "Barrier", "dual_simplex": "DualSimplex"}
_BUILDS = ("auto", "gpu", "cpu")
# Upper limit on COO entries generated per GPU chunk in _build_lp_gpu. The
# chunk is also capped by free GPU memory (about 64 bytes per entry).
BUILD_CHUNK_ENTRIES = 200_000_000
# build="auto" assembles on the GPU from this many stored entries on
# (about n = 215 for undirected graphs at edge density 0.5).
GPU_BUILD_MIN_NONZEROS = 10_000_000


def gpu_available() -> bool:
    """True when cuOpt could be imported."""
    return _CUOPT_IMPORT_ERROR is None


def cupy_available() -> bool:
    """True when CuPy could be imported (needed for build="gpu")."""
    return _CUPY_IMPORT_ERROR is None


def gpu_status() -> dict:
    """Diagnostic information about the GPU backend."""
    return {
        "available": gpu_available(),
        "import_error": None if gpu_available() else repr(_CUOPT_IMPORT_ERROR),
        "cupy_available": cupy_available(),
        "cupy_import_error": (None if cupy_available()
                              else repr(_CUPY_IMPORT_ERROR)),
        "presolve_var_threshold": PRESOLVE_VAR_THRESHOLD,
    }


def _require_cuopt() -> None:
    if not gpu_available():
        raise ImportError(
            "cuOpt is not available. Install it with "
            "`pip install --extra-index-url https://pypi.nvidia.com cuopt-cu13` "
            "on Linux with an NVIDIA GPU."
        ) from _CUOPT_IMPORT_ERROR


def _data_model(c, A_ub, b_ub, A_eq, b_eq):
    """Build the cuOpt DataModel for min c'x s.t. A_ub x <= b_ub,
    A_eq x = b_eq, x >= 0 (A_ub and b_ub may be None).

    The DataModel keeps its own host copy of the constraint matrix, so the
    caller can free its matrix once this returns.
    """
    # cuOpt wants contiguous float64 vectors.
    c = np.ascontiguousarray(c, dtype=np.float64).ravel()
    n = c.size
    # The graph GOSPA LP has equality constraints only (A_ub is None).
    if A_ub is None:
        A_ub = sps.csr_matrix((0, n))
        b_ub = np.zeros(0)
    b_ub = np.ascontiguousarray(b_ub, dtype=np.float64).ravel()
    b_eq = np.ascontiguousarray(b_eq, dtype=np.float64).ravel()

    # cuOpt uses ranged rows, lower <= G x <= upper, so both blocks go into
    # one matrix:
    #   inequality rows: -inf <= A_ub x <= b_ub
    #   equality rows:   b_eq <= A_eq x <= b_eq
    # Without inequality rows, use A_eq directly: vstack would copy the whole
    # matrix, which is the largest object in host memory.
    if A_ub.shape[0] == 0:
        G = A_eq.tocsr()
    else:
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
    # Drop the local reference; a copy made by vstack is freed here.
    del G
    dm.set_constraint_lower_bounds(lower)
    dm.set_constraint_upper_bounds(upper)
    # x >= 0, no upper bound (W <= 1 already follows from the row sums).
    dm.set_variable_lower_bounds(np.zeros(n))
    dm.set_variable_upper_bounds(np.full(n, np.inf))
    dm.set_objective_coefficients(c)
    dm.set_maximize(False)
    return dm


def _solve_data_model(dm, c, *, presolve, tolerance, time_limit, method):
    """Solve a DataModel from _data_model. Returns (x, c @ x)."""
    if method not in _METHODS:
        raise ValueError(f"method must be one of {sorted(_METHODS)}")
    c = np.ascontiguousarray(c, dtype=np.float64).ravel()
    # Presolve gives tighter agreement with the CPU result but is slow on
    # large models, so "auto" enables it only for small ones.
    if presolve == "auto":
        presolve_on = c.size <= PRESOLVE_VAR_THRESHOLD
    else:
        presolve_on = bool(presolve)

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
    dm = _data_model(c, A_ub, b_ub, A_eq, b_eq)
    return _solve_data_model(dm, c, presolve=presolve, tolerance=tolerance,
                             time_limit=time_limit, method=method)


def _mismatch_chunk(xi, xk, xv, yl, yj, yv, j0, j1, n_x, sign):
    """COO entries (on the GPU) of the rows of vec(A_X W_r - W_r A_Y) that
    belong to the mismatch columns j0 <= j < j1, times sign.

    Row (i, j) is row i + j*n_x of the block (column-major, as in
    graphGOSPA._edge_difference_matrix), and W(i, j) is variable i + j*(n_x+1).
      (A_X W_r)(i,j) = sum_k A_X(i,k) W(k,j): one entry per nonzero (i,k)
                       of A_X, in column k + j*(n_x+1)
      (W_r A_Y)(i,j) = sum_l W(i,l) A_Y(l,j): one entry per nonzero (l,j)
                       of A_Y and per i, in column i + l*(n_x+1)
    Row indices are local to the chunk (row 0 is row j0*n_x of the block).
    """
    xp = _cp
    J = xp.arange(j0, j1, dtype=xp.int32)
    # A_X W_r: every nonzero of A_X, repeated for each column j of the chunk
    r1 = (xi[None, :] + (J - j0)[:, None] * n_x).ravel()
    c1 = (xk[None, :] + J[:, None] * (n_x + 1)).ravel()
    v1 = xp.broadcast_to(sign * xv[None, :], (J.size, xv.size)).ravel()
    # W_r A_Y: the nonzeros (l, j) of A_Y with j in the chunk, for every row i
    sel = (yj >= j0) & (yj < j1)
    l, j, v = yl[sel], yj[sel], yv[sel]
    I = xp.arange(n_x, dtype=xp.int32)
    r2 = (I[None, :] + (j - j0)[:, None] * n_x).ravel()
    c2 = (I[None, :] + l[:, None] * (n_x + 1)).ravel()
    v2 = xp.broadcast_to(-sign * v[:, None], (v.size, n_x)).ravel()
    return (xp.concatenate((r1, r2)), xp.concatenate((c1, c2)),
            xp.concatenate((v1, v2)))


def _build_lp_gpu(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, directed):
    """Same LP as graphGOSPA._build_lp, assembled on the GPU with CuPy.

    The returned Aeq is a SciPy CSR matrix equal to the one from _build_lp
    (same rows, same column order, duplicates summed, zeros removed). The
    edge-mismatch rows are generated on the GPU in chunks of mismatch
    columns, converted to CSR there and copied into preallocated host
    arrays. Host memory therefore holds only the final matrix plus one chunk.
    """
    if not cupy_available():
        raise ImportError("build='gpu' needs CuPy") from _CUPY_IMPORT_ERROR
    xp = _cp
    n_x = X_attr.shape[0]
    n_y = Y_attr.shape[0]
    DAB = _cpu.locCostComp(X_attr, Y_attr, c, p)
    nxny = n_x * n_y
    WLen = (n_x + 1) * (n_y + 1)
    edge_weight = epsilon ** p / 4 if directed else epsilon ** p / 2
    nBlocks = 2 if directed else 1
    nSlack = nBlocks * nxny
    nParam = WLen + 2 * nSlack

    f = np.empty(nParam)
    f[:WLen] = DAB.reshape(WLen, order="F")
    f[WLen:] = edge_weight
    beq = np.concatenate((np.ones(n_x + n_y), np.zeros(nSlack)))

    X = np.asarray(X_adj, dtype=np.float64)
    Y = np.asarray(Y_adj, dtype=np.float64)
    # Block 1 is vec(A_X W_r - W_r A_Y). Block 2 (directed graphs only) is
    # vec(W_r A_Y^T - A_X^T W_r), which equals block 1 built from the
    # transposed adjacency matrices, times -1.
    blocks = [(X, Y, 1.0)]
    if directed:
        blocks.append((X.T, Y.T, -1.0))

    # Upper bound on the number of stored entries. It is exact unless a
    # self-loop puts two entries of one row in the same column.
    nnzX = np.count_nonzero(X)
    nnzY = np.count_nonzero(Y)
    nAssign = n_y * (n_x + 1) + n_x * (n_y + 1)
    bound = nAssign + nBlocks * (n_y * nnzX + n_x * nnzY + 2 * nxny)
    if bound > _INT32_MAX or nParam > _INT32_MAX:
        raise ValueError(
            f"LP has up to {bound} nonzeros, more than cuOpt's int32 limit "
            f"({_INT32_MAX})")

    nRows = n_y + n_x + nSlack
    data = np.empty(bound, dtype=np.float64)
    indices = np.empty(bound, dtype=np.int32)
    indptr = np.empty(nRows + 1, dtype=np.int32)

    # Assignment rows, in the order of graphGOSPA._assignment_constraints.
    # Row j < n_y sums column j of W, which is n_x+1 consecutive variables.
    # Row n_y+i sums row i of W, i.e. variable i + j*(n_x+1) for every j.
    k1 = n_y * (n_x + 1)
    indices[:k1] = np.arange(k1, dtype=np.int32)
    idx2 = (np.arange(n_x, dtype=np.int32)[:, None]
            + np.arange(n_y + 1, dtype=np.int32)[None, :] * (n_x + 1))
    indices[k1:nAssign] = idx2.ravel()
    data[:nAssign] = 1.0
    indptr[:n_y + 1] = np.arange(0, k1 + 1, n_x + 1)
    indptr[n_y + 1:n_y + n_x + 1] = k1 + np.arange(1, n_x + 1) * (n_y + 1)
    pos = nAssign
    row = n_y + n_x

    for b, (A, B, sign) in enumerate(blocks):
        xi, xk = np.nonzero(A)
        yl, yj = np.nonzero(B)
        xi_d = xp.asarray(xi, dtype=xp.int32)
        xk_d = xp.asarray(xk, dtype=xp.int32)
        xv_d = xp.asarray(A[xi, xk])
        yl_d = xp.asarray(yl, dtype=xp.int32)
        yj_d = xp.asarray(yj, dtype=xp.int32)
        yv_d = xp.asarray(B[yl, yj])
        # Entries per mismatch column j: nnz(A_X) + n_x*nnz(column j of A_Y)
        # + 2*n_x slack entries. The chunk is sized for the densest column.
        col_max = int(np.bincount(yj, minlength=n_y).max()) if yj.size else 0
        per_col = nnzX + n_x * col_max + 2 * n_x
        free = xp.cuda.Device().mem_info[0]
        budget = max(1, min(BUILD_CHUNK_ENTRIES, free // 64))
        step = max(1, budget // max(per_col, 1))
        for j0 in range(0, n_y, step):
            j1 = min(n_y, j0 + step)
            r, cidx, v = _mismatch_chunk(xi_d, xk_d, xv_d, yl_d, yj_d, yv_d,
                                         j0, j1, n_x, sign)
            # Slack columns of every row: s_plus with -1, s_minus with +1
            nr = (j1 - j0) * n_x
            rr = xp.arange(nr, dtype=xp.int32)
            slack = b * nxny + j0 * n_x + rr
            r = xp.concatenate((r, rr, rr))
            cidx = xp.concatenate((cidx, WLen + slack, WLen + nSlack + slack))
            v = xp.concatenate((v, xp.full(nr, -1.0), xp.full(nr, 1.0)))
            M = _csp.coo_matrix((v, (r, cidx)), shape=(nr, nParam))
            del r, cidx, v, rr, slack
            # Sum duplicates (a self-loop in both graphs puts two entries in
            # one column) and drop entries that cancel to zero. Zeros are
            # filtered on the COO arrays: cupyx's csr eliminate_zeros()
            # (CuPy 14.0.1) also drops negative values.
            M.sum_duplicates()
            keep = M.data != 0
            if not bool(keep.all()):
                M = _csp.coo_matrix((M.data[keep], (M.row[keep], M.col[keep])),
                                    shape=M.shape)
            del keep
            M = M.tocsr()
            if not M.has_sorted_indices:
                M.sort_indices()
            nnz = int(M.nnz)
            data[pos:pos + nnz] = M.data.get()
            indices[pos:pos + nnz] = M.indices.get()
            indptr[row + 1:row + nr + 1] = M.indptr[1:].get() + pos
            pos += nnz
            row += nr
            del M
        del xi_d, xk_d, xv_d, yl_d, yj_d, yv_d
    indptr[0] = 0
    # Give the GPU memory back before cuOpt allocates its own.
    xp.get_default_memory_pool().free_all_blocks()

    Aeq = sps.csr_matrix((data[:pos], indices[:pos], indptr),
                         shape=(nRows, nParam))
    return f, Aeq, beq, DAB, WLen, edge_weight


def _builder(build, n_x, n_y, X_adj, Y_adj, directed):
    # Pick the LP assembly: "gpu" (CuPy), "cpu" (SciPy) or "auto". "auto"
    # uses the GPU from GPU_BUILD_MIN_NONZEROS stored entries on; below that
    # the fixed cost of the first CuPy calls is larger than the gain.
    if build not in _BUILDS:
        raise ValueError(f"build must be one of {_BUILDS}")
    if build == "auto":
        nBlocks = 2 if directed else 1
        est = nBlocks * (n_y * np.count_nonzero(X_adj)
                         + n_x * np.count_nonzero(Y_adj))
        build = ("gpu" if cupy_available() and est >= GPU_BUILD_MIN_NONZEROS
                 else "cpu")
    return _build_lp_gpu if build == "gpu" else _cpu._build_lp


def _graph_GOSPA_gpu(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, directed,
                     presolve, tolerance, time_limit, method, build):
    # Same steps as graphGOSPA._graph_GOSPA. The difference is that the
    # host copy of the LP matrix is released once cuOpt's DataModel has
    # its own copy, so the matrix is held twice only for a moment.
    _require_cuopt()
    if method not in _METHODS:
        raise ValueError(f"method must be one of {sorted(_METHODS)}")
    X_attr, Y_attr, X_adj, Y_adj = _cpu._check_inputs(
        X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon)
    n_x = X_attr.shape[0]
    n_y = Y_attr.shape[0]
    if n_x == 0 or n_y == 0:
        return _cpu._empty_result(n_x, n_y, c, p)

    build_lp = _builder(build, n_x, n_y, X_adj, Y_adj, directed)
    f, Aeq, beq, DAB, WLen, edge_weight = build_lp(
        X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, directed)
    dm = _data_model(f, None, None, Aeq, beq)
    del Aeq
    x, obj = _solve_data_model(dm, f, presolve=presolve, tolerance=tolerance,
                               time_limit=time_limit, method=method)
    del dm
    return _cpu._components_from_solution(x, obj, DAB, n_x, n_y, WLen,
                                          edge_weight, p)


def LP_graph_GOSPA(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, *,
                   presolve="auto", tolerance=1e-8, time_limit=1800,
                   method="pdlp", build="auto"):
    """GPU version of graphGOSPA.LP_graph_GOSPA (undirected graphs)."""
    return _graph_GOSPA_gpu(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon,
                            False, presolve, tolerance, time_limit, method,
                            build)


def LP_graph_GOSPA_directed(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon, *,
                            presolve="auto", tolerance=1e-8, time_limit=1800,
                            method="pdlp", build="auto"):
    """GPU version of graphGOSPA.LP_graph_GOSPA_directed."""
    return _graph_GOSPA_gpu(X_attr, Y_attr, X_adj, Y_adj, c, p, epsilon,
                            True, presolve, tolerance, time_limit, method,
                            build)
