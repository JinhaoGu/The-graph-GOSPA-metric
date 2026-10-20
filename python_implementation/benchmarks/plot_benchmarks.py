"""Draw the figures in GPU_BENCHMARKS.md from the measured results.

Usage: python plot_benchmarks.py
Reads results/*.json next to this file and writes figures/*.png.
The numbers were measured on an NVIDIA A100-SXM4-80GB (Google Colab) with
cuOpt 26.10, CuPy 14.0.1 and SciPy 1.16.3; see GPU_BENCHMARKS.md.
"""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullLocator

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)


def load(name):
    with open(os.path.join(HERE, "results", name + ".json"), encoding="utf-8") as f:
        return json.load(f)


DATA = {k: load(k) for k in ("scaling", "versions", "solver_modes", "kron_fix", "max_n")}
SMALL = load("small_graphs")

# Validated categorical palette (fixed order), light surface
S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
MUTED = "#b5b3ab"
SURF, INK, INK2, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
MARK = ["o", "s", "^", "D", "v", "P"]
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "DejaVu Sans", "font.size": 10, "text.color": INK,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "axes.titlesize": 11, "axes.titleweight": "medium", "axes.titlelocation": "left", "axes.titlepad": 10,
    "legend.frameon": False, "legend.fontsize": 9,
})
GRAPH = {"Undirected": 0, "Directed": 1, "Complete": 2}


def fmt_s(v, _=None):
    return f"{v:,.0f} s" if v >= 10 else (f"{v:.1f} s" if v >= 1 else f"{v:.2f} s")


def finish(fig, name, note=None):
    if note:
        fig.text(0.01, 0.005, note, fontsize=8, color=INK2, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.04 if note else 0, 1, 1))
    path = os.path.join(OUT, name + ".png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def end_labels(ax, items, x_pad=0.012, min_gap_frac=0.055, log=False):
    # items: [(label, x, y, color)]; nudge apart so labels never touch
    import numpy as np
    lo, hi = ax.get_ylim()
    tr = (lambda v: np.log10(v)) if log else (lambda v: v)
    inv = (lambda v: 10 ** v) if log else (lambda v: v)
    span = tr(hi) - tr(lo)
    pts = sorted(((tr(y), lbl, x, c) for lbl, x, y, c in items))
    placed = []
    for ty, lbl, x, c in pts:
        if placed and ty - placed[-1][0] < min_gap_frac * span:
            ty = placed[-1][0] + min_gap_frac * span
        placed.append((ty, lbl, x, c))
    xr = ax.get_xlim()
    for ty, lbl, x, c in placed:
        ax.annotate(lbl, (x, inv(ty)), xytext=(6, 0), textcoords="offset points", va="center",
                    fontsize=9, color=INK, annotation_clip=False)


rows = DATA["scaling"]
by = lambda g: sorted([r for r in rows if r["graph"] == g], key=lambda r: r["n"])
paths = []

# 1. Total time vs n
fig, ax = plt.subplots(figsize=(7.2, 4.2))
items = []
for g, i in GRAPH.items():
    rr = by(g)
    ax.plot([r["n"] for r in rr], [r["total_s"] for r in rr], color=S[i], lw=2, marker=MARK[i], ms=6,
            mec=SURF, mew=1.2, label=g, solid_capstyle="round")
    items.append((f"{g}  {fmt_s(rr[-1]['total_s'])}", rr[-1]["n"], rr[-1]["total_s"], S[i]))
ax.set_yscale("log")
ax.yaxis.set_major_formatter(FuncFormatter(fmt_s))
ax.yaxis.set_minor_locator(NullLocator())
ax.set_xlabel("Nodes per graph (n)")
ax.set_title("Total time per graph pair")
ax.set_xlim(60, 1350)
end_labels(ax, items, log=True)
ax.legend(loc="upper left")
paths.append(finish(fig, "time_vs_n", "A100 80 GB, default settings (PDLP, build=\"auto\"). Edge density 0.5; complete graphs density 1. Log scale."))

# 2. Build vs solve, undirected
fig, ax = plt.subplots(figsize=(7.2, 4.0))
rr = by("Undirected")
xs = range(len(rr))
b = [r["build_s"] for r in rr]
s = [r["solve_s"] for r in rr]
ax.bar(xs, b, width=0.62, color=S[3], label="Build LP", edgecolor=SURF, linewidth=1.5)
ax.bar(xs, s, width=0.62, bottom=b, color=S[0], label="cuOpt model + solve", edgecolor=SURF, linewidth=1.5)
ax.set_xticks(list(xs), [f"{r['n']}" + ("\nCPU" if r["build_on"] == "CPU" else "") for r in rr])
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f} s"))
ax.set_xlabel("Nodes per graph (n), undirected")
ax.set_title("Where the time goes")
ax.grid(axis="x", visible=False)
ax.legend(loc="upper left")
for x, r in zip(xs, rr):
    if r["n"] in (700, 1100):
        ax.annotate(f"{r['build_s']/r['total_s']:.0%} build", (x, r["total_s"]), xytext=(0, 4), textcoords="offset points",
                    ha="center", fontsize=8, color=INK2)
paths.append(finish(fig, "build_vs_solve", "The LP is built on the CPU up to n = 200 and on the GPU (CuPy) from n = 300."))

# 3. Memory: VRAM and RAM side by side (same x, separate y)
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.0), sharex=True)
for ax, key, title in ((axs[0], "vram_gb", "GPU memory (peak VRAM)"), (axs[1], "ram_gb", "Host memory (peak RAM)")):
    items = []
    for g, i in GRAPH.items():
        rr = by(g)
        ax.plot([r["n"] for r in rr], [r[key] for r in rr], color=S[i], lw=2, marker=MARK[i], ms=6, mec=SURF, mew=1.2, label=g)
        items.append((f"{rr[-1][key]:.0f} GB", rr[-1]["n"], rr[-1][key], S[i]))
    if key == "vram_gb":
        ax.axhline(80, color=INK2, lw=1, ls=(0, (4, 4)))
        ax.text(1240, 81.5, "80 GB (A100)", fontsize=8, color=INK2, ha="right")
        ax.set_ylim(0, 92)
    else:
        ax.set_ylim(0, 92)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f} GB"))
    ax.set_title(title)
    ax.set_xlabel("Nodes per graph (n)")
    ax.set_xlim(60, 1250)
    end_labels(ax, items)
axs[0].legend(loc="upper left", bbox_to_anchor=(0.0, 0.8))
paths.append(finish(fig, "memory_vs_n", "Same runs as the time figure. Both axes use the same scale, so host RAM reads as about 0.4x the VRAM."))

# 4. What each change did
V = ["Original", "#3 L1 + PDLP", "#4 Memory fixes", "#5 GPU assembly"]
VC = [MUTED, S[1], S[2], S[0]]
ver = DATA["versions"]
ns = sorted({r["n"] for r in ver})
fig, axs = plt.subplots(1, 3, figsize=(11, 4.0))
for ax, key, title, unit in ((axs[0], "total_s", "Total time", "s"), (axs[1], "ram_gb", "Host RAM", "GB"), (axs[2], "vram_gb", "GPU memory", "GB")):
    w = 0.2
    top = max(r[key] for r in ver if r[key] is not None)
    for vi, v in enumerate(V):
        for ni, n in enumerate(ns):
            r = next(r for r in ver if r["version"] == v and r["n"] == n)
            x = ni + (vi - 1.5) * w
            if r[key] is None:
                ax.text(x, top * 0.02, "did not fit\n(> 80 GB)", rotation=90, ha="center", va="bottom", fontsize=7.5, color=INK2)
            else:
                ax.bar(x, r[key], width=w * 0.92, color=VC[vi], edgecolor=SURF, linewidth=1, label=v if ni == 0 else None)
    ax.set_xticks(range(len(ns)), [f"n = {n:,}" for n in ns])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _, u=unit: f"{v:,.0f} {u}"))
    ax.set_title(title)
    ax.grid(axis="x", visible=False)
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0.04, 1, 0.92))
fig.text(0.01, 0.005, "Undirected graphs, density 0.5, same instances for every version. Original = main before PR #3 (concurrent solver, old LP).",
         fontsize=8, color=INK2)
p = os.path.join(OUT, "changes_by_pr.png"); fig.savefig(p, dpi=160); plt.close(fig); paths.append(p)

# 5. Many small graphs
MODES = [("gpu_seq", "GPU, one at a time"), ("gpu_thr2", "GPU, 2 threads"), ("gpu_thr4", "GPU, 4 threads"),
         ("gpu_thr8", "GPU, 8 threads"), ("cpu_seq", "CPU, one process"), ("cpu_proc12", "CPU, 12 processes")]
fig, ax = plt.subplots(figsize=(7.6, 4.4))
items = []
for i, (m, lbl) in enumerate(MODES):
    rr = sorted([r for r in SMALL if r.get("mode") == m and "per_pair_ms" in r], key=lambda r: r["n"])
    if not rr:
        continue
    ax.plot([r["n"] for r in rr], [r["per_pair_ms"] for r in rr], color=S[i], lw=2, marker=MARK[i], ms=6, mec=SURF, mew=1.2,
            label=lbl, ls="-" if m.startswith("gpu") else (0, (5, 2)))
ax.set_yscale("log")
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f} ms" if v >= 1 else f"{v:.1f} ms"))
ax.yaxis.set_minor_locator(NullLocator())
ax.set_xticks(sorted({r["n"] for r in SMALL if "per_pair_ms" in r}))
ax.set_xlabel("Nodes per graph (n)")
ax.set_title("Many small graph pairs: time per pair (lower is better)")
ax.set_xlim(15, 105)
ax.legend(loc="upper left", ncol=2)
paths.append(finish(fig, "small_graphs", "Undirected, 32 pairs per size (fewer for slow CPU runs). GPU runs exclude the first call's start-up. Dashed = CPU (HiGHS). Log scale."))

# 6. Solver mode
sm = DATA["solver_modes"]
fig, ax = plt.subplots(figsize=(6.0, 3.6))
for j, (meth, col) in enumerate((("Concurrent", MUTED), ("PDLP only", S[0]))):
    vals = [next(r["solve_s"] for r in sm if r["n"] == n and r["method"] == meth) for n in (300, 700)]
    xs = [k + (j - 0.5) * 0.36 for k in range(2)]
    ax.bar(xs, vals, width=0.34, color=col, edgecolor=SURF, label=meth)
    for x, v in zip(xs, vals):
        ax.annotate(fmt_s(v), (x, v), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8.5, color=INK)
ax.set_xticks([0, 1], ["n = 300", "n = 700"])
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f} s"))
ax.set_title("Solve time: concurrent mode vs PDLP only")
ax.grid(axis="x", visible=False)
ax.legend(loc="upper left")
paths.append(finish(fig, "solver_mode", "Concurrent also starts dual simplex and barrier on the CPU and waits for them; PDLP always finished first."))

# 7. kron fix (complete graph n = 100)
kf = {r["code"]: r for r in DATA["kron_fix"]}
fig, axs = plt.subplots(1, 2, figsize=(9.0, 2.9))
lab = ["Before the fix", "After the fix"]
axs[0].barh([1, 0], [600, kf["after_fix"]["total_s"]], color=[MUTED, S[0]], height=0.5)
axs[0].set_xscale("log"); axs[0].set_xlim(0.5, 2000)
axs[0].set_yticks([1, 0], lab)
axs[0].xaxis.set_major_formatter(FuncFormatter(fmt_s)); axs[0].xaxis.set_minor_locator(NullLocator())
axs[0].text(0.6, 1, " did not finish in 600 s", va="center", fontsize=8.5, color=INK)
axs[0].text(kf["after_fix"]["total_s"], 0, f"  {kf['after_fix']['total_s']:.1f} s", va="center", fontsize=8.5, color=INK)
axs[0].set_title("Time, complete graph n = 100 (log)"); axs[0].grid(axis="y", visible=False)
axs[1].barh([1, 0], [kf["before_fix"]["stored_entries_m"], kf["after_fix"]["stored_entries_m"]], color=[MUTED, S[0]], height=0.5)
axs[1].set_yticks([1, 0], lab)
axs[1].xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f} M"))
for y, v in ((1, kf["before_fix"]["stored_entries_m"]), (0, kf["after_fix"]["stored_entries_m"])):
    axs[1].text(v, y, f"  {v:,.1f} M", va="center", fontsize=8.5, color=INK)
axs[1].set_xlim(0, 125)
axs[1].set_title("Stored entries in the LP matrix"); axs[1].grid(axis="y", visible=False)
paths.append(finish(fig, "kron_fix", "SciPy's kron returned BSR blocks that stored zeros for dense graphs. The real number of nonzeros is 2.0 M."))

# 8. Largest graph per GPU memory
mx = DATA["max_n"]
gbs = sorted({r["gpu_gb"] for r in mx})
fig, ax = plt.subplots(figsize=(7.2, 3.8))
for j, (g, i) in enumerate(GRAPH.items()):
    vals = [next(r["max_n"] for r in mx if r["gpu_gb"] == gb and r["graph"] == g) for gb in gbs]
    xs = [k + (j - 1) * 0.26 for k in range(len(gbs))]
    ax.bar(xs, vals, width=0.24, color=S[i], edgecolor=SURF, label=g)
    if j == 0:
        for x, v in zip(xs, vals):
            ax.annotate(f"{v:,}", (x, v), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8, color=INK)
ax.set_xticks(range(len(gbs)), [f"{gb} GB" for gb in gbs])
ax.set_xlabel("GPU memory")
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
ax.set_title("Largest graph that fits (nodes per graph)")
ax.grid(axis="x", visible=False)
ax.legend(loc="upper left")
paths.append(finish(fig, "max_n_per_gpu", "From VRAM fits; measured up to n = 1100 undirected, 800 directed, 700 complete. Larger values are extrapolated."))

for p in paths:
    print(os.path.relpath(p, HERE))
