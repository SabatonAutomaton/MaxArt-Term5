# -*- coding: utf-8 -*-
from pathlib import Path
from pickle import TRUE

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent / "MO_3"
OUTPUT_FILE = SCRIPT_DIR / "result.svg"
PLOT_INFO = ROOT_DIR / "plot_info.txt"

USE_MANUAL_SCALE = 0
X_MIN = None
X_MAX = None
Y_MIN = None
Y_MAX = None
AUTO_PADDING = 0.25
GRID_DENSITY = 500
CONTOUR_LEVELS = 18
DEFAULT_METHOD = 1
DEFAULT_PROBLEM = 1


def objective(x, y, function_type):
    if function_type == 1:
        return 100 * (y - x) ** 2 + (1 - x) ** 2
    if function_type == 2:
        return 100 * (y - x * x) ** 2 + (1 - x) ** 2
    if function_type == 3:
        term1 = 2 / (1 + ((x - 1) / 2) ** 2 + (y - 2) ** 2)
        term2 = 1 / (1 + ((x - 3) / 3) ** 2 + ((y - 1) / 3) ** 2)
        return -(term1 + term2)
    if function_type == 4:
        return 2 * (x - y) ** 2 + 14 * (y - 3) ** 2
    raise ValueError(f"Unsupported function type: {function_type}")


def boundary_line(constraint_type, x):
    if constraint_type == 1:
        return x + 1.0
    if constraint_type == 2:
        return -x
    return None


def read_points(path, x_idx, y_idx):
    iterations = []
    xs = []
    ys = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) <= max(x_idx, y_idx):
                continue
            try:
                iterations.append(int(float(parts[0])))
                xs.append(float(parts[x_idx]))
                ys.append(float(parts[y_idx]))
            except ValueError:
                continue
    return np.array(iterations), np.array(xs), np.array(ys)


def read_plot_info():
    if not PLOT_INFO.exists():
        return DEFAULT_METHOD, DEFAULT_PROBLEM

    try:
        parts = PLOT_INFO.read_text(encoding="utf-8").split()
        return int(parts[0]), int(parts[1])
    except Exception:
        return DEFAULT_METHOD, DEFAULT_PROBLEM


def choose_input(method_type, problem_type):
    if method_type == 1:
        return ROOT_DIR / "out_penalty_traj.txt", problem_type, 4, 1, 2
    return ROOT_DIR / "out_barrier_traj.txt", 1, 4, 1, 2


def auto_bounds(xs, ys):
    x_min, x_max = np.min(xs), np.max(xs)
    y_min, y_max = np.min(ys), np.max(ys)
    x_span = x_max - x_min if x_max != x_min else 1.0
    y_span = y_max - y_min if y_max != y_min else 1.0
    pad = max(1.0, AUTO_PADDING * max(x_span, y_span))
    return x_min - pad, x_max + pad, y_min - pad, y_max + pad


def contour_levels(z, function_type):
    values = np.asarray(z, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.linspace(0.0, 1.0, CONTOUR_LEVELS)

    lo = float(np.quantile(values, 0.08))
    hi = float(np.quantile(values, 0.92))
    if np.isclose(lo, hi):
        lo = float(np.min(values))
        hi = float(np.max(values))
        if np.isclose(lo, hi):
            hi = lo + 1.0

    if function_type in (1, 2, 4) and lo > 0.0:
        lo = max(lo, 1e-6)
        hi = max(hi, lo * 1.01)
        return np.geomspace(lo, hi, CONTOUR_LEVELS)

    return np.linspace(lo, hi, CONTOUR_LEVELS)


method_type, problem_type = read_plot_info()
input_file, constraint_type, function_type, x_idx, y_idx = choose_input(method_type, problem_type)

iterations, xs, ys = read_points(input_file, x_idx, y_idx)
if len(xs) == 0:
    raise RuntimeError(f"Trajectory file is empty: {input_file}")

xmin, xmax, ymin, ymax = auto_bounds(xs, ys)
if USE_MANUAL_SCALE:
    if X_MIN is not None:
        xmin = X_MIN
    if X_MAX is not None:
        xmax = X_MAX
    if Y_MIN is not None:
        ymin = Y_MIN
    if Y_MAX is not None:
        ymax = Y_MAX

x = np.linspace(xmin, xmax, GRID_DENSITY)
y = np.linspace(ymin, ymax, GRID_DENSITY)
X, Y = np.meshgrid(x, y)
Z = objective(X, Y, function_type)

plt.figure(figsize=(8, 6))
levels = contour_levels(Z, function_type)
contours = plt.contour(X, Y, Z, levels=levels, colors="0.35")
plt.clabel(contours, inline=True, fontsize=8)

ax = plt.gca()
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

if constraint_type == 1:
    infeasible = (Y < X + 1.0).astype(float)
    plt.contourf(X, Y, infeasible, levels=[0.5, 1.5], colors=["tab:red"], alpha=0.12)
    bx = np.linspace(xmin, xmax, GRID_DENSITY)
    plt.plot(bx, boundary_line(constraint_type, bx), "--", color="tab:red", linewidth=2)
elif constraint_type == 2:
    infeasible = (np.abs(X + Y) > 1e-9).astype(float)
    plt.contourf(X, Y, infeasible, levels=[0.5, 1.5], colors=["tab:red"], alpha=0.08)
    bx = np.linspace(xmin, xmax, GRID_DENSITY)
    plt.plot(bx, boundary_line(constraint_type, bx), "--", color="tab:red", linewidth=2)

plt.plot(xs, ys, marker="o", color="tab:blue", linewidth=2)
plt.scatter([xs[0]], [ys[0]], color="tab:green", s=70, zorder=3)
plt.scatter([xs[-1]], [ys[-1]], color="tab:purple", s=70, zorder=3)
plt.annotate(r"$x_0$", (xs[0], ys[0]), textcoords="offset points", xytext=(-10, 10), fontsize=10)
plt.annotate(rf"$x_{{{iterations[-1]}}}$", (xs[-1], ys[-1]), textcoords="offset points", xytext=(10, -12), fontsize=10)

plt.title("Trajectory of optimization")
plt.xlabel("x")
plt.ylabel("y")
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.grid(True, alpha=0.25)
plt.tight_layout()
plt.savefig(OUTPUT_FILE, format="svg")
plt.close()
