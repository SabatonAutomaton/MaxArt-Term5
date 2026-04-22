import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def objective(x, y):
    return 2 * (x - y) ** 2 + 14 * (y - 3) ** 2


def read_points(path):
    iterations = []
    xs = []
    ys = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                iterations.append(int(float(parts[0])))
                xs.append(float(parts[1]))
                ys.append(float(parts[2]))
            except ValueError:
                continue
    return np.array(iterations), np.array(xs), np.array(ys)


def boundary_line(constraint_type, x):
    if constraint_type == 1:
        return x + 1.0
    if constraint_type == 2:
        return -x
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--method", type=int, default=1)
    parser.add_argument("--constraint", type=int, default=0)
    parser.add_argument("--title", default="Trajectory of optimization")
    args = parser.parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output)
    iters, xs, ys = read_points(input_file)

    if len(xs) == 0:
        raise RuntimeError("Trajectory file is empty")

    x_min, x_max = xs.min(), xs.max()
    y_min, y_max = ys.min(), ys.max()
    x_pad = max(1.0, 0.25 * (x_max - x_min if x_max != x_min else 1.0))
    y_pad = max(1.0, 0.25 * (y_max - y_min if y_max != y_min else 1.0))

    x = np.linspace(x_min - x_pad, x_max + x_pad, 500)
    y = np.linspace(y_min - y_pad, y_max + y_pad, 500)
    X, Y = np.meshgrid(x, y)
    Z = objective(X, Y)

    plt.figure(figsize=(8, 6))
    levels = np.geomspace(max(Z.min(), 1e-6), max(Z.max(), 1e-6), 12)
    contour = plt.contour(X, Y, Z, levels=levels, colors="0.35")
    plt.clabel(contour, inline=True, fontsize=8)

    if args.constraint in (1, 2):
        bx = np.linspace(x.min(), x.max(), 400)
        by = boundary_line(args.constraint, bx)
        plt.plot(bx, by, "--", linewidth=2, color="tab:red")

    plt.plot(xs, ys, marker="o", linewidth=2, color="tab:blue")
    plt.scatter([xs[0]], [ys[0]], color="tab:green", s=70, zorder=3)
    plt.scatter([xs[-1]], [ys[-1]], color="tab:purple", s=70, zorder=3)

    plt.title(args.title)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, format="svg")
    plt.close()


if __name__ == "__main__":
    main()
