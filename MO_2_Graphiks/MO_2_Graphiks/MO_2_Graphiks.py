import numpy as np
import matplotlib.pyplot as plt


# =========================
# ФУНКЦИЯ (можно менять)
# =========================
def func(x, y):
    #return 100*(y-x)*(y-x) + (1-x)*(1-x)
    # return 100*(y-x*x)*(y-x*x) + (1-x)*(1-x)
     term1 = 2 / (1 + ((x - 1) / 2)**2 + ((y - 2) / 1)**2)
     term2 = 1 / (1 + ((x - 3) / 3)**2 + ((y - 1) / 3)**2)
     return term1 + term2
    # пример другой функции:
    # return x**2 + 2*y**2 + np.sin(x*y)


# =========================
# НАСТРОЙКИ
# =========================

INPUT_FILE = "out2.txt"
OUTPUT_FILE = "result.svg"

# если None → диапазон определяется автоматически
X_RANGE = None  # например (-5, 5)
Y_RANGE = None  # например (-5, 5)

GRID_DENSITY = 400
CONTOUR_LEVELS = 8


# =========================
# ЧТЕНИЕ ФАЙЛА
# =========================

def read_points(filename):
    iterations = []
    xs = []
    ys = []

    with open(filename, "r") as f:
        for line in f:
            if line.strip() == "":
                continue

            parts = line.split()
            it = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])

            iterations.append(it)
            xs.append(x)
            ys.append(y)

    return np.array(iterations), np.array(xs), np.array(ys)


# =========================
# АВТОМАТИЧЕСКИЙ ДИАПАЗОН
# =========================

def auto_range(values, padding_ratio=0.3):
    vmin = np.min(values)
    vmax = np.max(values)

    span = vmax - vmin
    if span == 0:
        span = 1

    padding = span * padding_ratio

    return vmin - padding, vmax + padding


# =========================
# ОСНОВНАЯ ПРОГРАММА
# =========================

def main():

    iterations, xs, ys = read_points(INPUT_FILE)

    # диапазон
    if X_RANGE is None:
        xmin, xmax = auto_range(xs)
    else:
        xmin, xmax = X_RANGE

    if Y_RANGE is None:
        ymin, ymax = auto_range(ys)
    else:
        ymin, ymax = Y_RANGE

    # сетка
    x = np.linspace(xmin, xmax, GRID_DENSITY)
    y = np.linspace(ymin, ymax, GRID_DENSITY)

    X, Y = np.meshgrid(x, y)
    Z = func(X, Y)

    # =========================
    # ГРАФИК
    # =========================

    plt.figure(figsize=(8, 6))

    # линии уровня
    z_min = np.min(Z)
    z_max = np.max(Z)

    levels = np.geomspace(z_min + 1e-6, z_max, CONTOUR_LEVELS)

    contours = plt.contour(X, Y, Z, levels=levels)
    plt.clabel(contours, inline=True, fontsize=8)

    # траектория
    plt.plot(xs, ys, marker='o')

    # подписи итераций
    # подпись первой точки
    plt.text(xs[0], ys[0], f"x_{iterations[0]}",
             fontsize=10, ha='right', va='bottom')

    # подпись последней точки
    plt.text(xs[-1], ys[-1], f"x_{iterations[-1]}",
             fontsize=10, ha='left', va='top')

    plt.title("Trajectory of Optimization Algorithm")
    plt.xlabel("x")
    plt.ylabel("y")

    plt.grid(True)

    # сохраняем SVG
    plt.savefig(OUTPUT_FILE, format="svg")

    print(f"График сохранён в {OUTPUT_FILE}")

    plt.show()


# =========================

if __name__ == "__main__":
    main()
