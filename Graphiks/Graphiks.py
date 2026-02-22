# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


def parabola(x):
    return (x - 3)**2


def read_input(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip() != ""]

    a0 = float(lines[0])
    b0 = float(lines[1])
    x_values = [float(x) for x in lines[2:]]

    return a0, b0, x_values


def main():
    #filename = "dichotomy.txt"
    #filename = "goldenRatio.txt"
    #filename = "fibonacci.txt"
    a0, b0, x_points = read_input(filename)

    x = np.linspace(a0, b0, 1000)
    y = parabola(x)

    plt.figure(figsize=(8, 10))  # вытянутый вверх график

    plt.plot(x, y, label="y = (x - 3)^2")

    colors = ['red', 'blue', 'green', 'orange',
              'purple', 'brown', 'pink', 'gray']

    xtick_positions = []
    xtick_labels = []

    # ---- границы ----
    for boundary, label in [(a0, "a_0"), (b0, "b_0")]:
        yb = parabola(boundary)
        plt.plot([boundary, boundary], [0, yb],
                 linestyle='--',
                 color='black',
                 linewidth=1)

        xtick_positions.append(boundary)
        xtick_labels.append(label)

    # ---- точки ----
    index = 1
    for i in range(0, len(x_points), 2):
        pair = x_points[i:i+2]
        color = colors[(i // 2) % len(colors)]

        for xp in pair:
            yp = parabola(xp)
            plt.scatter(xp, yp, color=color, s=60)
            plt.plot([xp, xp], [0, yp],
                     linestyle='--',
                     color=color,
                     linewidth=1)

            if index <= 4:
                xtick_positions.append(xp)
                xtick_labels.append(f"x_{index}")

            index += 1

    # подписи на оси X
    plt.xticks(xtick_positions, xtick_labels, rotation=45)

    # оси
    plt.axhline(0, color='black', linewidth=1)
    plt.axvline(0, color='black', linewidth=1)

    # ---- границы X с небольшим отступом ----
    dx = b0 - a0
    plt.xlim(a0 - dx*0.05, b0 + dx*0.05)

    # немного воздуха сверху
    plt.ylim(bottom=0, top=max(y)*1.05)

    plt.title("График y = (x - 3)^2")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()