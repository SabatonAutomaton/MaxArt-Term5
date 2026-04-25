import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

# ==================== КОЭФФИЦИЕНТЫ ФУНКЦИИ (4-й вариант) ====================
# Соответствуют вашему Function.h
C = [4, 9, 1, 7, 5, 6]
a = [7, -9, 6, -8, -10, -2]
b = [9, -1, 5, -2, -8, -4]

def calculate_function(x, y):
    """
    Расчёт значения функции f(x,y) = sum(Ci / (1 + (x-ai)^2 + (y-bi)^2))
    Полностью соответствует вашей C++ функции
    """
    sum_val = 0
    for i in range(6):
        denominator = 1 + (x - a[i])**2 + (y - b[i])**2
        sum_val += C[i] / denominator
    return sum_val

# ==================== ФУНКЦИИ ДЛЯ ОТРИСОВКИ ====================

def plot_contours_only(title, x_min=-10, x_max=10, y_min=-10, y_max=10):
    """
    Рисует ТОЛЬКО изолинии функции (ландшафт)
    """
    grid_points = 300
    x_grid = np.linspace(x_min, x_max, grid_points)
    y_grid = np.linspace(y_min, y_max, grid_points)
    X, Y = np.meshgrid(x_grid, y_grid)
    
    # Вычисляем Z для сетки
    Z = np.zeros_like(X)
    for i in range(grid_points):
        for j in range(grid_points):
            Z[i, j] = calculate_function(X[i, j], Y[i, j])
    
    plt.figure(figsize=(10, 8))
    
    # Изолинии (контуры)
    contour_levels = 30
    contours = plt.contour(X, Y, Z, levels=contour_levels, 
                           cmap='viridis', linewidths=1.5)
    plt.clabel(contours, inline=True, fontsize=8, fmt='%.4f')
    
    # Заливка между изолиниями для наглядности
    plt.contourf(X, Y, Z, levels=contour_levels, cmap='viridis', alpha=0.6)
    
    # Границы области поиска
    plt.plot([x_min, x_max, x_max, x_min, x_min],
             [y_min, y_min, y_max, y_max, y_min], 'r-', linewidth=2.5, label='Search domain')
    
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel("X", fontsize=12)
    plt.ylabel("Y", fontsize=12)
    plt.axis('equal')
    plt.xlim(x_min - 0.5, x_max + 0.5)
    plt.ylim(y_min - 0.5, y_max + 0.5)
    
    # Цветовая шкала
    cbar = plt.colorbar(label='f(x, y)')
    
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right')
    
    output_filename = "function_contours.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_filename}")
    plt.show()


def plot_file_with_contours(filename, title, x_min=-10, x_max=10, y_min=-10, y_max=10):
    """
    Рисует график с изолиниями функции, точками и траекторией улучшений
    """
    # Типы точек
    R_x, R_y = [], []      # Случайные точки (Random)
    P_x, P_y = [], []      # Все точки (Points) - ложные + все пробы
    B_x, B_y = [], []      # Улучшения (Better)
    L_x, L_y = [], []      # Локальные минимумы (Local)
    S_x, S_y = [], []      # Стартовые точки (Start)

    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) == 3:
                    t, x, y = parts
                    x, y = float(x), float(y)
                    
                    if t == 'R':
                        R_x.append(x); R_y.append(y)
                    elif t == 'P':
                        P_x.append(x); P_y.append(y)
                    elif t == 'B':
                        B_x.append(x); B_y.append(y)
                    elif t == 'L':
                        L_x.append(x); L_y.append(y)
                    elif t == 'S':
                        S_x.append(x); S_y.append(y)
    except FileNotFoundError:
        print(f"File {filename} not found, skipping...")
        return

    print(f"Loaded from {filename}: R={len(R_x)}, P={len(P_x)}, B={len(B_x)}, L={len(L_x)}, S={len(S_x)}")

    # Создаём фигуру
    plt.figure(figsize=(12, 10))

    # ===== 1. ИЗОЛИНИИ ФУНКЦИИ =====
    grid_points = 200
    x_grid = np.linspace(x_min, x_max, grid_points)
    y_grid = np.linspace(y_min, y_max, grid_points)
    X, Y = np.meshgrid(x_grid, y_grid)
    
    # Вычисляем Z для сетки
    Z = np.zeros_like(X)
    for i in range(grid_points):
        for j in range(grid_points):
            Z[i, j] = calculate_function(X[i, j], Y[i, j])
    
    # Рисуем изолинии (20 уровней)
    contour_levels = 20
    contours = plt.contour(X, Y, Z, levels=contour_levels, 
                           cmap='viridis', alpha=0.7, linewidths=1)
    plt.clabel(contours, inline=True, fontsize=7, fmt='%.3f')
    
    # Полупрозрачная заливка
    plt.contourf(X, Y, Z, levels=contour_levels, cmap='viridis', alpha=0.3)

    # ===== 2. ВСЕ ТОЧКИ (ложные + случайные) - серые, мелкие =====
    all_points_x = []
    all_points_y = []
    if P_x:
        all_points_x.extend(P_x)
        all_points_y.extend(P_y)
    if R_x:
        all_points_x.extend(R_x)
        all_points_y.extend(R_y)
    
    if all_points_x:
        plt.scatter(all_points_x, all_points_y, s=3, label="All evaluated points", 
                   alpha=0.25, color='gray', zorder=1)

    # ===== 3. СТАРТОВАЯ ТОЧКА =====
    if S_x:
        plt.scatter(S_x, S_y, marker='s', label="Start point", s=180,
                   color='green', edgecolor='black', zorder=5, linewidth=2)

    # ===== 4. ТРАЕКТОРИЯ УЛУЧШЕНИЙ =====
    path_x = []
    path_y = []
    if S_x:
        path_x.append(S_x[0])
        path_y.append(S_y[0])
    if B_x:
        path_x.extend(B_x)
        path_y.extend(B_y)
    
    if len(path_x) > 1:
        plt.plot(path_x, path_y, marker='o', label="Improvement path", 
                color='blue', linewidth=2, markersize=6, zorder=4, alpha=0.9)
        
        # Добавляем стрелки направления
        for i in range(len(path_x)-1):
            plt.annotate('', xy=(path_x[i+1], path_y[i+1]), 
                        xytext=(path_x[i], path_y[i]),
                        arrowprops=dict(arrowstyle='->', color='blue', alpha=0.6, lw=1.5))

    # ===== 5. ЛОКАЛЬНЫЕ МИНИМУМЫ =====
    if L_x:
        plt.scatter(L_x, L_y, marker='*', label="Local minima", 
                   color='red', s=300, zorder=6, edgecolor='black', linewidth=1.5)
        
        # Подписываем значения минимумов
        for i, (lx, ly) in enumerate(zip(L_x, L_y)):
            f_val = calculate_function(lx, ly)
            plt.annotate(f'f={f_val:.5f}', xy=(lx, ly), 
                        xytext=(15, 10), textcoords='offset points',
                        fontsize=9, color='red', fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # ===== 6. ГРАНИЦЫ ОБЛАСТИ =====
    plt.plot([x_min, x_max, x_max, x_min, x_min],
             [y_min, y_min, y_max, y_max, y_min], 
             'k-', linewidth=2.5, label="Search domain")
    
    # ===== 7. НАСТРОЙКИ =====
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel("X", fontsize=12)
    plt.ylabel("Y", fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    plt.xlim(x_min - 0.5, x_max + 0.5)
    plt.ylim(y_min - 0.5, y_max + 0.5)
    
    # Цветовая шкала
    cbar = plt.colorbar(label='f(x, y)')
    
    # Легенда
    plt.legend(loc='upper right', framealpha=0.95, fontsize=10)
    
    # Сохранение
    output_filename = filename.replace(".txt", "_contours.png")
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_filename}")
    plt.show()


def plot_points_only(filename, title, x_min=-10, x_max=10, y_min=-10, y_max=10):
    """
    Рисует ТОЛЬКО точки (без изолиний) для быстрого просмотра
    """
    R_x, R_y = [], []
    P_x, P_y = [], []
    B_x, B_y = [], []
    L_x, L_y = [], []
    S_x, S_y = [], []

    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) == 3:
                    t, x, y = parts
                    x, y = float(x), float(y)
                    
                    if t == 'R':
                        R_x.append(x); R_y.append(y)
                    elif t == 'P':
                        P_x.append(x); P_y.append(y)
                    elif t == 'B':
                        B_x.append(x); B_y.append(y)
                    elif t == 'L':
                        L_x.append(x); L_y.append(y)
                    elif t == 'S':
                        S_x.append(x); S_y.append(y)
    except FileNotFoundError:
        print(f"File {filename} not found, skipping...")
        return

    plt.figure(figsize=(10, 8))

    # Все точки
    all_x = []
    all_y = []
    if P_x:
        all_x.extend(P_x); all_y.extend(P_y)
    if R_x:
        all_x.extend(R_x); all_y.extend(R_y)
    
    if all_x:
        plt.scatter(all_x, all_y, s=2, label="All points", alpha=0.2, color='gray')

    if S_x:
        plt.scatter(S_x, S_y, marker='s', label="Start", s=150, color='green', edgecolor='black')

    path_x = []
    path_y = []
    if S_x:
        path_x.append(S_x[0]); path_y.append(S_y[0])
    if B_x:
        path_x.extend(B_x); path_y.extend(B_y)
    if len(path_x) > 1:
        plt.plot(path_x, path_y, marker='o', label="Improvement path", color='blue', linewidth=2)

    if L_x:
        plt.scatter(L_x, L_y, marker='*', label="Local minima", color='red', s=250)

    plt.plot([x_min, x_max, x_max, x_min, x_min],
             [y_min, y_min, y_max, y_max, y_min], 'k-', linewidth=2)

    plt.title(title, fontsize=14)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right')
    
    output_filename = filename.replace(".txt", "_points.png")
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_filename}")
    plt.show()


# ==================== ОСНОВНОЙ КОД ====================

if __name__ == "__main__":
    X_MIN, X_MAX = -10, 10
    Y_MIN, Y_MAX = -10, 10

    # 1. Сначала рисуем изолинии самой функции
    print("\n=== Plotting function contours ===")
    plot_contours_only("Function Landscape f(x,y)", X_MIN, X_MAX, Y_MIN, Y_MAX)

    # 2. Файлы для отрисовки
    files = [
        ("simple.txt", "Simple Random Search"),
        ("alg1.txt", "Algorithm 1 (Multistart)"),
        ("alg2.txt", "Algorithm 2 (Sequential Improvement)"),
        ("alg3.txt", "Algorithm 3 (Gradient Descent)")
    ]

    for filename, title in files:
        print(f"\n=== Processing {filename} ===")
        plot_file_with_contours(filename, title, X_MIN, X_MAX, Y_MIN, Y_MAX)
        plot_points_only(filename, f"{title} (points only)", X_MIN, X_MAX, Y_MIN, Y_MAX)
    
    print("\n=== All plots completed! ===")