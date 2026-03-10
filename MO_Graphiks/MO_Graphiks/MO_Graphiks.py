import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# ГЛОБАЛЬНЫЕ МАТЕМАТИЧЕСКИЕ КОНСТАНТЫ
# ============================================================
PI = np.pi
TAU = 2 * np.pi
E = np.e
INF = np.inf

PI_2 = PI / 2
PI_4 = PI / 4
PI_3 = PI / 3
PI_6 = PI / 6
PI_8 = PI / 8

SQRT2 = np.sqrt(2)
SQRT3 = np.sqrt(3)
SQRT2_2 = SQRT2 / 2
SQRT3_2 = SQRT3 / 2

DEG_TO_RAD = PI / 180
RAD_TO_DEG = 180 / PI
# ============================================================


def read_input_file(filename='input.txt'):
    """
    Чтение параметров из файла input.txt
    Формат файла:
    x0 x1
    a b
    ellipse_params (необязательно: a_ell b_ell)
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()

        lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]

        x0, x1 = map(float, lines[0].split())
        a, b = map(float, lines[1].split())
        
        ellipse_params = None
        line_idx = 2
        
        if len(lines) > 2:
            parts = lines[2].split()
            if len(parts) == 2 and all(p.replace('-', '').replace('.', '').isdigit() for p in parts):
                ellipse_params = tuple(map(float, parts))
                line_idx = 3

        return x0, x1, a, b, ellipse_params

    except Exception as e:
        print(f"Ошибка при чтении файла: {e}")
        print("Использую значения по умолчанию")
        return -1.0, 1.0, 0.0, 2.0, (1.0, 1.0)


def choose_function():
    """
    Меню выбора типа функции
    """
    print("\n" + "="*50)
    print("ВЫБОР ТИПА ФУНКЦИИ")
    print("="*50)
    print("1. y = C * tg(x)")
    print("2. y = C * cos(x)")
    print("3. y = (x - C)^3")
    print("4. y = C * (x^2 - 2x)")
    print("5. y = C * sin(x - π/4)")
    print("6. y = e^(x + C)")
    print("="*50)
    
    while True:
        try:
            choice = int(input("Введите номер функции (1-6): "))
            if 1 <= choice <= 6:
                return choice
            else:
                print("Ошибка: введите число от 1 до 6")
        except ValueError:
            print("Ошибка: введите целое число")


def get_function_name(choice):
    """
    Возвращает название функции для заголовка
    """
    names = {
        1: "y = C \\cdot tg x",
        2: "y = C \\cdot \\cos x",
        3: "y = (x - C)^3",
        4: "y = C \\cdot (x^2 - 2x)",
        5: "y = C \\cdot \\sin(x - \\frac{\\pi}{4})",
        6: "y = e^{x+C}"
    }
    return names.get(choice, "y = f(x, C)")


def right_function(C, x, choice):
    """
    Вычисление функции в зависимости от выбора
    """
    if choice == 1:
        # y = C * tg(x)
        return C * np.tan(x)
    
    elif choice == 2:
        # y = C * cos(x)
        return C * np.cos(x)
    
    elif choice == 3:
        # y = (x - C)^3
        return (x - C)**3
    
    elif choice == 4:
        # y = C * (x^2 - 2x)
        return C * (x**2 - 2*x)
    
    elif choice == 5:
        # y = C * sin(x - π/4)
        return C * np.sin(x - PI/4)
    
    elif choice == 6:
        # y = e^(x + C)
        return np.exp(x + C)
    
    else:
        # По умолчанию
        return (x - C)**3


def get_function_label(C, choice):
    """
    Возвращает форматированную подпись для функции
    """
    if choice == 1:
        return f'$y = {C:g} \\cdot tg x$'
    elif choice == 2:
        return f'$y = {C:g} \\cdot \\cos x$'
    elif choice == 3:
        return f'$y = (x - {C:g})^3$'
    elif choice == 4:
        return f'$y = {C:g} \\cdot (x^2 - 2x)$'
    elif choice == 5:
        return f'$y = {C:g} \\cdot \\sin(x - \\frac{{\\pi}}{{4}})$'
    elif choice == 6:
        return f'$y = e^{{x + {C:g}}}$'
    else:
        return f'$y = f(x, {C:g})$'


def calculate_C_values(a, b):
    """
    Вычисление 5 значений C
    """
    C_pos = [a, (a+b)/2, b]
    C_neg = [-c for c in C_pos if c != 0]

    C_values = sorted(C_neg + C_pos)
    return C_values


def f(C, x, choice, ellipse_params=None):
    """
    Основная функция с проверкой области
    """
    y = right_function(C, x, choice)
    
    if ellipse_params is not None:
        a_ell, b_ell = ellipse_params
        # Проверка нахождения внутри эллипса/окружности
        mask = (x**2 / a_ell**2 + y**2 / b_ell**2) <= 1
        return np.where(mask, y, np.nan)
    
    return y


def calculate_y_limits(x0, x1, C_values, choice, ellipse_params=None, padding=0.15):
    """
    Вычисление границ y0 и y1 с учетом области
    """
    x = np.linspace(x0, x1, 1000)
    all_y = []

    for C in C_values:
        y = right_function(C, x, choice)
        
        if ellipse_params is not None:
            a_ell, b_ell = ellipse_params
            mask = (x**2 / a_ell**2 + y**2 / b_ell**2) <= 1
            if np.any(mask):
                all_y.extend(y[mask])
        else:
            all_y.extend(y)

    if not all_y:
        if ellipse_params:
            return -ellipse_params[1], ellipse_params[1]
        return -1.0, 1.0

    y_min = np.min(all_y)
    y_max = np.max(all_y)

    y_range = y_max - y_min

    y0 = y_min - padding * y_range
    y1 = y_max + padding * y_range

    return y0, y1


def place_labels_safe(target_positions, min_dist, y_range):
    """
    Размещение подписей без изменения порядка.
    """
    positions = target_positions.copy()

    for i in range(len(positions)):
        if abs(positions[i]) < 0.04 * y_range:
            positions[i] = 0.04 * y_range

    for i in range(1, len(positions)):
        if positions[i-1] - positions[i] < min_dist:
            positions[i] = positions[i-1] - min_dist

    return positions


def find_intersections(x0, x1, C_values, a, choice, ellipse_params=None, num_points=10000):
    """
    Находит точки пересечения функции с C = a и другими функциями
    """
    x = np.linspace(x0, x1, num_points)
    y_a = f(a, x, choice, ellipse_params)
    
    intersections = []
    
    for C in C_values:
        if abs(C - a) < 1e-10:
            continue
            
        y_other = f(C, x, choice, ellipse_params)
        
        valid_mask = ~(np.isnan(y_a) | np.isnan(y_other))
        if not np.any(valid_mask):
            continue
            
        x_valid = x[valid_mask]
        y_a_valid = y_a[valid_mask]
        y_other_valid = y_other[valid_mask]
        
        diff = y_a_valid - y_other_valid
        sign_changes = np.where(np.diff(np.signbit(diff)))[0]
        
        for idx in sign_changes:
            x1_idx = x_valid[idx]
            x2_idx = x_valid[idx + 1]
            y1_diff = diff[idx]
            y2_diff = diff[idx + 1]
            
            if abs(y2_diff - y1_diff) > 1e-10:
                x_intersect = x1_idx - y1_diff * (x2_idx - x1_idx) / (y2_diff - y1_diff)
            else:
                x_intersect = (x1_idx + x2_idx) / 2
            
            if x0 <= x_intersect <= x1:
                is_duplicate = False
                for existing_x in intersections:
                    if abs(existing_x - x_intersect) < 1e-6 * (x1 - x0):
                        is_duplicate = True
                        break
                
                if not is_duplicate:
                    intersections.append(x_intersect)
    
    return sorted(intersections)


def find_closest_point_in_domain(C, x_boundary, x0, x1, choice, ellipse_params, direction='left'):
    """
    Находит ближайшую точку внутри области от границы
    """
    if direction == 'left':
        x_search = np.linspace(x_boundary, x1, 200)
    else:
        x_search = np.linspace(x0, x_boundary, 200)
    
    y_search = f(C, x_search, choice, ellipse_params)
    valid_idx = np.where(~np.isnan(y_search))[0]
    
    if len(valid_idx) > 0:
        if direction == 'left':
            return x_search[valid_idx[0]], y_search[valid_idx[0]]
        else:
            return x_search[valid_idx[-1]], y_search[valid_idx[-1]]
    
    return None, None


def draw_graphs(x0, x1, a, b, choice, ellipse_params):
    """
    Отрисовка графиков с настраиваемыми параметрами отступов
    """
    
    # ============================================================
    # ГЛОБАЛЬНЫЕ ПАРАМЕТРЫ ДЛЯ НАСТРОЙКИ ВСЕХ ПОДПИСЕЙ И ОТСТУПОВ
    # ============================================================

    # ===== ЛЕВЫЕ ПОДПИСИ (C = ...) =====
    left_label_x_offset = -0.2         # смещение по X: отрицательное = левее, положительное = правее
    left_label_y_offset = 0.22           # смещение по Y: положительное = выше, отрицательное = ниже
    left_label_fontsize = 10              # размер шрифта

    # ===== ПРАВЫЕ ПОДПИСИ (y = ...) =====
    right_label_x_offset = 0.35           # смещение по X: положительное = правее, отрицательное = левее
    right_label_y_offset = 0.13           # смещение по Y: положительное = выше, отрицательное = ниже
    right_label_fontsize = 9               # размер шрифта

    # ===== РАССТОЯНИЕ МЕЖДУ ПОДПИСЯМИ =====
    label_vertical_spacing = 0.18         # меньше = плотнее, больше = разреженнее

    # ===== ПОДПИСИ ГРАНИЦ ОБЛАСТИ =====
    # x₀
    x0_label_x_offset = -0.02             # смещение x₀ по X (отрицательное = левее)
    x0_label_y_offset = -0.02             # смещение x₀ по Y (отрицательное = ниже)
    # x₁
    x1_label_x_offset = 0.02              # смещение x₁ по X (положительное = правее)
    x1_label_y_offset = -0.02             # смещение x₁ по Y (отрицательное = ниже)
    # y₀
    y0_label_x_offset = -0.02             # смещение y₀ по X (отрицательное = левее)
    y0_label_y_offset = -0.02             # смещение y₀ по Y (отрицательное = ниже)
    # y₁
    y1_label_x_offset = -0.02             # смещение y₁ по X (отрицательное = левее)
    y1_label_y_offset = 0.02              # смещение y₁ по Y (положительное = выше)

    # ===== ПОДПИСИ ОСЕЙ X И Y =====
    x_axis_label_offset = 0.05            # смещение подписи "x" от конца оси (положительное = правее)
    y_axis_label_offset = 0.05            # смещение подписи "y" от конца оси (положительное = выше)

    # ===== ЦИФРЫ НА ОСИ X =====
    x_tick_x_offset_pos = 0.015           # смещение положительных цифр по X (вправо)
    x_tick_x_offset_neg = -0.015          # смещение отрицательных цифр по X (влево)
    x_tick_y_offset = -0.03               # смещение цифр по Y (отрицательное = ниже оси)

    # ===== ЦИФРЫ НА ОСИ Y =====
    y_tick_y_offset_pos = 0.015           # смещение положительных цифр по Y (вверх)
    y_tick_y_offset_neg = -0.015          # смещение отрицательных цифр по Y (вниз)
    y_tick_x_offset = -0.03               # смещение цифр по X (отрицательное = левее оси)

    # ===== ПОДПИСЬ НУЛЯ =====
    zero_x_offset = -0.02                  # смещение нуля по X (отрицательное = левее)
    zero_y_offset = -0.02                  # смещение нуля по Y (отрицательное = ниже)

    # ===== СТРЕЛКИ ОСЕЙ =====
    arrow_head_width = 0.02                # ширина стрелки (доля от размера области)
    arrow_head_length = 0.02               # длина стрелки (доля от размера области)

    # ===== ПУНКТИРНЫЕ ЛИНИИ ГРАНИЦ =====
    dash_line_width = 1.5                  # толщина пунктира
    dash_alpha = 0.6                       # прозрачность (0-1, где 1 = непрозрачный)

    # ===== ОТСТУПЫ ДЛЯ РИСОВАНИЯ ГРАФИКОВ =====
    plot_x_left_padding = 0.15              # отступ слева (доля от ширины области)
    plot_x_right_padding = 0.3             # отступ справа (доля от ширины области)
    plot_y_bottom_padding = 0.2            # отступ снизу (доля от высоты области)
    plot_y_top_padding = 0.2               # отступ сверху (доля от высоты области)

    # ===== ИТОГОВЫЕ ОТСТУПЫ ГРАФИКА =====
    final_x_padding = 0.1                   # итоговый отступ по X (доля от ширины области)
    final_y_padding = 0.1                   # итоговый отступ по Y (доля от высоты области)

    # ===== ПУЧКИ (ЖИРНЫЕ ТОЧКИ) =====
    bundle_markersize = 8                    # размер точек пучков
    bundle_tolerance = 0.001                 # точность определения пересечений (меньше = точнее)

    # ============================================================
    # КОНЕЦ ГЛОБАЛЬНЫХ ПАРАМЕТРОВ
    # ============================================================

    C_values = calculate_C_values(a, b)
    print(f"Значения C: {[round(c, 3) for c in C_values]}")

    y0, y1 = calculate_y_limits(x0, x1, C_values, choice, ellipse_params, padding=0.15)

    print(f"Границы области: x∈[{x0:.3f}, {x1:.3f}], y∈[{y0:.3f}, {y1:.3f}]")
    if ellipse_params:
        print(f"Область: x²/{ellipse_params[0]}² + y²/{ellipse_params[1]}² ≤ 1")

    x_range = x1 - x0
    y_range = y1 - y0

    x_extended = np.linspace(x0 - plot_x_left_padding*x_range, 
                             x1 + plot_x_right_padding*x_range, 10000)

    fig, ax = plt.subplots(figsize=(12, 8))

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(axis='both', which='both', length=0)

    for spine in ax.spines.values():
        spine.set_visible(False)

    x_axis_start = min(x0 - plot_x_left_padding*x_range, 0)
    x_axis_end = max(x1 + plot_x_right_padding*x_range, 0)

    y_axis_start = min(y0 - plot_y_bottom_padding*y_range, 0)
    y_axis_end = max(y1 + plot_y_top_padding*y_range, 0)

    # Рисуем оси
    ax.arrow(
        x_axis_start, 0,
        x_axis_end - x_axis_start, 0,
        head_width=arrow_head_width*y_range,
        head_length=arrow_head_length*x_range,
        fc='black', ec='black', linewidth=1,
        length_includes_head=True, clip_on=False
    )

    ax.arrow(
        0, y_axis_start,
        0, y_axis_end - y_axis_start,
        head_width=arrow_head_width*x_range,
        head_length=arrow_head_length*y_range,
        fc='black', ec='black', linewidth=1,
        length_includes_head=True, clip_on=False
    )

    # Рисуем границы области
    if ellipse_params:
        a_ell, b_ell = ellipse_params
        theta = np.linspace(0, 2*PI, 200)
        ellipse_x = a_ell * np.cos(theta)
        ellipse_y = b_ell * np.sin(theta)
        ax.plot(ellipse_x, ellipse_y, color='gray', linestyle='--', 
               linewidth=dash_line_width, alpha=dash_alpha, zorder=1)
        
        x0_display, x1_display = -a_ell, a_ell
        y0_display, y1_display = -b_ell, b_ell
    else:
        if abs(x0) > 1e-10:
            ax.axvline(x=x0, color='gray', linestyle='--', 
                      linewidth=dash_line_width, alpha=dash_alpha, zorder=1)
        if abs(x1) > 1e-10:
            ax.axvline(x=x1, color='gray', linestyle='--', 
                      linewidth=dash_line_width, alpha=dash_alpha, zorder=1)
        if abs(y0) > 1e-10:
            ax.axhline(y=y0, color='gray', linestyle='--', 
                      linewidth=dash_line_width, alpha=dash_alpha, zorder=1)
        if abs(y1) > 1e-10:
            ax.axhline(y=y1, color='gray', linestyle='--', 
                      linewidth=dash_line_width, alpha=dash_alpha, zorder=1)
        
        x0_display, x1_display = x0, x1
        y0_display, y1_display = y0, y1

    # Масштабные координаты на осях
    x_tick_step = 0.5
    x_ticks = []
    tick = np.floor(min(x0_display, 0) / x_tick_step) * x_tick_step

    while tick <= max(x1_display, 0) + 0.001:
        if abs(tick) > 1e-10:
            x_ticks.append(tick)
        tick += x_tick_step

    for tick in x_ticks:
        if x_axis_start <= tick <= x_axis_end:
            ax.plot(tick, 0, 'ko', markersize=3)
            if tick > 0:
                ax.text(tick + x_tick_x_offset_pos*x_range, x_tick_y_offset*y_range, 
                       f'{tick:g}', fontsize=9, ha='center', va='top')
            else:
                ax.text(tick + x_tick_x_offset_neg*x_range, x_tick_y_offset*y_range, 
                       f'{tick:g}', fontsize=9, ha='center', va='top')

    y_tick_step = 0.5
    y_ticks = []
    tick = np.floor(min(y0_display, 0) / y_tick_step) * y_tick_step

    while tick <= max(y1_display, 0) + 0.001:
        if abs(tick) > 1e-10:
            y_ticks.append(tick)
        tick += y_tick_step

    for tick in y_ticks:
        if y_axis_start <= tick <= y_axis_end:
            ax.plot(0, tick, 'ko', markersize=3)
            if tick > 0:
                ax.text(y_tick_x_offset*x_range, tick + y_tick_y_offset_pos*y_range, 
                       f'{tick:g}', fontsize=9, ha='right', va='center')
            else:
                ax.text(y_tick_x_offset*x_range, tick + y_tick_y_offset_neg*y_range, 
                       f'{tick:g}', fontsize=9, ha='right', va='center')

    ax.text(zero_x_offset*x_range, zero_y_offset*y_range, '0', 
           fontsize=9, ha='right', va='top', zorder=30)

    # Подписываем границы области
    ax.text(x0_display + x0_label_x_offset*x_range, x0_label_y_offset*y_range, '$x_0$', 
           fontsize=11, ha='right', va='top', fontweight='bold', zorder=30)
    ax.text(x1_display + x1_label_x_offset*x_range, x1_label_y_offset*y_range, '$x_1$', 
           fontsize=11, ha='left', va='top', fontweight='bold', zorder=30)

    ax.text(y0_label_x_offset*x_range, y0_display + y0_label_y_offset*y_range, '$y_0$', 
           fontsize=11, ha='right', va='top', fontweight='bold', zorder=30)
    ax.text(y1_label_x_offset*x_range, y1_display + y1_label_y_offset*y_range, '$y_1$', 
           fontsize=11, ha='right', va='bottom', fontweight='bold', zorder=30)

    # Подписываем оси
    ax.text(x_axis_end + x_axis_label_offset*x_range, -0.01*y_range, '$x$', 
           fontsize=14, fontweight='bold')
    ax.text(-0.01*x_range, y_axis_end + y_axis_label_offset*y_range, '$y$', 
           fontsize=14, fontweight='bold')

    # Рисуем функции
    for C in C_values:
        y_extended = f(C, x_extended, choice, ellipse_params)
        ax.plot(x_extended, y_extended, color='black', linewidth=1.5, 
               linestyle='-', zorder=5)

    # Левые подписи
    left_pairs = []
    for C in C_values:
        y_at_x0 = f(C, x0_display, choice, ellipse_params)
        
        if not np.isnan(y_at_x0):
            left_pairs.append((C, y_at_x0))
        else:
            x_valid, y_valid = find_closest_point_in_domain(C, x0_display, x0_display, x1_display, 
                                                            choice, ellipse_params, direction='left')
            if x_valid is not None:
                left_pairs.append((C, y_valid))

    if left_pairs:
        left_pairs.sort(key=lambda p: p[1], reverse=True)
        targets = [p[1] for p in left_pairs]
        min_dist = label_vertical_spacing * y_range
        left_positions = place_labels_safe(targets, min_dist, y_range)

        for (C, _), y_pos in zip(left_pairs, left_positions):
            ax.text(
                x0_display + left_label_x_offset*x_range,
                y_pos + left_label_y_offset*y_range,
                f'$C = {C:g}$',
                fontsize=left_label_fontsize,
                ha='right', va='center', zorder=20
            )

    # Правые подписи
    right_pairs = []
    for C in C_values:
        y_at_x1 = f(C, x1_display, choice, ellipse_params)
        
        if not np.isnan(y_at_x1):
            right_pairs.append((C, y_at_x1))
        else:
            x_valid, y_valid = find_closest_point_in_domain(C, x1_display, x0_display, x1_display, 
                                                            choice, ellipse_params, direction='right')
            if x_valid is not None:
                right_pairs.append((C, y_valid))

    if right_pairs:
        right_pairs.sort(key=lambda p: p[1], reverse=True)
        targets = [p[1] for p in right_pairs]
        right_positions = place_labels_safe(targets, min_dist, y_range)

        for (C, _), y_pos in zip(right_pairs, right_positions):
            label_text = get_function_label(C, choice)
            ax.text(
                x1_display + right_label_x_offset*x_range,
                y_pos + right_label_y_offset*y_range,
                label_text,
                fontsize=right_label_fontsize,
                ha='left', va='center', zorder=20
            )

    # Рисование пучков
    for C in C_values:
        if abs(C - a) < 1e-10:
            intersections = find_intersections(x0_display, x1_display, C_values, a, choice, ellipse_params)
            
            print(f"Найдено пересечений для C = {a}: {len(intersections)}")
            
            if intersections:
                for x_intersect in intersections:
                    y_intersect = f(a, x_intersect, choice, ellipse_params)
                    if not np.isnan(y_intersect):
                        ax.plot(x_intersect, y_intersect, 'ko', markersize=bundle_markersize, 
                               markeredgecolor='black', markerfacecolor='black', zorder=10)
                
                at_x0 = any(abs(x - x0_display) < bundle_tolerance * x_range for x in intersections)
                at_x1 = any(abs(x - x1_display) < bundle_tolerance * x_range for x in intersections)
                
                if not at_x0:
                    y_at_x0 = f(a, x0_display, choice, ellipse_params)
                    if not np.isnan(y_at_x0):
                        ax.plot(x0_display, y_at_x0, 'ko', markersize=bundle_markersize, 
                               markeredgecolor='black', markerfacecolor='black', zorder=10)
                
                if not at_x1:
                    y_at_x1 = f(a, x1_display, choice, ellipse_params)
                    if not np.isnan(y_at_x1):
                        ax.plot(x1_display, y_at_x1, 'ko', markersize=bundle_markersize, 
                               markeredgecolor='black', markerfacecolor='black', zorder=10)
            else:
                y_at_x0 = f(a, x0_display, choice, ellipse_params)
                y_at_x1 = f(a, x1_display, choice, ellipse_params)
                if not np.isnan(y_at_x0):
                    ax.plot(x0_display, y_at_x0, 'ko', markersize=bundle_markersize, 
                           markeredgecolor='black', markerfacecolor='black', zorder=10)
                if not np.isnan(y_at_x1):
                    ax.plot(x1_display, y_at_x1, 'ko', markersize=bundle_markersize, 
                           markeredgecolor='black', markerfacecolor='black', zorder=10)
            
            break

    # Заголовок
    function_name = get_function_name(choice)
    if ellipse_params:
        a_ell, b_ell = ellipse_params
        if a_ell == b_ell:
            domain_text = f"x^2 + y^2 \\leq {a_ell**2}"
        else:
            domain_text = f"\\frac{{x^2}}{{{a_ell**2}}} + \\frac{{y^2}}{{{b_ell**2}}} \\leq 1"
        ax.set_title(f'${function_name},\\quad {domain_text}$', 
                    fontsize=16, fontweight='bold', pad=20, y=1.02)
    else:
        ax.set_title(f'${function_name}$', fontsize=16, fontweight='bold', pad=20, y=1.02)

    ax.set_xlim(x_axis_start - final_x_padding*x_range, x_axis_end + final_x_padding*x_range)
    ax.set_ylim(y_axis_start - final_y_padding*y_range, y_axis_end + final_y_padding*y_range)

    ax.grid(True, alpha=0.1, linestyle='-', linewidth=0.3)

    plt.tight_layout()

    return fig, ax


def main():
    # Читаем параметры из файла
    x0, x1, a, b, ellipse_params = read_input_file()
    
    # Выбираем тип функции
    choice = choose_function()

    print("\n" + "="*50)
    print("ЗАГРУЖЕННЫЕ ПАРАМЕТРЫ")
    print("="*50)
    print(f"Границы области по x: [{x0:.3f}, {x1:.3f}]")
    print(f"Параметры для C: a = {a}, b = {b}")
    if ellipse_params:
        print(f"Параметры области: a_ell = {ellipse_params[0]}, b_ell = {ellipse_params[1]}")
    print(f"Выбрана функция: {get_function_name(choice)}")
    print("="*50)

    # Строим график
    fig, ax = draw_graphs(x0, x1, a, b, choice, ellipse_params)

    # Сохраняем
    filename = f'graph_func{choice}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.savefig(f'graph_func{choice}.pdf', bbox_inches='tight')

    plt.show()

    print(f"\nГрафики успешно построены и сохранены как graph_func{choice}.png/pdf!")


if __name__ == "__main__":
    main()