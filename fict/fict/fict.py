def transform_fict(input_filename, output_filename):
    with open(input_filename, 'r') as f:
        content = f.read().strip()
    
    # Разбиваем на числа
    numbers = list(map(int, content.split()))
    
    # Проверяем, что количество чисел кратно 5
    if len(numbers) % 5 != 0:
        raise ValueError("Количество чисел в файле не кратно 5")
    
    # Читаем области
    old_regions = []
    for i in range(0, len(numbers), 5):
        region = numbers[i:i+5]
        old_regions.append(region)
    
    # Создаем новые области
    new_regions = []
    
    for region in old_regions:
        type_val, xa, xb, ya, yb = region
        
        # Новые координаты границ
        xa_new = xa * 2
        xb_new = xb * 2
        ya_new = ya * 2
        yb_new = yb * 2
        
        # Создаем 4 подобласти
        # Левая нижняя
        new_regions.append([type_val, xa_new, xa_new+1, ya_new, ya_new+1])
        # Правая нижняя
        new_regions.append([type_val, xa_new+1, xb_new, ya_new, ya_new+1])
        # Левая верхняя
        new_regions.append([type_val, xa_new, xa_new+1, ya_new+1, yb_new])
        # Правая верхняя
        new_regions.append([type_val, xa_new+1, xb_new, ya_new+1, yb_new])
    
    # Записываем результат
    with open(output_filename, 'w') as f:
        for region in new_regions:
            f.write(' '.join(map(str, region)) + ' ')

# Использование
transform_fict('fict.txt', 'fict_new.txt')

def transform_kray(input_filename, output_filename):
    with open(input_filename, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    
    for line in lines:
        # Пропускаем пустые строки
        if not line.strip():
            new_lines.append('\n')
            continue
            
        # Разбиваем строку на числа
        numbers = list(map(int, line.strip().split()))
        
        if len(numbers) != 6:
            print(f"Предупреждение: строка '{line.strip()}' имеет не 6 чисел, пропускаем")
            continue
        
        # Первые два числа - идентификаторы краевого условия
        type1, type2 = numbers[0], numbers[1]
        # Остальные 4 числа - координаты (xa, xb, ya, yb)
        xa, xb, ya, yb = numbers[2], numbers[3], numbers[4], numbers[5]
        
        # Новые координаты границ (умножаем на 2)
        xa_new = xa * 2
        xb_new = xb * 2
        ya_new = ya * 2
        yb_new = yb * 2
        
        # Обрабатываем в зависимости от типа линии
        if xa == xb:  # Вертикальная линия
            # Создаем 2 подобласти (по вертикали)
            new_lines.append(f"{type1} {type2} {xa_new} {xa_new+1} {ya_new} {ya_new+1}\n")
            new_lines.append(f"{type1} {type2} {xa_new} {xa_new+1} {ya_new+1} {yb_new}\n")
            
        elif ya == yb:  # Горизонтальная линия
            # Создаем 2 подобласти (по горизонтали)
            new_lines.append(f"{type1} {type2} {xa_new} {xa_new+1} {ya_new} {ya_new+1}\n")
            new_lines.append(f"{type1} {type2} {xa_new+1} {xb_new} {ya_new} {ya_new+1}\n")
            
        else:  # Прямоугольная область (хотя в краевых обычно линии)
            # Создаем 4 подобласти
            new_lines.append(f"{type1} {type2} {xa_new} {xa_new+1} {ya_new} {ya_new+1}\n")
            new_lines.append(f"{type1} {type2} {xa_new+1} {xb_new} {ya_new} {ya_new+1}\n")
            new_lines.append(f"{type1} {type2} {xa_new} {xa_new+1} {ya_new+1} {yb_new}\n")
            new_lines.append(f"{type1} {type2} {xa_new+1} {xb_new} {ya_new+1} {yb_new}\n")
    
    # Записываем результат
    with open(output_filename, 'w') as f:
        f.writelines(new_lines)

# Использование
transform_kray('kray.txt', 'kray_new.txt')