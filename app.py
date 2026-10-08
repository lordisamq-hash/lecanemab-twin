import matplotlib.pyplot as plt
import streamlit as st

# Заголовок в интерфейсе Streamlit
st.header("📐 Технологическая схема производства по ОСТ 64-02-003-2002")

# Настройка шрифта для ГОСТ
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]

# Данные для материального баланса
stages = [
    {
        "id": "ВСП 1",
        "name": "Приготовление питательных сред\nи буферных растворов",
        "desc": "Подготовка WFI, сред, буферов, маннита",
    },
    {
        "id": "ТП 1",
        "name": "Подготовка и масштабирование\nинокулята клеток CHO",
        "desc": "Криопробирка WCB -> шейкеры -> сид-реактор 20 л",
    },
    {
        "id": "ТП 2",
        "name": "Доливно-периодический биосинтез\n(Fed-Batch)",
        "desc": "Биореактор Single-Use 200 л\nТитр: 1,392 г/л | Сырой белок: 278,40 г",
    },
    {
        "id": "ТП 3",
        "name": "Первичное осветление и разделение",
        "desc": "Глубинная фильтрация / Центрифугирование\nПотери: 5,0% | Выход: 264,48 г",
    },
    {
        "id": "ТП 4",
        "name": "Аффинная хроматография\n(Выделение на Protein A)",
        "desc": "Удаление HCP/HCD | Инактивация вирусов (pH 3,5)\nПотери: 15,0% | Выход: 224,81 г",
    },
    {
        "id": "ТП 5",
        "name": "Тангенциальная ультрафильтрация\nи полировка",
        "desc": "Концентрирование и смена буфера\nПотери: 5,0% | Выход: 213,57 г",
    },
    {
        "id": "ТП 6",
        "name": "Вирусная нанофильтрация\nи стерилизующий розлив",
        "desc": "Фильтр 20 нм | Розлив (Класс А)\nВнесение маннита (массовая доля 1,2 : 1)",
    },
    {
        "id": "УМО 7",
        "name": "Сублимационная сушка\nи укупорка",
        "desc": "Лиофилизация под вакуумом (флаконы 10R)\nПотери: 2,0% | Чистый белок: 209,24 г",
    },
    {
        "id": "ГП",
        "name": "ГОТОВЫЙ ПРОДУКТ\n(Лиофилизат леканемаба)",
        "desc": "Выход серии: ~418 флаконов по 500 мг\nСквозной Yield очистки: 75,16%",
    },
]

# Создаем фигуру под пропорции страницы А4 (вертикальная ориентация)
fig, ax = plt.subplots(figsize=(10, 14))

# Начальные координаты для отрисовки блоков сверху вниз
start_y = 100
box_height = 6
box_width = 55
box_x = 22  # Центрирование блоков по оси X

y_coords = []

# Отрисовка прямоугольников стадий (по ОСТ 64-02-003-2002)
for i, stage in enumerate(stages):
    current_y = start_y - (i * 11)
    y_coords.append(current_y)

    # Строгий черный контур, белый фон для чертежного стиля
    rect = plt.Rectangle(
        (box_x, current_y),
        box_width,
        box_height,
        facecolor="white",
        edgecolor="black",
        linewidth=1.5,
    )
    ax.add_patch(rect)

    # Индекс стадии (ВСП 1, ТП 2 и т.д.) в левом верхнем углу блока
    ax.text(
        box_x + 1,
        current_y + box_height - 1.5,
        stage["id"],
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="top",
    )

    # Основное название стадии
    ax.text(
        box_x + box_width / 2,
        current_y + box_height / 2 + 0.8,
        stage["name"],
        fontsize=11,
        fontweight="semibold",
        ha="center",
        va="center",
    )

    # Технологическое описание и параметры материального баланса
    ax.text(
        box_x + box_width / 2,
        current_y + 1.2,
        stage["desc"],
        fontsize=9.5,
        fontstyle="italic",
        ha="center",
        va="center",
    )

# Отрисовка линий потоков и точек КТП
for i in range(len(stages) - 1):
    y_top = y_coords[i]
    y_bottom = y_coords[i + 1] + box_height
    x_center = box_x + box_width / 2

    # Линии основного потока (жирные начиная со стадии ТП 2)
    lw = 2.5 if i >= 1 else 1.2
    ax.annotate(
        "",
        xy=(x_center, y_bottom),
        xytext=(x_center, y_top),
        arrowprops=dict(
            arrowstyle="->", color="black", linewidth=lw, shrinkA=0, shrinkB=0
        ),
    )

    # Исправлено: Отрисовка боковых линий отходов/потерь для стадий Downstream (ТП 3, ТП 4, ТП 5, УМО 7)
    if i in:
        y_arrow = y_top - 2.5
        ax.annotate(
            "",
            xy=(box_x + box_width + 8, y_arrow - 2),
            xytext=(box_x + box_width, y_arrow),
            arrowprops=dict(
                arrowstyle="->",
                color="black",
                linewidth=1.0,
                connectionstyle="angle,angleA=0,angleB=-90,rad=0",
            ),
        )
        ax.text(
            box_x + box_width + 1,
            y_arrow + 0.5,
            "Потери / Отход",
            fontsize=8.5,
            ha="left",
        )

# Добавление точек КТП
ktp_points = [
    {"idx": 1, "label": "КТП 1.1\nЖизнеспособность"},
    {"idx": 2, "label": "КТП 2.1\nТитр (1,392 г/л)"},
    {"idx": 7, "label": "КТП 7.1\nВыходной контроль\n(SE-HPLC >= 95%)"},
]

for ktp in ktp_points:
    y_p = y_coords[ktp["idx"]] + box_height / 2
    x_p = box_x + box_width

    circle = plt.Circle(
        (x_p, y_p), 1.0, facecolor="white", edgecolor="black", linewidth=1.2, zorder=5
    )
    ax.add_patch(circle)

    ax.text(
        x_p + 2,
        y_p,
        ktp["label"],
        fontsize=9,
        va="center",
        ha="left",
        bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.8),
    )

# Настройки отображения графика
ax.set_xlim(0, 100)
ax.set_ylim(0, 110)
ax.axis("off")

# Вывод графика в интерфейс Streamlit
st.pyplot(fig)

# Подрисуночная подпись текстом Streamlit
st.caption(
    "**Рисунок 3.1** – Технологическая схема опытно-промышленного производства биосимиляра леканемаба (ОПР-200) по ОСТ 64-02-003-2002"
)
