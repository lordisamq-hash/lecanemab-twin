# -*- coding: utf-8 -*-
"""
Веб-интерфейс (Streamlit) для цифрового двойника линии производства 
биосимиляра леканемаба (по патенту US 8,025,878 B2)
Исправленная версия: Фикс графика Моно + Расчет в Рублях
"""

import streamlit as st
import math
import pandas as pd

# Настройка страницы браузера
st.set_page_config(page_title="Цифровой двойник: Леканемаб", layout="wide")

def simulate_cho_batch_with_history(v_reactor, days):
    """Кинетическая модель Моно для 1 батча клеток CHO с фиксом шагов истории"""
    mu_max = 0.035        # Максимальная скорость роста клеток (1/ч)
    Ks = 0.5              # Константа Моно (г/л глюкозы)
    Y_xs = 0.25           # Выход биомассы по субстрату (г клеток / г глюкозы)
    m_s = 0.002           # Потребление глюкозы на поддержание жизни (г/г клеток * ч)
    
    # Скорректированные коэффициенты Людекинга-Пирета для выхода ~1.5 г/л
    alpha = 0.025        
    beta = 0.0030        
    
    X = 0.2               # Начальная концентрация клеток (г/л)
    S = 55.0              # Начальная концентрация глюкозы в среде DMEM (г/л)
    Product = 0.0         # Начальная концентрация антитела (г/л)
    
    dt = 0.1              # Шаг интегрирования (часы)
    total_steps = int(days * 24 / dt)
    steps_per_day = int(24 / dt)
    
    history = []
    
    for step in range(total_steps + 1):
        # Запись точек строго по целым дням без плавающей запятой
        if step % steps_per_day == 0:
            current_day = step // steps_per_day
            history.append({
                "День": current_day,
                "Клетки CHO (г/л)": round(X, 2),
                "Глюкоза (г/л)": round(S, 2),
                "Леканемаб (г/л)": round(Product, 3)
            })
            
        if S <= 0.01:
            S = 0.0
            mu = 0.0
            r_p = 0.0
            r_s = 0.0
        else:
            mu = mu_max * (S / (Ks + S))
            r_p = alpha * (mu * X) + beta * X
            r_s = (mu * X / Y_xs) + (m_s * X)
        
        X += mu * X * dt
        S -= r_s * dt
        Product += r_p * dt
        
    df_history = pd.DataFrame(history)
    return round(Product, 3), round(X, 2), round(S, 2), df_history

# Заголовок веб-интерфейса
st.title("📊 Цифровой двойник опытно-промышленной линии")
st.subheader("Технико-экономическая модель производства биосимиляра леканемаба (US 8,025,878 B2)")
st.markdown("---")

# Разделение экрана на две колонки
col_inputs, col_results = st.columns([1, 1.2])

with col_inputs:
    st.header("⚙️ Входные параметры")
    
    st.subheader("Главные параметры")
    target_pure_protein_kg = st.number_input("Целевой объем чистого белка в год, кг", min_value=0.05, max_value=10.0, value=1.5, step=0.05)
    reactor_volume_l = st.number_input("Рабочий объем биореактора CHO, л", min_value=10.0, max_value=2000.0, value=200.0, step=10.0)
    batch_days = st.number_input("Длительность цикла культивирования, дней", min_value=7, max_value=21, value=14, step=1)
    
    with st.expander("🧬 Настройка технологических потерь (Downstream)"):
        loss_centrifugation = st.slider("Потери на Центрифугировании", 0.01, 0.15, 0.05, 0.01)
        loss_chromatography = st.slider("Потери на Хроматографии Protein A", 0.05, 0.30, 0.15, 0.01)
        loss_ultrafiltration = st.slider("Потери на Ультрафильтрации", 0.01, 0.15, 0.05, 0.01)
        loss_lyophilization = st.slider("Потери при Розливе/Лиофилизации", 0.01, 0.10, 0.02, 0.01)
        
    with st.expander("💰 Настройка стоимости расходников (в Рублях, ₽)"):
        cost_media_per_l = st.number_input("Стоимость 1 л среды DMEM, ₽", value=2300.0)
        cost_protein_a_resin_l = st.number_input("Стоимость 1 л смолы Protein A, ₽", value=870000.0)
        cost_mannitol_per_kg = st.number_input("Стоимость 1 кг Маннита, ₽", value=4100.0)
        cost_vial_finish = st.number_input("Стоимость флакона + финишных работ, ₽", value=320.0)

# --- РАСЧЕТНАЯ ЧАСТЬ МОДЕЛИ ---
titer_g_l, final_cells, final_sugar, df_plots = simulate_cho_batch_with_history(reactor_volume_l, batch_days)

total_yield = (1 - loss_centrifugation) * (1 - loss_chromatography) * (1 - loss_ultrafiltration) * (1 - loss_lyophilization)

target_pure_g = target_pure_protein_kg * 1000
required_raw_protein_g = target_pure_g / total_yield
protein_per_batch_g = reactor_volume_l * titer_g_l

number_of_batches = math.ceil(required_raw_protein_g / protein_per_batch_g)
actual_raw_protein_g = number_of_batches * protein_per_batch_g
actual_pure_protein_g = actual_raw_protein_g * total_yield

dose_per_vial_g = 0.5  
number_of_vials = math.floor(actual_pure_protein_g / dose_per_vial_g) if dose_per_vial_g > 0 else 0
mannitol_per_vial_g = dose_per_vial_g * 1.2  
total_mannitol_kg = (number_of_vials * mannitol_per_vial_g) / 1000 if number_of_vials > 0 else 0

total_media_volume = number_of_batches * reactor_volume_l
media_total_cost = total_media_volume * cost_media_per_l
column_volume_l = reactor_volume_l * 0.05
resin_amortization_cost = (column_volume_l * cost_protein_a_resin_l) * (number_of_batches / 50)
mannitol_total_cost = total_mannitol_kg * cost_mannitol_per_kg
vials_total_cost = number_of_vials * cost_vial_finish

total_opex_materials = media_total_cost + resin_amortization_cost + mannitol_total_cost + vials_total_cost
cost_per_vial = total_opex_materials / number_of_vials if number_of_vials > 0 else 0.0

with col_results:
    st.header("📈 Результаты моделирования")
    
    st.subheader("Технические показатели")
    c1, c2, c3 = st.columns(3)
    c1.metric("Титр на выходе (CHO)", f"{titer_g_l} г/л")
    c2.metric("Эффективность очистки (Yield)", f"{total_yield*100:.2f}%")
    c3.metric("Циклов реактора в год", f"{number_of_batches} шт.")
    
    c4, c5, c6 = st.columns(3)
    c4.metric("Готовая продукция (Флаконы)", f"{number_of_vials:,} шт.")
    c5.metric("Расход Маннита в год", f"{total_mannitol_kg:.2f} кг")
    c6.metric("Остаточная глюкоза среды", f"{final_sugar} г/л")
    
    # --- ОТЛАЖЕННЫЙ БЛОК ГРАФИКА МОНО ---
    st.markdown("---")
    st.subheader("📉 Кинетические кривые Моно (Динамика биореактора)")
    
    # Индексируем по дням для корректного графика
    chart_data = df_plots.set_index("День")
    st.line_chart(chart_data)
    st.caption("График отображает плавный рост клеток (г/л), падение концентрации глюкозы (г/л) и наработку белка леканемаба (г/л) с шагом в 1 день.")
    
    st.markdown("---")
    st.subheader("Финансовые показатели (Материальный OPEX)")
    
    financial_data = {
        "Категория затрат": [
            "Питательные среды (DMEM + добавки)", 
            "Амортизация хроматографической смолы Protein A", 
            "Упаковка и вспомогательные вещества (Флаконы + Маннит)",
            "ИТОГО ЗАТРАТ НА СЫРЬЕ В ГОД"
        ],
        "Сумма в год (₽)": [
            f"{media_total_cost:,.2f} ₽",
            f"{resin_amortization_cost:,.2f} ₽",
            f"{mannitol_total_cost + vials_total_cost:,.2f} ₽",
            f"{total_opex_materials:,.2f} ₽"
        ]
    }
    st.table(pd.DataFrame(financial_data))
    
    st.success(f"### 🎯 Расчетная себестоимость одного флакона биосимиляра: {cost_per_vial:,.2f} ₽")
