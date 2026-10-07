# -*- coding: utf-8 -*-
"""
Веб-интерфейс (Streamlit) для цифрового двойника линии производства леканемаба
Финальная инженерная версия: Только кинетика Моно, Масштабирование и Материальный баланс (Без экономики)
"""

import streamlit as st
import math
import pandas as pd

st.set_page_config(page_title="Цифровой двойник: Леканемаб", layout="wide")

def simulate_fed_batch_mono(days, mu_max, Ks, Kd, q_p, scale_drop_factor, reactor_volume_l):
    """Кинетическая модель Моно (Рунге-Кутта 4-го порядка) с учетом масштабного фактора"""
    Xmax = 16.0            
    Y_xs = 0.45            
    m_s = 0.001            
    F_feed_rate = 0.15     
    
    X = 0.4                
    S = 6.0                
    P = 0.0                
    dt = 0.1               
    total_hours = days * 24      
    
    history = []
    
    def equations(X_val, S_val, P_val, current_t):
        if current_t >= 48.0 and S_val < 4.0:
            F_feed = F_feed_rate
            dVdt = 0.0125  
        else:
            F_feed = 0.0
            dVdt = 0.0
            
        if S_val <= 0.001:
            S_val = 0.0
            mu = 0.0
            dXdt = -Kd * X_val
            dSdt = 0.0
            dPdt = 0.0
        else:
            mu = mu_max * (S_val / (Ks + S_val)) * (1.0 - (X_val / Xmax))
            dXdt = (mu - Kd) * X_val
            dSdt = - (1.0 / Y_xs) * mu * X_val - m_s * X_val + F_feed - (S_val * dVdt / 1.0)
            dPdt = q_p * X_val
            
        return dXdt, dSdt, dPdt

    for step in range(int(total_hours / dt) + 1):
        t = step * dt
        if step % int(24 / dt) == 0:
            history.append({
                "День": int(t // 24),
                "Клетки X (x10⁹ кл/л)": round(X, 2),
                "Глюкоза S (г/л)": round(S, 2),
                "Леканемаб P (г/л)": round(P, 3)
            })
            
        k1_x, k1_s, k1_p = equations(X, S, P, t)
        k2_x, k2_s, k2_p = equations(X + 0.5*dt*k1_x, S + 0.5*dt*k1_s, P + 0.5*dt*k1_p, t + 0.5*dt)
        k3_x, k3_s, k3_p = equations(X + 0.5*dt*k2_x, S + 0.5*dt*k2_s, P + 0.5*dt*k2_p, t + 0.5*dt)
        k4_x, k4_s, k4_p = equations(X + dt*k3_x, S + dt*k3_s, P + dt*k3_p, t + dt)
        
        X += (dt / 6.0) * (k1_x + 2*k2_x + 2*k3_x + k4_x)
        S += (dt / 6.0) * (k1_s + 2*k2_s + 2*k3_s + k4_s)
        P += (dt / 6.0) * (k1_p + 2*k2_p + 2*k3_p + k4_p)
        
    if reactor_volume_l > 10.0:
        scale_loss = (scale_drop_factor / 100.0) * math.log10(reactor_volume_l / 10.0)
        P_scaled = P * (1.0 - scale_loss)
    else:
        P_scaled = P
        
    return round(P_scaled, 3), round(X, 2), round(S, 2), pd.DataFrame(history)

st.title("📊 Цифровой двойник опытно-промышленной линии")
st.subheader("Технологический регламент материального баланса по стандарту ОСТ 64-02-003-2002")
st.markdown("---")

# Теперь делим экран на 3 колонки (без ценников)
col_main, col_mono, col_losses = st.columns(3)

with col_main:
    st.markdown("### 🏢 Масштаб завода")
    target_pure_protein_kg = st.number_input("Годовой план по чистому белку, кг", min_value=0.01, max_value=50.0, value=1.5, step=0.05)
    reactor_volume_l = st.number_input("Объем серии (биореактора), л", min_value=10.0, max_value=2000.0, value=200.0, step=10.0)
    batch_days = st.number_input("Дни культивирования", min_value=1, max_value=30, value=10, step=1)
    scale_drop_factor = st.slider("Стресс-фактор масштаба, % падения титра (рекомендовано: 5-15%)", 0, 30, 10)

with col_mono:
    st.markdown("### 🧬 Кинетика Моно (Глава 2)")
    mu_max = st.number_input("Скорость роста μmax, 1/ч (документ: 0.04)", value=0.040, format="%.3f")
    Ks = st.number_input("Константа Насыщения Ks, г/л (документ: 0.5)", value=0.50, format="%.2f")
    Kd = st.number_input("Скорость гибели Kd, 1/ч (документ: 0.004)", value=0.004, format="%.3f")
    q_p = st.number_input("Синтез mAb qp, г/(10⁹кл*ч) (документ: 0.0012)", value=0.0012, format="%.4f")

with col_losses:
    st.markdown("### 🧪 Потери очистки (Downstream)")
    loss_centrifugation = st.slider("Потери ТП 3 (Осветление)", 0.01, 0.15, 0.05, 0.01)
    loss_chromatography = st.slider("Потери ТП 4 (Protein A)", 0.05, 0.30, 0.15, 0.01)
    loss_ultrafiltration = st.slider("Потери ТП 5 (Ультрафильтр)", 0.01, 0.15, 0.05, 0.01)
    loss_lyophilization = st.slider("Потери УМО 7 (Лиофилизация)", 0.01, 0.10, 0.02, 0.01)

st.markdown("---")

# --- ВЫЧИСЛЕНИЯ МОДЕЛИ ---
titer_g_l, final_cells, final_sugar, df_plots = simulate_fed_batch_mono(batch_days, mu_max, Ks, Kd, q_p, scale_drop_factor, reactor_volume_l)

total_yield = (1 - loss_centrifugation) * (1 - loss_chromatography) * (1 - loss_ultrafiltration) * (1 - loss_lyophilization)
protein_per_batch_g = reactor_volume_l * titer_g_l
target_pure_g = target_pure_protein_kg * 1000

number_of_batches = math.ceil(target_pure_g / (protein_per_batch_g * total_yield)) if protein_per_batch_g > 0 else 0
actual_raw_protein_g = number_of_batches * protein_per_batch_g
actual_pure_protein_g = actual_raw_protein_g * total_yield

dose_per_vial_g = 0.5  
number_of_vials = math.floor(actual_pure_protein_g / dose_per_vial_g) if dose_per_vial_g > 0 else 0
total_mannitol_kg = (number_of_vials * (dose_per_vial_g * 1.2)) / 1000

# --- ИНЖЕНЕРНЫЙ ОТЧЕТ И ГРАФИКИ ---
col_graph, col_tables = st.columns(2)

with col_graph:
    st.subheader("📉 Кинетика Fed-Batch цикла (Модель Моно)")
    chart_data = df_plots.set_index("День")
    st.line_chart(chart_data)
    st.caption("Удельные кривые биосинтеза на 1 л объема по методу Рунге-Кутты 4-го порядка.")

with col_tables:
    st.subheader("📋 Сводный материальный баланс серии (ОСТ 64-02-003-2002)")
    metrics_data = {
        "Параметр материального баланса": [
            "Выход белка из реактора (Титр с учетом масштаба)",
            f"Фактический сырой белок за 1 цикл ({reactor_volume_l} л)",
            "Общая технологическая эффективность очистки (Yield)",
            "Необходимое число циклов (батчей) в год",
            "Всего выпущено флаконов готового продукта (доза 500 мг)",
            "Годовой расход Маннита (патентная пропорция 1.2:1)"
        ],
        "Значение": [
            f"{titer_g_l:.3f} г/л",
            f"{protein_per_batch_g:.2f} г",
            f"{total_yield*100:.2f}%",
            f"{number_of_batches} батчей/год",
            f"{number_of_vials:,} шт.",
            f"{total_mannitol_kg:.2f} кг"
        ]
    }
    st.table(pd.DataFrame(metrics_data))
