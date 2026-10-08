import streamlit as st
import math
import pandas as pd

st.set_page_config(page_title="Масштабирование Моно: Леканемаб", layout="wide")

def simulate_fed_batch_mono(days, mu_max, Ks, Kd, q_p, scale_drop_factor, reactor_volume_l):
    """Кинетическая модель Моно (Рунге-Кутта 4-го порядка) с учетом масштабного фактора"""
    Xmax = 16.0            # Максимальная емкость среды, x10^9 кл/л
    Y_xs = 0.45            # Выход биомассы по субстрату, x10^9 кл/г
    m_s = 0.001            # Расход на поддержание жизни, г/(10^9 кл * ч)
    F_feed_rate = 0.15     # Скорость подачи подпитки, г/(л * ч)
    
    # Стартовые удельные концентрации на 1 литр
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
        
    # --- ИНЖЕНЕРНЫЙ ЭФФЕКТ МАСШТАБИРОВАНИЯ ---
    if reactor_volume_l > 10.0:
        scale_loss = (scale_drop_factor / 100.0) * math.log10(reactor_volume_l / 10.0)
        P_scaled = P * (1.0 - scale_loss)
    else:
        P_scaled = P
        
    return round(P_scaled, 3), round(X, 2), round(S, 2), pd.DataFrame(history)

# --- ИНТЕРФЕЙС STREAMLIT ---
st.title(" Цифровой двойник опытно-промышленной линии")
st.subheader("Моделирование кинетики Fed-Batch биосинтеза с учетом масштабирования")
st.markdown("---")

col_inputs, col_graph = st.columns([1, 1.3])

with col_inputs:
    st.markdown("###  Масштаб аппарата")
    reactor_volume_l = st.number_input("Объем серии (биореактора), л", min_value=10.0, max_value=2000.0, value=200.0, step=10.0)
    batch_days = st.number_input("Дни культивирования (модель Моно: 10 дней)", min_value=1, max_value=30, value=10, step=1)
    scale_drop_factor = st.slider("Стресс-фактор масштаба, % падения титра (ориентир: 5-15%)", 0, 30, 10)
    
    st.markdown("###  Биокинетика штамма CHO (Глава 2)")
    mu_max_val = st.number_input("Скорость роста μmax, 1/ч (документ: 0.04)", value=0.040, format="%.3f")
    Ks_val = st.number_input("Константа насыщения Ks, г/л (документ: 0.5)", value=0.50, format="%.2f")
    Kd_val = st.number_input("Скорость гибели Kd, 1/ч (документ: 0.004)", value=0.004, format="%.3f")
    q_p_val = st.number_input("Синтез mAb qp, г/(10⁹кл*ч) (документ: 0.0012)", value=0.0012, format="%.4f")

# Вычисление модели
titer_g_l, final_cells, final_sugar, df_plots = simulate_fed_batch_mono(batch_days, mu_max_val, Ks_val, Kd_val, q_p_val, scale_drop_factor, reactor_volume_l)

# Расчет абсолютной массы сырого белка на батч
protein_per_batch_g = reactor_volume_l * titer_g_l

with col_graph:
    st.subheader(" Кинетические кривые (Динамика биореактора)")
    chart_data = df_plots.set_index("День")
    st.line_chart(chart_data)
    
    st.markdown("### Результаты серии")
    metrics_data = {
        "Параметр материального потока": [
            "Удельный выход белка из реактора (Титр с учетом масштаба)",
            f"Суммарная масса сырого белка на серию (для реактора {reactor_volume_l} л)",
            "Финальная плотность жизнеспособных клеток CHO",
            "Остаточная глюкоза в культуральной среде"
        ],
        "Значение": [
            f"{titer_g_l:.3f} г/л",
            f"{protein_per_batch_g:.2f} г",
            f"{final_cells} × 10⁹ кл/л",
            f"{final_sugar} г/л"
        ]
    }
    st.table(pd.DataFrame(metrics_data))
