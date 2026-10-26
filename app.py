import streamlit as st
import math
import pandas as pd

st.set_page_config(page_title="Цифровой двойник: леканемаб", layout="wide")


def simulate_fed_batch(
    days, mu_max, Ks, Kd, q_p,
    reactor_volume_l,
    scale_loss_per_decade_pct=10.0,
    S_feed=200.0, feed_rate=0.172,
    S0=6.0, X0=0.4, P0=0.0,
    t_feed_start=48.0,
    dt=0.1,
):
    """Fed-Batch с ростом объёма.
    Продукция: Luedeking-Piret при α = 0 (β = qp).
    Гибель: Kd первого порядка в уравнении dX/dt.
    """
    Xmax = 16.0
    Y_xs = 0.45
    m_s = 0.001

    F_vol = feed_rate / S_feed if S_feed > 0 else 0.0

    V, X, S, P = 1.0, X0, S0, P0
    total_hours = days * 24
    history = []

    def deriv(X, S, P, V, t):
        F = F_vol if t >= t_feed_start else 0.0
        D = F / V
        if S <= 0.0:
            return (-Kd - D) * X, D * S_feed, -D * P, F
        mu = mu_max * (S / (Ks + S)) * max(0.0, 1.0 - X / Xmax)
        dX = (mu - Kd - D) * X
        dS = -(mu * X) / Y_xs - m_s * X + D * (S_feed - S)
        dP = q_p * X - D * P           # Luedeking-Piret, α = 0, β = qp
        return dX, dS, dP, F

    n_steps = int(total_hours / dt) + 1
    for step in range(n_steps):
        t = step * dt
        mu_now = (mu_max * S / (Ks + S) * max(0.0, 1.0 - X / Xmax)) if S > 0 else 0.0

        if step % int(round(24 / dt)) == 0:
            history.append({
                "Сутки (t)": int(round(t / 24)),
                "Клетки CHO, ×10⁹ кл/л (X)": round(X, 3),
                "Глюкоза, г/л (S)": round(max(S, 0.0), 3),
                "Леканемаб, г/л (P)": round(P, 4),
                "Скорость роста, ч⁻¹ (μ)": round(mu_now, 5),
                "Кратность объёма (V/V₀)": round(V, 4),
            })

        k1 = deriv(X, S, P, V, t)
        k2 = deriv(X + 0.5*dt*k1[0], S + 0.5*dt*k1[1],
                   P + 0.5*dt*k1[2], V + 0.5*dt*k1[3], t + 0.5*dt)
        k3 = deriv(X + 0.5*dt*k2[0], S + 0.5*dt*k2[1],
                   P + 0.5*dt*k2[2], V + 0.5*dt*k2[3], t + 0.5*dt)
        k4 = deriv(X + dt*k3[0], S + dt*k3[1],
                   P + dt*k3[2], V + dt*k3[3], t + dt)
        for i, (a, b, c, d) in enumerate(zip(k1, k2, k3, k4)):
            inc = (dt/6.0) * (a + 2*b + 2*c + d)
            if i == 0: X += inc
            elif i == 1: S += inc
            elif i == 2: P += inc
            else: V += inc

        X, S, P = max(X, 0.0), max(S, 0.0), max(P, 0.0)

    if reactor_volume_l > 10.0:
        decades = math.log10(reactor_volume_l / 10.0)
        loss = (scale_loss_per_decade_pct / 100.0) * decades
        P_s = P * max(0.0, 1.0 - loss)
    else:
        P_s = P

    return {
        "titer": P_s,
        "titer_raw": P,
        "mass_g": P_s * V * reactor_volume_l,
        "V_final_l": V * reactor_volume_l,
        "X": X, "S": S, "V": V,
        "loss_pct": (1.0 - P_s / P) * 100.0 if P > 0 else 0.0,
        "F_vol": F_vol,
        "history": pd.DataFrame(history),
    }


def material_balance_chain(mass_raw_g, vial_mg=500.0):
    stages = [
        ("Осветление и глубинная фильтрация", 0.95),
        ("Аффинная хроматография (Protein A)", 0.85),
        ("Тангенциальная ультрафильтрация и полировка", 0.95),
        ("Формуляция, розлив и лиофилизация", 0.98),
    ]
    rows = [{"Стадия": "Вход: неочищенный белок",
             "Выход стадии, %": 100.0,
             "Масса, г": round(mass_raw_g, 2)}]
    current = mass_raw_g
    for name, y in stages:
        current *= y
        rows.append({"Стадия": name,
                     "Выход стадии, %": round(y * 100, 1),
                     "Масса, г": round(current, 2)})
    n_vials = int(current / (vial_mg / 1000.0))
    return pd.DataFrame(rows), current, n_vials


# ============================================================
# ИНТЕРФЕЙС
# ============================================================
st.title("Цифровой двойник опытно-промышленной линии")
st.subheader("Моделирование кинетики Fed-Batch с учётом масштабирования")

col_in, col_out = st.columns([1, 1.3])

with col_in:
    st.markdown("### Параметры масштаба")
    V_l = st.number_input("Объём биореактора, л (V)",
                          min_value=10.0, max_value=2000.0,
                          value=200.0, step=10.0)
    days = st.number_input("Длительность культивирования, сутки (t)",
                           min_value=1, max_value=14, value=10, step=1)
    loss_pct = st.slider(
        "Потеря титра при масштабировании, %/декада (k)",
        min_value=0.0, max_value=30.0, value=10.0, step=0.5,
        help="Принято 10 % как консервативное значение в литературном "
             "диапазоне 5–15 % [Xing et al., 2009]. "
             "При V = 200 л фактическая потеря ≈ 13 %."
    )

    st.markdown("### Параметры подпитки")
    S_feed = st.number_input("Концентрация субстрата в подпитке, г/л (S_feed)",
                             min_value=50.0, max_value=1000.0,
                             value=200.0, step=10.0)
    feed_rate = st.number_input(
        "Скорость подпитки по глюкозе, г/(л·ч) (feed_rate)",
        min_value=0.0, max_value=2.0, value=0.172,
        format="%.3f", step=0.001,
        help="Подача глюкозы на 1 л рабочего объёма. "
             "Документ: 0,15–0,172 г/(л·ч)."
    )

    st.markdown("### Кинетика штамма CHO")
    mu = st.number_input("Максимальная удельная скорость роста, 1/ч (μmax)",
                         min_value=0.001, max_value=0.1, value=0.040,
                         format="%.3f", step=0.001,
                         help="Уравнение Моно: μ = μmax·S/(Ks+S)·(1–X/Xmax).")
    Ks = st.number_input("Константа насыщения по субстрату, г/л (Ks)",
                         min_value=0.01, max_value=5.0, value=0.50,
                         format="%.2f", step=0.01)
    Kd = st.number_input(
        "Константа гибели клеток, 1/ч (Kd)",
        min_value=0.000, max_value=0.05, value=0.004,
        format="%.3f", step=0.001,
        help="Спонтанный апоптоз: dX/dt = (μ – Kd – D)·X. "
             "Отдельное уравнение для мёртвых клеток Xd не вводится: "
             "при Kd = 0,004 ч⁻¹ и 240 ч вклад в баланс < 5 %."
    )
    qp = st.number_input(
        "Удельная продуктивность, г/(10⁹ кл·ч) (qp)",
        min_value=0.0001, max_value=0.01, value=0.0012,
        format="%.4f", step=0.0001,
        help="Luedeking-Piret при α = 0: dP/dt = qp·X – D·P. "
             "Синтез не ассоциирован с ростом (стандарт для IgG1-CHO)."
    )

res = simulate_fed_batch(
    days=int(days), mu_max=mu, Ks=Ks, Kd=Kd, q_p=qp,
    reactor_volume_l=V_l,
    scale_loss_per_decade_pct=loss_pct,
    S_feed=S_feed, feed_rate=feed_rate,
)

with col_out:
    st.subheader("Кинетические кривые")
    st.caption(
        f"Удельный объёмный поток подпитки: "
        f"F/V = {res['F_vol']:.5f} л/(л·ч) "
        f"(= feed_rate / S_feed = {feed_rate:.3f} / {S_feed:.0f}).  "
        f"Модель: Моно + Fed-Batch + Kd + Luedeking-Piret (α = 0)."
    )
    st.line_chart(
        res["history"].set_index("Сутки (t)")[
            ["Клетки CHO, ×10⁹ кл/л (X)",
             "Глюкоза, г/л (S)",
             "Леканемаб, г/л (P)"]]
    )

    st.subheader("Результаты серии")
    st.table(pd.DataFrame({
        "Параметр": ["Титр без масштаба", "Потеря при масштабировании",
                     "Титр с учётом масштаба",
                     "Объём в конце культивирования",
                     "Масса неочищенного белка за серию",
                     "X финальная", "S финальная",
                     "Kd (коэффициент гибели)"],
        "Обозначение": ["P", "ΔP", "P_scaled",
                        "V_final", "m", "X_final", "S_final", "Kd"],
        "Значение": [f"{res['titer_raw']:.3f}",
                     f"{res['loss_pct']:.1f}",
                     f"{res['titer']:.3f}",
                     f"{res['V_final_l']:.1f}",
                     f"{res['mass_g']:.1f}",
                     f"{res['X']:.2f}",
                     f"{res['S']:.2f}",
                     f"{Kd:.3f}"],
        "Ед. изм.": ["г/л", "%", "г/л", "л", "г",
                     "×10⁹ кл/л", "г/л", "ч⁻¹"],
    }))
    st.caption(
        "Kd входит в уравнение dX/dt = (μ – Kd – D)·X. "
        "Отдельное уравнение для мёртвых клеток Xd не вводится — "
        "см. раздел 2.2 главы 2 «Что не включено в модель и почему»."
    )

# ============================================================
# МАТЕРИАЛЬНЫЙ БАЛАНС
# ============================================================
st.subheader("Материальный баланс серии")
st.caption(
    "Цепочка выходов по стадиям Downstream. Вход — масса неочищенного "
    "белка из результатов серии выше."
)
balance_df, final_mass, n_vials = material_balance_chain(res["mass_g"])
st.table(balance_df)

# ============================================================
# ИТОГО
# ============================================================
st.subheader("Итого")
k_isp = final_mass / res["mass_g"] if res["mass_g"] > 0 else 0.0
k_rash = res["mass_g"] / final_mass if final_mass > 0 else 0.0

col1, col2, col3 = st.columns(3)
col1.metric("Масса неочищенного белка", f"{res['mass_g']:.1f} г")
col2.metric("Масса готового леканемаба", f"{final_mass:.1f} г")
col3.metric("Готовых флаконов (500 мг)", f"{n_vials}")

col1, col2, col3 = st.columns(3)
col1.metric("K_исп", f"{k_isp:.3f}")
col2.metric("K_расх", f"{k_rash:.2f}")
col3.metric("Объём в конце, л", f"{res['V_final_l']:.1f}")
