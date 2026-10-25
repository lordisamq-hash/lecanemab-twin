import streamlit as st
import math
import pandas as pd

st.set_page_config(page_title="Цифровой двойник: леканемаб", layout="wide")


def simulate_fed_batch(
    days, mu_max, Ks, Kd, q_p,
    reactor_volume_l,
    scale_loss_per_decade_pct=10.0,
    S_feed=200.0, F_vol=0.00075,
    S0=6.0, X0=0.4, P0=0.0,
    t_feed_start=48.0,
    dt=0.1, method="RK4",
):
    """Fed-Batch с ростом объёма. Концентрации — на 1 л начального объёма."""
    Xmax = 16.0
    Y_xs = 0.45
    m_s = 0.001

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
        dP = q_p * X - D * P
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

        if method == "RK4":
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
        else:
            dX, dS, dP, F = deriv(X, S, P, V, t)
            X += dt*dX; S += dt*dS; P += dt*dP; V += dt*F

        X, S, P = max(X, 0.0), max(S, 0.0), max(P, 0.0)

    if reactor_volume_l > 10.0:
        decades = math.log10(reactor_volume_l / 10.0)
        loss = (scale_loss_per_decade_pct / 100.0) * decades
        P_s = P * max(0.0, 1.0 - loss)
    else:
        P_s = P

    return {"titer": P_s, "titer_raw": P,
            "mass_g": P_s * V * reactor_volume_l,
            "X": X, "S": S, "V": V, "history": pd.DataFrame(history)}


@st.cache_data(show_spinner=False)
def compare_methods(_days, _mu, _Ks, _Kd, _qp, _V, _Sf, _Fv):
    ref_P = simulate_fed_batch(
        days=_days, mu_max=_mu, Ks=_Ks, Kd=_Kd, q_p=_qp,
        reactor_volume_l=_V, S_feed=_Sf, F_vol=_Fv,
        dt=0.01, method="RK4",
    )["titer_raw"]

    rows = []
    for method, dt in [("RK4", 0.01), ("RK4", 0.1), ("RK4", 0.2),
                       ("Euler", 0.01), ("Euler", 0.1),
                       ("Euler", 0.5), ("Euler", 1.0)]:
        r = simulate_fed_batch(
            days=_days, mu_max=_mu, Ks=_Ks, Kd=_Kd, q_p=_qp,
            reactor_volume_l=_V, S_feed=_Sf, F_vol=_Fv,
            dt=dt, method=method,
        )
        rows.append({
            "Метод интегрирования": method,
            "Шаг интегрирования, ч (Δt)": dt,
            "Титр без масштаба, г/л (P)": round(r["titer_raw"], 5),
            "Клетки, ×10⁹ кл/л (X)": round(r["X"], 4),
            "Отклонение от RK4 (Δt=0,01), %":
                round((r["titer_raw"] - ref_P) / ref_P * 100.0, 3),
        })
    return pd.DataFrame(rows)


@st.cache_data(show_spinner=False)
def compare_trajectories(_days, _mu, _Ks, _Kd, _qp, _V, _Sf, _Fv, _dt_cmp):
    """Две траектории (RK4 и Эйлер) при одинаковом шаге Δt."""
    kwargs = dict(
        days=_days, mu_max=_mu, Ks=_Ks, Kd=_Kd, q_p=_qp,
        reactor_volume_l=_V, S_feed=_Sf, F_vol=_Fv, dt=_dt_cmp,
    )
    r_rk4 = simulate_fed_batch(method="RK4", **kwargs)
    r_eul = simulate_fed_batch(method="Euler", **kwargs)

    df = pd.DataFrame({"Сутки (t)": r_rk4["history"]["Сутки (t)"]})
    for col in ["Клетки CHO, ×10⁹ кл/л (X)",
                "Глюкоза, г/л (S)",
                "Леканемаб, г/л (P)"]:
        df[f"{col} · RK4"] = r_rk4["history"][col].values
        df[f"{col} · Эйлер"] = r_eul["history"][col].values

    diffs = {
        "ΔP, %": (r_eul["titer_raw"] - r_rk4["titer_raw"])
                  / r_rk4["titer_raw"] * 100.0,
        "ΔX, %": (r_eul["X"] - r_rk4["X"]) / r_rk4["X"] * 100.0,
        "ΔS, г/л": r_eul["S"] - r_rk4["S"],
    }
    return df, diffs


# --- ИНТЕРФЕЙС ---
st.title("Цифровой двойник опытно-промышленной линии")
st.subheader("Моделирование кинетики Fed-Batch с учётом масштабирования")

col_in, col_out = st.columns([1, 1.3])

with col_in:
    st.markdown("### Параметры масштаба")
    V_l = st.number_input("Объём биореактора, л (V)",
                          min_value=10.0, max_value=2000.0,
                          value=200.0, step=10.0)
    days = st.number_input("Длительность культивирования, сутки (t)",
                           min_value=1, max_value=30, value=10, step=1)
    loss_pct = st.slider("Потеря титра при масштабировании, %/декада (scale_loss)",
                         min_value=0.0, max_value=30.0, value=10.0, step=0.5,
                         help="При V = 200 л и 10%/декада фактическая потеря ≈ 13%.")

    st.markdown("### Параметры подпитки")
    S_feed = st.number_input("Концентрация субстрата в подпитке, г/л (S_feed)",
                             min_value=50.0, max_value=1000.0,
                             value=200.0, step=10.0)
    F_vol = st.number_input("Удельный поток подпитки, л/(л·ч) (F/V)",
                            min_value=0.0, max_value=0.01, value=0.00075,
                            format="%.5f", step=0.0001)

    st.markdown("### Кинетика штамма CHO")
    mu = st.number_input("Максимальная удельная скорость роста, 1/ч (μmax)",
                         min_value=0.001, max_value=0.1, value=0.040,
                         format="%.3f", step=0.001)
    Ks = st.number_input("Константа насыщения по субстрату, г/л (Ks)",
                         min_value=0.01, max_value=5.0, value=0.50,
                         format="%.2f", step=0.01)
    Kd = st.number_input("Константа отмирания клеток, 1/ч (Kd)",
                         min_value=0.000, max_value=0.05, value=0.004,
                         format="%.3f", step=0.001)
    qp = st.number_input("Удельная продуктивность, г/(10⁹ кл·ч) (qp)",
                         min_value=0.0001, max_value=0.01, value=0.0012,
                         format="%.4f", step=0.0001)

    st.markdown("### Численная схема")
    method = st.radio("Метод интегрирования", ["RK4", "Euler"],
                      horizontal=True)
    dt = st.select_slider("Шаг интегрирования, ч (Δt)",
                          options=[0.01, 0.05, 0.1, 0.2], value=0.1)

res = simulate_fed_batch(
    days=int(days), mu_max=mu, Ks=Ks, Kd=Kd, q_p=qp,
    reactor_volume_l=V_l,
    scale_loss_per_decade_pct=loss_pct,
    S_feed=S_feed, F_vol=F_vol,
    dt=dt, method=method,
)

with col_out:
    st.subheader("Кинетические кривые")
    st.line_chart(
        res["history"].set_index("Сутки (t)")[
            ["Клетки CHO, ×10⁹ кл/л (X)",
             "Глюкоза, г/л (S)",
             "Леканемаб, г/л (P)"]]
    )

    st.subheader("Результаты серии")
    st.table(pd.DataFrame({
        "Параметр": ["Титр с учётом масштаба", "Титр без масштабирования",
                     "Масса целевого белка за серию", "Финальная плотность клеток",
                     "Остаточная глюкоза", "Кратность увеличения объёма"],
        "Обозначение": ["P_scaled", "P", "m", "X_final", "S_final", "V_final / V₀"],
        "Значение": [f"{res['titer']:.3f}", f"{res['titer_raw']:.3f}",
                     f"{res['mass_g']:.1f}", f"{res['X']:.2f}",
                     f"{res['S']:.2f}", f"{res['V']:.4f}"],
        "Ед. изм.": ["г/л", "г/л", "г", "×10⁹ кл/л", "г/л", "л/л₀"],
    }))

    st.subheader("Сравнение RK4 и Эйлера: траектории X(t), S(t), P(t)")
    st.caption(
        "Одинаковые параметры модели, одинаковый шаг Δt. Разница между "
        "методами проявляется в динамике роста и усиливается с ростом Δt. "
        "Для наглядности Δt по умолчанию увеличен до 0,5 ч."
    )

    dt_cmp = st.select_slider(
        "Шаг интегрирования для сравнения, ч (Δt)",
        options=[0.05, 0.1, 0.2, 0.5, 1.0], value=0.5,
    )

    df_cmp, diffs = compare_trajectories(
        int(days), mu, Ks, Kd, qp, V_l, S_feed, F_vol, dt_cmp,
    )

    st.line_chart(df_cmp.set_index("Сутки (t)"))

    st.table(pd.DataFrame({
        "Показатель расхождения (Эйлер vs RK4)": [
            "Конечный титр, % (ΔP)",
            "Финальная плотность клеток, % (ΔX)",
            "Остаточная глюкоза, г/л (ΔS)",
        ],
        "Значение": [
            f"{diffs['ΔP, %']:+.3f}",
            f"{diffs['ΔX, %']:+.3f}",
            f"{diffs['ΔS, г/л']:+.4f}",
        ],
    }))

    st.subheader("Проверка численного метода: финальные значения")
    st.caption(
        "Эталон — RK4 при Δt = 0,01 ч. Сравнение финальных значений "
        "при разных методах и шагах."
    )
    st.dataframe(
        compare_methods(int(days), mu, Ks, Kd, qp, V_l, S_feed, F_vol),
        use_container_width=True,
    )
