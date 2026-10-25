import streamlit as st
import math
import pandas as pd

st.set_page_config(page_title="Цифровой двойник: леканемаб", layout="wide")

def simulate_fed_batch(
    days, mu_max, Ks, Kd, q_p,
    reactor_volume_l,
    scale_loss_per_decade_pct=10.0,
    S_feed=400.0, F_vol=0.00075,
    S0=6.0, X0=0.4, P0=0.0,
    t_feed_start=48.0,
    dt=0.1, method="RK4",
):
    """Fed-Batch с ростом объёма. Концентрации — на 1 л начального объёма."""
    Xmax = 16.0
    Y_xs = 0.45          # 10^9 кл/г (X в 10^9 кл/л)
    m_s = 0.001          # г/(10^9 кл·ч)

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
                "День": int(round(t / 24)),
                "X, ×10⁹ кл/л": round(X, 2),
                "S, г/л": round(max(S, 0.0), 2),
                "P, г/л": round(P, 3),
                "μ, ч⁻¹": round(mu_now, 4),
                "V, л/л₀": round(V, 4),
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
        else:  # Euler
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


def compare_methods(**kw):
    rows = []
    for method, dt in [("RK4", 0.1), ("RK4", 0.05),
                       ("Euler", 0.1), ("Euler", 0.01)]:
        r = simulate_fed_batch(dt=dt, method=method, **kw)
        rows.append({"Метод": method, "dt, ч": dt,
                     "P, г/л": round(r["titer_raw"], 4),
                     "X, ×10⁹": round(r["X"], 4),
                     "V, л/л₀": round(r["V"], 4)})
    return pd.DataFrame(rows)


# --- UI ---
st.title("Цифровой двойник: Fed-Batch леканемаба")
col_in, col_out = st.columns([1, 1.3])

with col_in:
    st.markdown("### Масштаб")
    V_l = st.number_input("Номинальный объём, л", 10.0, 2000.0, 200.0, 10.0)
    days = st.number_input("Сутки", 1, 30, 10, 1)
    loss_pct = st.slider("Потеря титра, %/декада (10→13% при 200 л)",
                         0.0, 30.0, 10.0, 0.5)

    st.markdown("### Подпитка")
    S_feed = st.number_input("S_feed, г/л", value=400.0, step=10.0)
    F_vol = st.number_input("F_vol, л/(л·ч)", value=0.00075,
                             format="%.5f", step=0.0001)

    st.markdown("### Кинетика CHO")
    mu = st.number_input("μmax, 1/ч", value=0.040, format="%.3f")
    Ks = st.number_input("Ks, г/л", value=0.50, format="%.2f")
    Kd = st.number_input("Kd, 1/ч", value=0.004, format="%.3f")
    qp = st.number_input("qp, г/(10⁹кл·ч)", value=0.0012, format="%.4f")

    method = st.radio("Метод", ["RK4", "Euler"], horizontal=True)
    dt = st.select_slider("dt, ч", [0.01, 0.05, 0.1, 0.2], value=0.1)

res = simulate_fed_batch(days=int(days), mu_max=mu, Ks=Ks, Kd=Kd, q_p=qp,
                          reactor_volume_l=V_l,
                          scale_loss_per_decade_pct=loss_pct,
                          S_feed=S_feed, F_vol=F_vol,
                          dt=dt, method=method)

with col_out:
    st.subheader("Кинетические кривые")
    st.line_chart(res["history"].set_index("День")
                  [["X, ×10⁹ кл/л", "S, г/л", "P, г/л"]])

    st.table(pd.DataFrame({
        "Параметр": ["Титр (с масштабом)", "Титр (без масштаба)",
                     "Масса за серию, г", "X конечная",
                     "S конечная", "V_final/V₀"],
        "Значение": [f"{res['titer']:.3f} г/л", f"{res['titer_raw']:.3f} г/л",
                     f"{res['mass_g']:.1f} г", f"{res['X']:.2f}",
                     f"{res['S']:.2f}", f"{res['V']:.4f}"]}))

    st.subheader("Проверка численного метода")
    st.dataframe(compare_methods(days=int(days), mu_max=mu, Ks=Ks, Kd=Kd,
                                  q_p=qp, reactor_volume_l=V_l,
                                  S_feed=S_feed, F_vol=F_vol))
