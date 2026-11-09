import streamlit as st
import math
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

st.set_page_config(page_title="Расчётная модель: леканемаб", layout="wide")


def simulate_fed_batch(
    days, mu_max, Ks, Kd, q_p,
    V0_l,
    V_max_l=200.0,
    scale_loss_per_decade_pct=10.0,
    S_feed=300.0, feed_rate=0.172, S_trigger=4.0,
    S0=6.0, X0=0.4, P0=0.0, Xd0=0.0,
    dt=0.1,
):
    """Fed-Batch с ростом объёма и подпиткой по остаточной глюкозе.

    V0_l      — начальный объём культуральной жидкости (после инокуляции), л.
    V_max_l   — рабочий объём аппарата, л; по его достижении подпитка прекращается.
    S_trigger — порог включения подпитки по остаточной глюкозе, г/л.
    Продукция: Luedeking-Piret при alpha = 0 (beta = qp).
    Гибель: Kd первого порядка; погибшие клетки учитываются отдельным
    балансом Xd (лизис не моделируется).
    """
    Xmax = 16.0
    Y_xs = 0.45
    m_s = 0.001

    F_vol = feed_rate / S_feed if S_feed > 0 else 0.0
    V_ratio_max = V_max_l / V0_l if V0_l > 0 else 1.0

    V, X, Xd, S, P = 1.0, X0, Xd0, S0, P0
    total_hours = days * 24
    history = []
    # плотный трек для графиков: запись каждые 0,5 ч
    trace = {"t": [], "X": [], "Xd": [], "S": [], "P": [],
             "V_l": [], "mu": [], "viab": []}
    t_feed_start = None      # момент первого включения подпитки, ч
    t_feed_cutoff = None     # момент прекращения подпитки по объёму, ч

    def feeding_on(S_, V_):
        return (S_ < S_trigger) and (V_ < V_ratio_max)

    def deriv(X, Xd, S, P, V, t):
        F = F_vol if feeding_on(S, V) else 0.0
        D = F / V
        if S <= 0.0:
            return (-Kd - D) * X, Kd * X - D * Xd, D * S_feed, -D * P, F
        mu = mu_max * (S / (Ks + S)) * max(0.0, 1.0 - X / Xmax)
        dX = (mu - Kd - D) * X
        dXd = Kd * X - D * Xd
        dS = -(mu * X) / Y_xs - m_s * X + D * (S_feed - S)
        dP = q_p * X - D * P
        return dX, dXd, dS, dP, F

    n_steps = int(total_hours / dt) + 1
    for step in range(n_steps):
        t = step * dt
        on = feeding_on(S, V)
        if t_feed_start is None and on:
            t_feed_start = t
        if (t_feed_cutoff is None and t_feed_start is not None
                and not on and V >= V_ratio_max - 1e-9):
            t_feed_cutoff = t

        mu_now = (mu_max * S / (Ks + S) * max(0.0, 1.0 - X / Xmax)) if S > 0 else 0.0
        viab_now = 100.0 * X / (X + Xd) if (X + Xd) > 0 else 0.0

        if step % max(1, int(round(0.5 / dt))) == 0:
            trace["t"].append(t / 24.0)
            trace["X"].append(X)
            trace["Xd"].append(Xd)
            trace["S"].append(max(S, 0.0))
            trace["P"].append(P)
            trace["V_l"].append(V * V0_l)
            trace["mu"].append(mu_now)
            trace["viab"].append(viab_now)

        if step % int(round(24 / dt)) == 0:
            history.append({
                "Сутки (t)": int(round(t / 24)),
                "Живые клетки, ×10⁹ кл/л (X)": round(X, 3),
                "Погибшие клетки, ×10⁹ кл/л (Xd)": round(Xd, 3),
                "Жизнеспособность, % (V_viab)": round(viab_now, 1),
                "Глюкоза, г/л (S)": round(max(S, 0.0), 3),
                "Леканемаб, г/л (P)": round(P, 4),
                "Скорость роста, ч⁻¹ (μ)": round(mu_now, 5),
                "Кратность объёма (V/V₀)": round(V, 4),
                "Подпитка": "вкл." if on else "выкл.",
            })

        k1 = deriv(X, Xd, S, P, V, t)
        k2 = deriv(X + 0.5*dt*k1[0], Xd + 0.5*dt*k1[1], S + 0.5*dt*k1[2],
                   P + 0.5*dt*k1[3], V + 0.5*dt*k1[4], t + 0.5*dt)
        k3 = deriv(X + 0.5*dt*k2[0], Xd + 0.5*dt*k2[1], S + 0.5*dt*k2[2],
                   P + 0.5*dt*k2[3], V + 0.5*dt*k2[4], t + 0.5*dt)
        k4 = deriv(X + dt*k3[0], Xd + dt*k3[1], S + dt*k3[2],
                   P + dt*k3[3], V + dt*k3[4], t + dt)
        for i, (a, b, c, d) in enumerate(zip(k1, k2, k3, k4)):
            inc = (dt/6.0) * (a + 2*b + 2*c + d)
            if i == 0: X += inc
            elif i == 1: Xd += inc
            elif i == 2: S += inc
            elif i == 3: P += inc
            else: V += inc

        X, Xd, S, P = max(X, 0.0), max(Xd, 0.0), max(S, 0.0), max(P, 0.0)
        V = min(V, V_ratio_max)

    V_final_l = V * V0_l
    viability = 100.0 * X / (X + Xd) if (X + Xd) > 0 else 0.0

    # Поправка на масштабирование — принятое проектное допущение
    if V_final_l > 10.0:
        decades = math.log10(V_final_l / 10.0)
        loss = (scale_loss_per_decade_pct / 100.0) * decades
        P_s = P * max(0.0, 1.0 - loss)
    else:
        P_s = P

    return {
        "titer": P_s, "titer_raw": P, "mass_g": P_s * V_final_l,
        "V0_l": V0_l, "V_max_l": V_max_l, "V_final_l": V_final_l, "V_ratio": V,
        "X": X, "Xd": Xd, "S": S, "V": V, "viability": viability,
        "t_feed_start": t_feed_start, "t_feed_cutoff": t_feed_cutoff,
        "loss_pct": (1.0 - P_s / P) * 100.0 if P > 0 else 0.0,
        "F_vol": F_vol, "history": pd.DataFrame(history),
        "trace": trace, "days": days,
    }


# ----------------------------------------------------------------------------
# Оформление графиков
# ----------------------------------------------------------------------------
BLUE, VERM, TEAL, PURP = "#1f5fa8", "#c4461e", "#0f8a74", "#7a4ea3"
INK, MUTED, GRID = "#1a1a1a", "#5a5a5a", "#d8d8d4"

plt.rcParams.update({
    "font.size": 10,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


class _CommaFormatter(matplotlib.ticker.ScalarFormatter):
    """Десятичная запятая на осях, число знаков — как у штатного форматтера."""

    def __call__(self, x, pos=None):
        return super().__call__(x, pos).replace(".", ",")


def _style(ax, days, t_feed_days=None, label_feed=False, feed_y=0.94):
    """Единое оформление: сетка только по Y, две рамки, метка подпитки."""
    ax.grid(True, axis="y", color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlim(0, days)
    step = 1 if days <= 14 else (2 if days <= 28 else 5)
    ax.set_xticks(range(0, days + 1, step))
    ax.yaxis.set_major_formatter(_CommaFormatter())
    if t_feed_days is not None:
        ax.axvline(t_feed_days, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
        if label_feed:
            ax.annotate("включение подпитки",
                        xy=(t_feed_days, feed_y), xycoords=("data", "axes fraction"),
                        xytext=(6, 0), textcoords="offset points",
                        fontsize=8.6, color=MUTED, va="center")


def _marker_every(trace):
    """Маркеры раз в сутки независимо от шага записи."""
    t = trace["t"]
    if len(t) < 2:
        return 1
    return max(1, int(round(1.0 / (t[1] - t[0]))))


def fig_kinetics(res):
    """Рисунок: (а) клетки, (б) глюкоза и продукт."""
    tr, days = res["trace"], res["days"]
    tf = res["t_feed_start"] / 24.0 if res["t_feed_start"] is not None else None
    me = _marker_every(tr)

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 7.4), sharex=True)

    a1.plot(tr["t"], tr["X"], color=BLUE, lw=2, marker="o", ms=5,
            markevery=me, zorder=3, label="жизнеспособные клетки $X$")
    a1.plot(tr["t"], tr["Xd"], color=VERM, lw=2, ls=(0, (6, 3)), marker="s",
            ms=5, markevery=me, zorder=3, label="погибшие клетки $X_d$")
    a1.set_ylabel("Концентрация клеток,\n$\\times10^{9}$ кл/л")
    a1.set_ylim(bottom=0)
    a1.legend(frameon=False, loc="center right", fontsize=9.2, handlelength=2.6)
    a1.set_title("а", loc="left", fontsize=10, color=INK, pad=6)
    _style(a1, days, tf, label_feed=True)

    a2.plot(tr["t"], tr["S"], color=TEAL, lw=2, marker="^", ms=5.5,
            markevery=me, zorder=3, label="глюкоза $S$")
    a2.plot(tr["t"], tr["P"], color=PURP, lw=2, ls=(0, (1.5, 2)), marker="D",
            ms=4.5, markevery=me, zorder=3, label="леканемаб $P$")
    a2.set_ylabel("Концентрация, г/л")
    a2.set_xlabel("Время, сут")
    a2.set_ylim(bottom=0)
    a2.legend(frameon=False, loc="upper right", fontsize=9.2, handlelength=2.6)
    a2.set_title("б", loc="left", fontsize=10, color=INK, pad=6)
    _style(a2, days, tf)

    fig.tight_layout(h_pad=1.6)
    return fig


def fig_volume(res):
    """Рисунок: объём культуральной жидкости и рабочий объём аппарата."""
    tr, days = res["trace"], res["days"]
    tf = res["t_feed_start"] / 24.0 if res["t_feed_start"] is not None else None
    me = _marker_every(tr)
    V_max_l, V_fin = res["V_max_l"], tr["V_l"][-1]

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.plot(tr["t"], tr["V_l"], color=BLUE, lw=2, marker="o", ms=5,
            markevery=me, zorder=3)
    ax.axhline(V_max_l, color=VERM, lw=1.4, ls=(0, (5, 3)), zorder=2)
    ax.annotate(f"рабочий объём аппарата {V_max_l:.0f} л",
                xy=(0.02, V_max_l), xycoords=("axes fraction", "data"),
                xytext=(0, -5), textcoords="offset points",
                fontsize=9, color=VERM, va="top")
    ax.annotate(("%.1f л" % V_fin).replace(".", ","),
                xy=(days, V_fin), xytext=(-10, -16), textcoords="offset points",
                fontsize=9, color=INK, ha="right", va="top")
    ax.set_ylabel("Объём культуральной\nжидкости, л")
    ax.set_xlabel("Время, сут")
    span = max(V_max_l - res["V0_l"], 1.0)
    ax.set_ylim(res["V0_l"] - 0.12 * span, V_max_l + 0.12 * span)
    _style(ax, days, tf, label_feed=True, feed_y=0.70)

    fig.tight_layout()
    return fig


def fig_state(res):
    """Рисунок: жизнеспособность и удельная скорость роста (раздельные оси)."""
    tr, days = res["trace"], res["days"]
    tf = res["t_feed_start"] / 24.0 if res["t_feed_start"] is not None else None
    me = _marker_every(tr)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.1))

    a1.plot(tr["t"], tr["viab"], color=TEAL, lw=2, marker="o", ms=4.5,
            markevery=me, zorder=3)
    a1.axhline(95.0, color=VERM, lw=1.2, ls=(0, (5, 3)), zorder=2)
    a1.annotate("критерий 95 %", xy=(0.97, 95.0),
                xycoords=("axes fraction", "data"),
                xytext=(0, 4), textcoords="offset points",
                fontsize=8.6, color=VERM, va="bottom", ha="right")
    a1.set_ylabel("Жизнеспособность, %")
    a1.set_xlabel("Время, сут")
    a1.set_ylim(0, 105)
    _style(a1, days, tf)

    a2.plot(tr["t"], tr["mu"], color=PURP, lw=2, marker="D", ms=4,
            markevery=me, zorder=3)
    a2.set_ylabel("Удельная скорость\nроста $\\mu$, ч$^{-1}$")
    a2.set_xlabel("Время, сут")
    a2.set_ylim(bottom=0)
    _style(a2, days, tf)

    fig.tight_layout(w_pad=2.0)
    return fig


def material_balance_chain(mass_raw_g, vial_mg=500.0):
    stages = [
        ("Осветление и глубинная фильтрация", 0.95),
        ("Аффинная хроматография (Protein A)", 0.85),
        ("Тангенциальная ультрафильтрация и полировка", 0.95),
        ("Формуляция, розлив и лиофилизация", 0.98),
    ]
    rows = [{"Стадия": "Вход: неочищенный белок",
             "Выход стадии, %": 100.0, "Масса, г": round(mass_raw_g, 2)}]
    current = mass_raw_g
    for name, y in stages:
        current *= y
        rows.append({"Стадия": name, "Выход стадии, %": round(y * 100, 1),
                     "Масса, г": round(current, 2)})
    n_vials = int(current / (vial_mg / 1000.0))
    return pd.DataFrame(rows), current, n_vials


# ============================================================
st.title("Расчётная модель опытно-промышленной линии")
st.subheader("Моделирование кинетики Fed-Batch с учётом масштабирования")

col_in, col_out = st.columns([1, 1.3])

with col_in:
    st.markdown("### Параметры масштаба")
    V0_l = st.number_input(
        "Начальный объём культуральной жидкости, л (V₀)",
        min_value=10.0, max_value=2000.0, value=172.0, step=1.0,
        help="Объём после инокуляции, до начала подпитки."
    )
    V_max_l = st.number_input(
        "Рабочий объём аппарата, л (V_max)",
        min_value=10.0, max_value=2500.0, value=200.0, step=1.0,
        help="Предел заполнения аппарата: по его достижении подпитка "
             "прекращается. Для аппарата 250 л рабочий объём 200 л "
             "соответствует заполнению 80% (headspace 20% на пену и аэрацию)."
    )
    if V_max_l < V0_l:
        st.warning("Рабочий объём меньше начального: подпитка не включится.")

    days = st.number_input("Длительность культивирования, сутки (t)",
                           min_value=1, max_value=14, value=10, step=1)
    loss_pct = st.slider(
        "Потеря титра при масштабировании, %/декада (k)",
        min_value=0.0, max_value=30.0, value=10.0, step=0.5,
        help="Принятое проектное допущение: консервативная поправка на "
             "ухудшение перемешивания, переноса кислорода и удаления CO₂ "
             "при увеличении рабочего объёма. Экспериментально не "
             "подтверждена."
    )

    st.markdown("### Параметры подпитки")
    S_trigger = st.number_input(
        "Порог включения подпитки по глюкозе, г/л (S_trigger)",
        min_value=0.5, max_value=10.0, value=4.0, step=0.1,
        help="Подпитка включается, когда остаточная глюкоза опускается ниже "
             "этого значения, и ведётся с постоянной удельной скоростью."
    )
    S_feed = st.number_input("Концентрация субстрата в подпитке, г/л (S_feed)",
                             min_value=50.0, max_value=1000.0,
                             value=300.0, step=10.0,
                             help="Концентрированная подпитка: чем выше S_feed, "
                                  "тем меньше прирост объёма на ту же массу глюкозы.")
    feed_rate = st.number_input(
        "Скорость подпитки по глюкозе, г/(л·ч) (feed_rate)",
        min_value=0.05, max_value=0.30, value=0.172, format="%.3f", step=0.001,
        help="Принятое проектное допущение; экспериментально не обосновано."
    )
    st.caption(f"F/V = {feed_rate / S_feed if S_feed > 0 else 0:.5f} л/(л·ч).")

    st.markdown("### Кинетика штамма CHO")
    mu = st.number_input("Максимальная удельная скорость роста, 1/ч (μmax)",
                         min_value=0.001, max_value=0.1, value=0.040,
                         format="%.3f", step=0.001)
    Ks = st.number_input("Константа насыщения по субстрату, г/л (Ks)",
                         min_value=0.01, max_value=5.0, value=0.50,
                         format="%.2f", step=0.01)
    Kd = st.number_input(
        "Константа гибели клеток, 1/ч (Kd)",
        min_value=0.000, max_value=0.05, value=0.004, format="%.3f", step=0.001,
        help="dX/dt = (μ – Kd – D)·X; dXd/dt = Kd·X – D·Xd. Принята "
             "постоянной на весь процесс; лизис не моделируется."
    )
    qp = st.number_input(
        "Удельная продуктивность, г/(10⁹ кл·ч) (qp)",
        min_value=0.0001, max_value=0.01, value=0.0012, format="%.4f", step=0.0001,
        help="Luedeking-Piret при α = 0: dP/dt = qp·X – D·P."
    )

res = simulate_fed_batch(
    days=int(days), mu_max=mu, Ks=Ks, Kd=Kd, q_p=qp,
    V0_l=V0_l, V_max_l=V_max_l, scale_loss_per_decade_pct=loss_pct,
    S_feed=S_feed, feed_rate=feed_rate, S_trigger=S_trigger,
)

with col_out:
    st.subheader("Кинетические кривые")
    tfs = f"{res['t_feed_start']:.0f} ч" if res['t_feed_start'] is not None else "не включалась"
    st.caption(
        f"Подпитка включена с {tfs};  V₀ = {res['V0_l']:.0f} л → "
        f"V_final = {res['V_final_l']:.1f} л (предел {res['V_max_l']:.0f} л).  "
        f"Модель: Моно + Fed-Batch + Kd + Luedeking-Piret (α = 0)."
    )
    if res['t_feed_cutoff'] is not None:
        st.warning(
            f"Подпитка прекращена на {res['t_feed_cutoff']:.0f} ч по достижении "
            f"рабочего объёма, до конца культивирования. В оставшееся время "
            f"культура работает на остаточной глюкозе. Увеличьте S_feed или "
            f"V_max, либо сократите длительность."
        )

    for builder, cap in (
        (fig_kinetics, "Концентрация клеток (а) и субстрата с продуктом (б). "
                       "Маркеры нанесены раз в сутки, кривая рассчитана с "
                       "шагом 0,1 ч и записана каждые 0,5 ч."),
        (fig_volume, "Рост объёма культуральной жидкости за счёт подпитки."),
        (fig_state, "Жизнеспособность и удельная скорость роста."),
    ):
        f = builder(res)
        st.pyplot(f, use_container_width=True)
        st.caption(cap)
        plt.close(f)

    st.subheader("Результаты серии")
    st.table(pd.DataFrame({
        "Параметр": ["Титр без масштаба", "Потеря при масштабировании",
                     "Титр с учётом масштаба",
                     "Начальный объём культуральной жидкости",
                     "Рабочий объём аппарата",
                     "Финальный объём культуральной жидкости",
                     "Масса неочищенного белка за серию",
                     "Начало подпитки",
                     "X финальная", "Xd финальная",
                     "Жизнеспособность на выгрузке", "S финальная"],
        "Обозначение": ["P", "ΔP", "P_scaled", "V₀", "V_max", "V_final", "m",
                        "t_feed", "X_final", "Xd_final", "V_viab", "S_final"],
        "Значение": [f"{res['titer_raw']:.3f}", f"{res['loss_pct']:.1f}",
                     f"{res['titer']:.3f}", f"{res['V0_l']:.0f}",
                     f"{res['V_max_l']:.0f}", f"{res['V_final_l']:.1f}",
                     f"{res['mass_g']:.1f}",
                     f"{res['t_feed_start']:.0f}" if res['t_feed_start'] is not None else "—",
                     f"{res['X']:.2f}", f"{res['Xd']:.2f}",
                     f"{res['viability']:.1f}", f"{res['S']:.2f}"],
        "Ед. изм.": ["г/л", "%", "г/л", "л", "л", "л", "г", "ч",
                     "×10⁹ кл/л", "×10⁹ кл/л", "%", "г/л"],
    }))

    st.subheader("История расчёта по суткам")
    st.dataframe(res["history"], hide_index=True)

st.subheader("Материальный баланс серии")
st.caption("Цепочка выходов по стадиям Downstream. Вход — масса неочищенного "
           "белка из результатов серии.")
balance_df, final_mass, n_vials = material_balance_chain(res["mass_g"])
st.table(balance_df)

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
col3.metric("V_final, л", f"{res['V_final_l']:.1f}")
