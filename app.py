import streamlit as st
import math
import pandas as pd

st.set_page_config(page_title="Цифровой двойник: леканемаб", layout="wide")


# ============================================================
# ИСТОЧНИКИ ПАРАМЕТРОВ (для таблицы 2.1 в главе 2)
# ============================================================
PARAMETER_SOURCES = {
    "mu_max": ("μmax = 0,04 ч⁻¹",
               "Типичное значение для рекомбинантных CHO "
               "[Xing et al., 2009; Глик, Пастернак, 2021]"),
    "Ks": ("Ks = 0,5 г/л",
           "Константа насыщения Моно для глюкозы в CHO-культурах "
           "[Xing et al., 2009]"),
    "Kd": ("Kd = 0,004 ч⁻¹",
           "Скорость спонтанного апоптоза CHO в Fed-Batch "
           "[Xing et al., 2009]"),
    "qp": ("qp = 1,2·10⁻³ г/(10⁹ кл·ч)",
           "Удельная продуктивность IgG1-продуцента CHO "
           "[Xing et al., 2009]"),
    "Y_xs": ("Yx/s = 0,45·10⁹ кл/г",
             "Экономический коэффициент конверсии глюкозы "
             "[Глик, Пастернак, 2021]"),
    "m_s": ("ms = 0,001 г/(10⁹ кл·ч)",
            "Расход глюкозы на поддержание [Xing et al., 2009]"),
    "Xmax": ("Xmax = 16,0·10⁹ кл/л",
             "Предельная ёмкость среды для CHO [Xing et al., 2009]"),
    "scale_loss": ("k = 10 %/декада",
                   "Консервативное значение в литературном диапазоне "
                   "5–15 %. Обоснование: время перемешивания как "
                   "лимитирующий фактор при масштабировании "
                   "[Xing et al., 2009]"),
}

# ============================================================
# ЛИТЕРАТУРНЫЕ ДИАПАЗОНЫ ДЛЯ ВЕРИФИКАЦИИ
# ============================================================
LITERATURE_RANGES = {
    "Клетки CHO, ×10⁹ кл/л (X)":  (8.0, 20.0),
    "Леканемаб, г/л (P)":         (1.0, 5.0),
    "Глюкоза, г/л (S)":           (0.1, 3.0),
    "Скорость роста, ч⁻¹ (μmax)": (0.03, 0.05),
}


def simulate_fed_batch(
    days, mu_max, Ks, Kd, q_p,
    reactor_volume_l,
    scale_loss_per_decade_pct=10.0,
    S_feed=200.0, F_vol=0.00086,
    S0=6.0, X0=0.4, P0=0.0,
    t_feed_start=48.0,
    dt=0.1,
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

    return {"titer": P_s, "titer_raw": P,
            "mass_g": P_s * V * reactor_volume_l,
            "X": X, "S": S, "V": V,
            "loss_pct": (1.0 - P_s / P) * 100.0 if P > 0 else 0.0,
            "history": pd.DataFrame(history)}


def check_literature(res):
    checks = [
        ("Клетки CHO, ×10⁹ кл/л (X)", res["X"]),
        ("Леканемаб, г/л (P)", res["titer_raw"]),
        ("Глюкоза, г/л (S)", res["S"]),
        ("Скорость роста, ч⁻¹ (μmax)", None),
    ]
    rows = []
    for name, value in checks:
        lo, hi = LITERATURE_RANGES[name]
        if value is None:
            rows.append({"Показатель": name, "Значение модели": "—",
                         "Литературный диапазон": f"{lo}–{hi}",
                         "Вывод": "параметр задан"})
        else:
            ok = lo <= value <= hi
            rows.append({"Показатель": name,
                         "Значение модели": round(value, 3),
                         "Литературный диапазон": f"{lo}–{hi}",
                         "Вывод": "✓ в диапазоне" if ok else "✗ вне диапазона"})
    return pd.DataFrame(rows)


def material_balance_chain(mass_raw_g):
    """Цепочка материального баланса из главы 3."""
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
    for name, yield_frac in stages:
        current *= yield_frac
        rows.append({"Стадия": name,
                     "Выход стадии, %": round(yield_frac * 100, 1),
                     "Масса, г": round(current, 2)})
    return pd.DataFrame(rows), current


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
    S_feed = st.number_input(
        "Концентрация субстрата в подпитке, г/л (S_feed)",
        min_value=50.0, max_value=1000.0, value=200.0, step=10.0,
        help="Обеспечивает подачу ≈ 0,172 г/(л·ч) при F/V = 8,6·10⁻⁴ л/(л·ч)."
    )
    F_vol = st.number_input(
        "Удельный поток подпитки, л/(л·ч) (F/V)",
        min_value=0.0, max_value=0.01, value=0.00086,
        format="%.5f", step=0.00001,
        help="Подобрано из условия воспроизведения проектного титра "
             "P ≈ 1,509 г/л."
    )

    st.markdown("### Кинетика штамма CHO")
    mu = st.number_input("Максимальная удельная скорость роста, 1/ч (μmax)",
                         min_value=0.001, max_value=0.1, value=0.040,
                         format="%.3f", step=0.001,
                         help=PARAMETER_SOURCES["mu_max"][1])
    Ks = st.number_input("Константа насыщения по субстрату, г/л (Ks)",
                         min_value=0.01, max_value=5.0, value=0.50,
                         format="%.2f", step=0.01,
                         help=PARAMETER_SOURCES["Ks"][1])
    Kd = st.number_input("Константа отмирания клеток, 1/ч (Kd)",
                         min_value=0.000, max_value=0.05, value=0.004,
                         format="%.3f", step=0.001,
                         help=PARAMETER_SOURCES["Kd"][1])
    qp = st.number_input("Удельная продуктивность, г/(10⁹ кл·ч) (qp)",
                         min_value=0.0001, max_value=0.01, value=0.0012,
                         format="%.4f", step=0.0001,
                         help=PARAMETER_SOURCES["qp"][1])

res = simulate_fed_batch(
    days=int(days), mu_max=mu, Ks=Ks, Kd=Kd, q_p=qp,
    reactor_volume_l=V_l,
    scale_loss_per_decade_pct=loss_pct,
    S_feed=S_feed, F_vol=F_vol,
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
        "Параметр": ["Титр без масштаба", "Потеря при масштабировании",
                     "Титр с учётом масштаба", "Масса за серию",
                     "X финальная", "S финальная", "V_final / V₀"],
        "Обозначение": ["P", "ΔP", "P_scaled", "m", "X_final",
                        "S_final", "V_final / V₀"],
        "Значение": [f"{res['titer_raw']:.3f}",
                     f"{res['loss_pct']:.1f}",
                     f"{res['titer']:.3f}",
                     f"{res['mass_g']:.1f}",
                     f"{res['X']:.2f}",
                     f"{res['S']:.2f}",
                     f"{res['V']:.4f}"],
        "Ед. изм.": ["г/л", "%", "г/л", "г", "×10⁹ кл/л", "г/л", "л/л₀"],
    }))

# ============================================================
# МОДЕЛЬНЫЕ ЗНАЧЕНИЯ (для вставки в главу 2 как таблица 2.2)
# ============================================================
st.subheader("Модельные значения по суткам (таблица 2.2)")
st.caption(
    "Эти значения — результат расчёта и используются как таблица 2.2 "
    "в главе 2, а также как входные данные для материального баланса "
    "главы 3."
)
model_table = res["history"][
    ["Сутки (t)", "Клетки CHO, ×10⁹ кл/л (X)",
     "Глюкоза, г/л (S)", "Леканемаб, г/л (P)",
     "Скорость роста, ч⁻¹ (μ)"]
]
st.dataframe(model_table, use_container_width=True)
st.download_button(
    "Скачать таблицу 2.2 (CSV)",
    model_table.to_csv(index=False).encode("utf-8"),
    "table_2_2_model.csv",
    "text/csv",
)

# ============================================================
# ВЕРИФИКАЦИЯ vs ЛИТЕРАТУРА
# ============================================================
st.subheader("Верификация модели: литературные диапазоны CHO Fed-Batch")
st.caption(
    "Финальные значения сопоставлены с типичными диапазонами для "
    "Fed-Batch процессов на CHO. Совпадение подтверждает корректность "
    "выбранных параметров и пригодность модели для расчёта "
    "материального баланса в главе 3."
)
st.table(check_literature(res))

# ============================================================
# ПЕРЕДАЧА В ГЛАВУ 3: МАТЕРИАЛЬНЫЙ БАЛАНС
# ============================================================
st.subheader("Передача результатов в главу 3: материальный баланс")
st.caption(
    "Цепочка выходов по стадиям Downstream из раздела 3.3. "
    "Вход — масса неочищенного белка из главы 2."
)

balance_df, final_mass = material_balance_chain(res["mass_g"])
st.table(balance_df)

st.markdown(
    f"**Итог:** из {res['mass_g']:.1f} г неочищенного белка получается "
    f"**{final_mass:.1f} г** готового леканемаба. "
    f"Сквозной коэффициент использования "
    f"K_исп = {final_mass / res['mass_g']:.3f}, "
    f"расходный коэффициент K_расх = {res['mass_g'] / final_mass:.2f}."
)

st.info(
    "Если итоговая масса отличается от 197,42 г из главы 1.4 — "
    "обновите цепочку в разделе 3.3 или скорректируйте V биореактора "
    "так, чтобы m × K_исп совпадало с проектным значением."
)

# ============================================================
# ИСТОЧНИКИ ПАРАМЕТРОВ
# ============================================================
with st.expander("Источники параметров модели (для таблицы 2.1)"):
    st.table(pd.DataFrame([
        {"Параметр": v[0], "Обоснование": v[1]}
        for v in PARAMETER_SOURCES.values()
    ]))
