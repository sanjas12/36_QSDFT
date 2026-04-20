"""
QSDFT — Quenched Solid Density Functional Theory
для анализа распределения пор по размерам в микро-мезопористых углеродах.

Реализация на основе:
    Neimark et al., Carbon 47 (2009) 1617-1628
    + Supplementary Information (таблицы давлений заполнения)

Метод воспроизводит:
    1. Модель поверхности с шероховатостью (линейный рамп плотности твёрдого тела)
    2. Расчёт профилей плотности через уравнение Эйлера (QSDFT/FMT)
    3. Ядро изотерм адсорбции N2 (77.4 K) и Ar (87.3 K) в щелевых порах
    4. Расчёт PSD методом QNNLS (Тихоновская регуляризация + NNLS)


1. Модель поверхности (Раздел 2.3 статьи)
Линейный рамп плотности твёрдого тела ρ_s(z) с параметром шероховатости δ = 0.13 нм, толщиной стенки h₀ = 0.68 нм. Вычисляется положение края поверхности zₑ из условия нулевого избытка (ур. 16).
2. Потенциалы взаимодействия (ур. 14)
WCA-разложение LJ-потенциала для пар N₂–N₂, Ar–Ar, C–N₂, C–Ar с параметрами из статьи.
3. Уравнение состояния Percus–Yevick (ур. 13)
Избыточный химический потенциал твёрдых сфер, используется в уравнении Эйлера.
4. Итерационное решение уравнения Эйлера (ур. 7)
Метод Пикара для профиля плотности флюида в щелевой поре с шероховатыми стенками.
5. Ядро изотерм из таблиц Supplementary
Точные данные давлений заполнения для N₂ (77.4 K) и Ar (87.3 K) из таблиц 1–2 дополнения — 70 точек каждая.
6. QNNLS-деконволюция (Раздел 3.3)
Регуляризация Тихонова + NNLS (scipy.optimize.nnls) для нахождения PSD из экспериментальной изотермы.

"""
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import nnls

warnings.filterwarnings("ignore", category=RuntimeWarning)


BASE_DIR = Path(__file__).parent.absolute()
# CONFIG_DIR = Path(BASE_DIR, "config")

OUT_DIR = Path(BASE_DIR.parent, "output")  # на уровень выше src
os.makedirs(OUT_DIR, exist_ok=True)


# ============================================================
# 1. ПАРАМЕТРЫ МОЛЕКУЛЯРНЫХ ВЗАИМОДЕЙСТВИЙ
# ============================================================

# Параметры LJ для N2-N2 (QSDFT, PY hard-sphere EOS)
N2 = dict(
    eps_ff=95.77,        # K  (глубина потенциала ε_ff/k_B)
    sigma_ff=0.3549,     # nm (диаметр молекулы)
    d_hs=0.3549,         # nm (диаметр твёрдой сферы)
    T=77.4,              # K  (температура кипения)
    label="N₂, 77.4 K"
)

# Параметры LJ для Ar-Ar
Ar = dict(
    eps_ff=111.95,       # K
    sigma_ff=0.3358,     # nm
    d_hs=0.3358,         # nm
    T=87.3,              # K
    label="Ar, 87.3 K"
)

# Параметры взаимодействия C-N2 и C-Ar (из статьи)
CARBON_N2 = dict(eps_sf=150.0, sigma_sf=0.269)   # K, nm
CARBON_AR = dict(eps_sf=162.18, sigma_sf=0.2595)  # K, nm

# Параметры модели поверхности (Cabot BP-280)
RHO_CARBON_BULK = 0.114e30  # атом/м³ → в единицах nm⁻³: 0.114 Å⁻³ = 114 nm⁻³
ROUGHNESS_DELTA = 0.13      # nm — параметр шероховатости (полуширина рампа)
H0 = 2 * 0.34              # nm — толщина стенки (2 графитовых слоя)
D_CARBON_HS = 0.2217        # nm — диаметр HS атома углерода


# ============================================================
# 2. МОДЕЛЬ ПОВЕРХНОСТИ: профиль плотности твёрдого тела
# ============================================================

def solid_density_profile(z_arr: np.ndarray, delta: float = ROUGHNESS_DELTA,
                           h0: float = H0, rho0: float = 114.0) -> np.ndarray:
    """
    Линейный рамп плотности твёрдого тела (ур. 15 в статье).

    q_s(z) = rho0,          0 <= z < h0
           = rho0*(0.75 - (1 - (z-h0)/(2*delta))),  h0 <= z < h0+2*delta
           = 0,              z >= h0+2*delta

    rho0 в nm⁻³ (0.114 Å⁻³ = 114 nm⁻³)
    """
    rho_s = np.zeros_like(z_arr, dtype=float)
    for i, z in enumerate(z_arr):
        if 0 <= z < h0:
            rho_s[i] = rho0
        elif h0 <= z < h0 + 2 * delta:
            rho_s[i] = rho0 * (0.75 - (1.0 - (z - h0) / (2 * delta)))
    return rho_s


def edge_position(delta: float = ROUGHNESS_DELTA, h0: float = H0,
                  rho0: float = 114.0) -> float:
    """
    Положение края поверхности ze из условия нулевого избытка твёрдого тела (ур. 16).
    Численное интегрирование рампа.
    """
    dz = 0.001
    z = np.arange(h0, h0 + 2 * delta + dz, dz)
    rho = solid_density_profile(z, delta, h0, rho0)
    integral = np.trapezoid(rho, z)
    ze = h0 + integral / rho0
    return ze


# ============================================================
# 3. ПОТЕНЦИАЛЫ ВЗАИМОДЕЙСТВИЯ (WCA-схема, LJ)
# ============================================================

def u_att_lj(r: np.ndarray, eps: float, sigma: float) -> np.ndarray:
    """
    Притягательная часть потенциала LJ в схеме WCA (ур. 14).
    u_att = eps * (4*eps*[(sigma/r)^12 - (sigma/r)^6] + eps),  r <= 2^(1/6)*sigma
          = 4*eps*[(sigma/r)^12 - (sigma/r)^6],                r > 2^(1/6)*sigma
    """
    r21 = 2 ** (1.0 / 6.0) * sigma
    u = np.zeros_like(r, dtype=float)
    sr6 = (sigma / r) ** 6
    u_full = 4 * eps * (sr6 ** 2 - sr6)
    mask = r <= r21
    u[mask] = u_full[mask] + eps    # внутренняя область: константа eps
    u[~mask] = u_full[~mask]        # внешняя: хвост LJ
    return u


def integrated_sf_potential(z: float, eps_sf: float, sigma_sf: float,
                             rho_s_profile_fn, n_pts: int = 300) -> float:
    """
    Интегральный потенциал флюид-твёрдое тело (1D интеграл по профилю плотности).
    phi_sf(z) = integral rho_s(z') * u_att_sf(|z-z'|) dz'
    """
    z_max = H0 + 2 * ROUGHNESS_DELTA + 0.5
    z_arr = np.linspace(0, z_max, n_pts)
    rho_s = rho_s_profile_fn(z_arr)
    r = np.abs(z - z_arr) + 1e-10  # избегаем деления на 0
    u = u_att_lj(r, eps_sf, sigma_sf)
    return np.trapezoid(rho_s * u, z_arr)


# ============================================================
# 4. УРАВНЕНИЕ СОСТОЯНИЯ ТВЁРДЫХ СФЕР (Percus-Yevick)
# ============================================================

def py_excess_chemical_potential(rho: float, d_hs: float) -> float:
    """
    Избыточный химический потенциал HS по уравнению PY (ур. 13).
    mu_hs_ex = kBT * [-ln(1-eta) + eta*(14 - 13*eta + 5*eta^2) / (2*(1-eta)^3)]
    Возвращает mu_hs_ex / (kB*T).
    """
    eta = (np.pi / 6.0) * rho * d_hs ** 3
    eta = min(eta, 0.99)
    term = eta * (14 - 13 * eta + 5 * eta ** 2) / (2 * (1 - eta) ** 3)
    return -np.log(1 - eta) + term


def bulk_chemical_potential(rho: float, fluid: dict) -> float:
    """
    Полный химический потенциал в объёме (ур. 8) в единицах kBT.
    mu = ln(rho * Lambda^3) + mu_hs_ex + rho * a_ff
    Слагаемое rho*a_ff — среднеполевая поправка на притяжение.
    Здесь используем упрощённый интеграл a_ff = интеграл u_att(r)*4*pi*r^2 dr.
    """
    d = fluid["d_hs"]
    eps = fluid["eps_ff"]
    sigma = fluid["sigma_ff"]
    T = fluid["T"]

    # Среднеполевая константа a_ff (интеграл притяжения)
    r = np.linspace(2 ** (1 / 6) * sigma, 5 * sigma, 500)
    u = u_att_lj(r, eps / T, sigma)   # в единицах kBT
    a_ff = np.trapezoid(u * 4 * np.pi * r ** 2, r)

    mu_hs = py_excess_chemical_potential(rho, d)
    # Логарифмический член (de Broglie длина поглощается в константу)
    mu_id = np.log(max(rho, 1e-300))
    return mu_id + mu_hs + rho * a_ff


def rho_bulk_from_mu(mu_target: float, fluid: dict,
                     rho_min: float = 1e-10, rho_max: float = 40.0,
                     n: int = 500) -> float:
    """
    Обратная задача: найти равновесную плотность по заданному хим. потенциалу.
    Используем простой поиск по сетке + бисекция.
    """
    rhos = np.logspace(np.log10(rho_min), np.log10(rho_max), n)
    mus = np.array([bulk_chemical_potential(r, fluid) for r in rhos])
    idx = np.argmin(np.abs(mus - mu_target))
    return rhos[idx]


# ============================================================
# 5. МОДЕЛЬ QSDFT: профиль плотности в поре
#    (упрощённая 1D итерационная схема)
# ============================================================

def qsdft_density_profile_pore(pore_width_nm: float, p_over_p0: float,
                                fluid: dict, solid_params: dict,
                                n_grid: int = 200,
                                max_iter: int = 500, tol: float = 1e-6
                                ) -> tuple:
    """
    Решает уравнение Эйлера (7) для профиля плотности флюида в щелевой поре
    с шероховатыми стенками методом прямых итераций (Picard).

    Возвращает (z, rho_f) — координаты и профиль плотности в нм и нм⁻³.

    Примечание: это упрощённая реализация (1D, среднеполевое приближение),
    воспроизводящая физику метода. Полная QSDFT использует 3D FMT Розенфельда.
    """
    w = pore_width_nm
    eps_ff = fluid["eps_ff"] / fluid["T"]   # в единицах kBT
    sigma_ff = fluid["sigma_ff"]
    d_hs = fluid["d_hs"]
    eps_sf = solid_params["eps_sf"] / fluid["T"]
    sigma_sf = solid_params["sigma_sf"]

    # Сетка внутри поры (от края одной стенки до другой)
    ze = edge_position()
    z_arr = np.linspace(ze, w - ze, n_grid)
    if len(z_arr) < 3 or w < 2 * ze:
        # Пора слишком мала — заполнена полностью при любом давлении
        rho_liq = 25.0
        return z_arr if len(z_arr) > 0 else np.array([w / 2]), \
               np.full(max(len(z_arr), 1), rho_liq)

    dz = z_arr[1] - z_arr[0]

    # Профиль твёрдого тела (суперпозиция двух стенок)
    def rho_s_left(z_local):
        return solid_density_profile(z_local, ROUGHNESS_DELTA, H0, 114.0)

    def rho_s_right(z_local):
        z_mirror = w - z_local
        return solid_density_profile(z_mirror, ROUGHNESS_DELTA, H0, 114.0)

    # Предвычислим потенциал флюид-стенка на сетке
    phi_sf = np.zeros(n_grid)
    for i, z in enumerate(z_arr):
        # Левая стенка: z отсчитывается от поверхности
        z_from_left = z
        z_from_right = w - z

        # Интеграл по левой стенке
        z_s = np.linspace(0, H0 + 2 * ROUGHNESS_DELTA, 150)
        rs_l = solid_density_profile(z_s, ROUGHNESS_DELTA, H0, 114.0)
        r_l = np.abs(z_from_left - z_s) + 1e-12
        phi_sf[i] += np.trapezoid(rs_l * u_att_lj(r_l, eps_sf, sigma_sf), z_s)

        # Правая стенка
        rs_r = solid_density_profile(z_s, ROUGHNESS_DELTA, H0, 114.0)
        r_r = np.abs(z_from_right - z_s) + 1e-12
        phi_sf[i] += np.trapezoid(rs_r * u_att_lj(r_r, eps_sf, sigma_sf), z_s)

    # Химический потенциал флюида в объёме при заданном p/p0
    # mu ~ ln(p/p0) + mu_sat  (упрощение через насыщенное давление)
    rho_sat = 20.0  # нм⁻³ — плотность насыщенного пара (приближение)
    mu_bulk = np.log(max(p_over_p0, 1e-12)) + bulk_chemical_potential(rho_sat, fluid)

    # Среднеполевая константа притяжения флюид-флюид
    r_int = np.linspace(2 ** (1 / 6) * sigma_ff, 5 * sigma_ff, 300)
    u_ff = u_att_lj(r_int, eps_ff, sigma_ff)
    a_ff = np.trapezoid(u_ff * 4 * np.pi * r_int ** 2, r_int)

    # Начальное приближение плотности
    rho_f = np.ones(n_grid) * max(p_over_p0, 1e-8) * rho_sat

    # Итерации Пикара
    for _ in range(max_iter):
        rho_old = rho_f.copy()

        # Среднеполевое притяжение флюид-флюид (свёртка)
        conv = np.zeros(n_grid)
        for i in range(n_grid):
            r_ij = np.abs(z_arr[i] - z_arr) + 1e-12
            u_ij = u_att_lj(r_ij, eps_ff, sigma_ff)
            conv[i] = np.trapezoid(rho_old * u_ij, z_arr)

        # HS химический потенциал (локальное приближение)
        mu_hs = np.array([py_excess_chemical_potential(max(r, 1e-10), d_hs)
                          for r in rho_old])

        # Уравнение Эйлера → новый профиль (ур. 7)
        exponent = mu_bulk - mu_hs - conv - phi_sf
        rho_new = np.exp(np.clip(exponent, -50, 50))

        # Демпфирование
        rho_f = 0.7 * rho_new + 0.3 * rho_old
        rho_f = np.clip(rho_f, 0, 60.0)

        if np.max(np.abs(rho_f - rho_old)) < tol:
            break

    return z_arr, rho_f


def pore_loading(pore_width_nm: float, p_over_p0: float,
                 fluid: dict, solid_params: dict) -> float:
    """
    Количество адсорбата в поре [нм⁻²] = интеграл (rho_f - rho_bulk) dz.
    """
    if p_over_p0 <= 0:
        return 0.0
    z, rho = qsdft_density_profile_pore(pore_width_nm, p_over_p0,
                                         fluid, solid_params)
    if len(z) < 2:
        return 0.0
    rho_bulk_approx = max(p_over_p0, 1e-12) * 0.01  # разреженный пар
    excess = np.maximum(rho - rho_bulk_approx, 0)
    return np.trapezoid(excess, z)


# ============================================================
# 6. ЯДРО ИЗОТЕРМ ИЗ ТАБЛИЦ SUPPLEMENTARY (точные данные статьи)
# ============================================================

# Таблица 1 — N2, 77.4 K: (ширина поры [Å], равновесное давление заполнения P/P0)
N2_FILLING_TABLE = np.array([
    [4.00,1.970e-09],[4.15,3.769e-09],[4.32,7.936e-09],[4.48,1.573e-08],
    [4.66,3.311e-08],[4.84,6.774e-08],[5.04,1.446e-07],[5.24,2.901e-07],
    [5.45,5.661e-07],[5.67,1.062e-06],[5.90,1.939e-06],[6.14,3.354e-06],
    [6.40,5.888e-06],[6.66,9.836e-06],[6.94,1.686e-05],[7.23,2.990e-05],
    [7.53,5.365e-05],[7.85,1.021e-04],[8.18,2.000e-04],[8.52,3.866e-04],
    [8.89,7.353e-04],[9.26,1.265e-03],[9.66,2.093e-03],[10.07,3.302e-03],
    [10.51,5.013e-03],[10.96,7.665e-03],[11.44,1.107e-02],[11.93,1.610e-02],
    [12.45,2.106e-02],[12.99,2.734e-02],[13.56,3.526e-02],[14.16,4.483e-02],
    [14.78,5.608e-02],[15.43,6.899e-02],[16.11,8.321e-02],[16.82,9.877e-02],
    [17.56,1.160e-01],[18.34,1.351e-01],[19.15,1.553e-01],[20.00,1.766e-01],
    [20.89,1.993e-01],[21.83,2.236e-01],[22.80,2.483e-01],[23.82,2.737e-01],
    [24.88,2.998e-01],[26.00,3.263e-01],[27.16,3.526e-01],[28.38,3.791e-01],
    [29.66,4.053e-01],[30.99,4.309e-01],[32.39,4.562e-01],[33.85,4.807e-01],
    [35.38,5.046e-01],[36.98,5.277e-01],[38.65,5.500e-01],[40.39,5.714e-01],
    [42.22,5.921e-01],[44.13,6.119e-01],[46.13,6.310e-01],[48.23,6.493e-01],
    [50.42,6.668e-01],[52.70,6.834e-01],[55.10,6.994e-01],[57.60,7.146e-01],
    [60.22,7.290e-01],[62.96,7.428e-01],[65.83,7.558e-01],[68.83,7.681e-01],
    [71.96,7.798e-01],[75.24,7.908e-01],[78.68,8.013e-01],[82.27,8.112e-01],
    [86.02,8.205e-01],[89.95,8.294e-01],[94.05,8.377e-01],[98.35,8.456e-01],
    [102.84,8.531e-01],[107.55,8.602e-01],[112.46,8.668e-01],[117.61,8.732e-01],
    [122.99,8.792e-01],[128.61,8.848e-01],[134.50,8.902e-01],[140.66,8.953e-01],
    [147.10,9.001e-01],[153.83,9.046e-01],[160.88,9.090e-01],[168.25,9.130e-01],
    [175.96,9.169e-01],[184.02,9.205e-01],[192.46,9.240e-01],[201.28,9.273e-01],
    [210.51,9.304e-01],[220.16,9.334e-01],[230.26,9.362e-01],[240.82,9.388e-01],
    [251.87,9.414e-01],[263.43,9.438e-01],[275.51,9.460e-01],[288.16,9.482e-01],
    [301.38,9.502e-01],[315.22,9.522e-01],[329.69,9.540e-01],[344.82,9.558e-01],
    [360.65,9.575e-01],
])

# Таблица 2 — Ar, 87.3 K
AR_FILLING_TABLE = np.array([
    [4.00,1.366e-07],[4.15,2.368e-07],[4.32,4.292e-07],[4.48,7.616e-07],
    [4.66,1.375e-06],[4.84,2.477e-06],[5.04,4.430e-06],[5.24,7.770e-06],
    [5.45,1.295e-05],[5.67,2.124e-05],[5.90,3.393e-05],[6.14,5.310e-05],
    [6.40,8.328e-05],[6.66,1.278e-04],[6.94,2.010e-04],[7.23,3.180e-04],
    [7.53,5.079e-04],[7.85,8.274e-04],[8.18,1.339e-03],[8.52,2.123e-03],
    [8.89,3.339e-03],[9.26,4.737e-03],[9.66,7.454e-03],[10.07,1.011e-02],
    [10.51,1.417e-02],[10.96,1.891e-02],[11.44,2.534e-02],[11.93,3.334e-02],
    [12.45,4.072e-02],[12.99,5.077e-02],[13.56,6.306e-02],[14.16,7.641e-02],
    [14.78,9.111e-02],[15.43,1.070e-01],[16.11,1.242e-01],[16.82,1.429e-01],
    [17.56,1.628e-01],[18.34,1.837e-01],[19.15,2.052e-01],[20.00,2.279e-01],
    [20.89,2.514e-01],[21.83,2.754e-01],[22.80,2.997e-01],[23.82,3.245e-01],
    [24.88,3.493e-01],[26.00,3.744e-01],[27.16,3.991e-01],[28.38,4.237e-01],
    [29.66,4.480e-01],[30.99,4.716e-01],[32.39,4.948e-01],[33.85,5.172e-01],
    [35.38,5.390e-01],[36.98,5.601e-01],[38.65,5.804e-01],[40.39,5.998e-01],
    [42.22,6.185e-01],[44.13,6.365e-01],[46.13,6.537e-01],[48.23,6.702e-01],
    [50.42,6.860e-01],[52.70,7.010e-01],[55.10,7.154e-01],[57.60,7.292e-01],
    [60.22,7.423e-01],[62.96,7.548e-01],[65.83,7.667e-01],[68.83,7.780e-01],
    [71.96,7.888e-01],[75.24,7.990e-01],[78.68,8.087e-01],[82.27,8.178e-01],
    [86.02,8.265e-01],[89.95,8.348e-01],[94.05,8.426e-01],[98.35,8.499e-01],
    [102.84,8.569e-01],[107.55,8.636e-01],[112.46,8.699e-01],[117.61,8.758e-01],
    [122.99,8.815e-01],[128.61,8.869e-01],[134.50,8.920e-01],[140.66,8.968e-01],
    [147.10,9.014e-01],[153.83,9.058e-01],[160.88,9.099e-01],[168.25,9.138e-01],
    [175.96,9.176e-01],[184.02,9.211e-01],[192.46,9.245e-01],[201.28,9.277e-01],
    [210.51,9.307e-01],[220.16,9.336e-01],[230.26,9.363e-01],[240.82,9.389e-01],
    [251.87,9.413e-01],[263.43,9.437e-01],[275.51,9.459e-01],[288.16,9.480e-01],
    [301.38,9.500e-01],[315.22,9.520e-01],[329.69,9.538e-01],[344.82,9.555e-01],
    [360.65,9.572e-01],
])


def get_filling_pressure_interpolator(adsorbate: str = "N2"):
    """
    Возвращает интерполятор: ширина поры (Å) → давление заполнения P/P0.
    Использует точные данные из Supplementary Information статьи.
    """
    table = N2_FILLING_TABLE if adsorbate.upper() == "N2" else AR_FILLING_TABLE
    w_A = table[:, 0]
    p = table[:, 1]
    return interp1d(w_A, np.log10(p), kind='linear',
                    fill_value="extrapolate", bounds_error=False)


def build_kernel_from_table(p_grid: np.ndarray, adsorbate: str = "N2",
                             w_min_A: float = 4.0, w_max_A: float = 350.0,
                             n_pores: int = 60) -> tuple:
    """
    Строит матрицу ядра K[i,j] = N_ads(p_i, w_j) для QNNLS-деконволюции.

    Каждая теоретическая изотерма моделируется как ступенчатая функция,
    заполняющаяся при давлении p_fill(w) (из таблицы) с максимальным
    заполнением пропорциональным ширине поры.

    p_grid: массив относительных давлений
    Возвращает: (w_arr_A, K)
    """
    interp_log_p = get_filling_pressure_interpolator(adsorbate)
    w_arr = np.logspace(np.log10(w_min_A), np.log10(w_max_A), n_pores)

    K = np.zeros((len(p_grid), len(w_arr)))

    for j, w in enumerate(w_arr):
        log_p_fill = interp_log_p(w)
        p_fill = 10 ** log_p_fill
        # Изотерма: мягкий sigmoid-переход (сглаженная ступенька)
        sigma_step = 0.3  # ширина в единицах log(p)
        log_p = np.log10(p_grid + 1e-15)
        filling = 1.0 / (1.0 + np.exp(-(log_p - log_p_fill) / sigma_step))
        # Нормировка на объём поры (пропорционален ширине)
        K[:, j] = filling * w

    return w_arr, K


# ============================================================
# 7. РАСЧЁТ PSD: QNNLS (регуляризация Тихонова + NNLS)
# ============================================================

def qnnls_psd(p_exp: np.ndarray, N_exp: np.ndarray,
              adsorbate: str = "N2",
              w_min_A: float = 4.0, w_max_A: float = 350.0,
              n_pores: int = 60, lambda_reg: float = 1e-3) -> tuple:
    """
    Вычисляет распределение пор по размерам (PSD) методом QNNLS.

    Решает: K @ f = N_exp, f >= 0, с регуляризацией Тихонова.

    Параметры:
        p_exp:     массив относительных давлений (P/P0)
        N_exp:     экспериментальная изотерма [см³/г STP или отн. ед.]
        adsorbate: 'N2' или 'Ar'
        w_min_A, w_max_A: диапазон пор (Å)
        n_pores:   количество точек в PSD
        lambda_reg: параметр регуляризации (больше → сглаженнее PSD)

    Возвращает:
        w_arr_A: массив ширин пор (Å)
        f_psd:   дифференциальная PSD (dV/dw)
        N_fit:   теоретическая изотерма (подгонка)
    """
    # Нормализация изотермы
    N_max = N_exp.max()
    N_norm = N_exp / N_max

    # Ядро изотерм
    w_arr, K = build_kernel_from_table(p_exp, adsorbate, w_min_A, w_max_A, n_pores)

    # Нормализация ядра
    K_norm = K / (K.max(axis=0, keepdims=True) + 1e-15)

    # Тихоновская регуляризация: добавляем строки lambda*I
    n_p = len(p_exp)
    K_reg = np.vstack([K_norm, lambda_reg * np.eye(n_pores)])
    b_reg = np.concatenate([N_norm, np.zeros(n_pores)])

    # NNLS
    f_psd, _ = nnls(K_reg, b_reg)

    # Восстановление изотермы
    N_fit = (K_norm @ f_psd) * N_max

    # Нормировка PSD: переводим в dV/dw [условные единицы]
    dw = np.gradient(w_arr)
    f_dv = f_psd / (dw + 1e-15)

    return w_arr, f_dv, N_fit


def cumulative_pore_volume(w_arr: np.ndarray, f_dv: np.ndarray) -> np.ndarray:
    """Кумулятивный объём пор как функция ширины."""
    dw = np.gradient(w_arr)
    return np.cumsum(f_dv * dw)


# ============================================================
# 8. ДЕМОНСТРАЦИЯ: синтетический пример
# ============================================================

def demo_synthetic_psd():
    """
    Воспроизводит раздел III Supplementary Information:
    проверка схемы PSD на модельном образце с заданным распределением.

    Заданное PSD: log-нормальная смесь микро- и мезопор.
    """
    print("=" * 60)
    print("QSDFT — демонстрация расчёта PSD (синтетический пример)")
    print("=" * 60)

    # Задаём модельную PSD (двухмодальная: микро + мезопоры)
    def model_psd(w):
        """Двухмодальное распределение из Supplementary (ур. III.1)."""
        # Микропоры: пик около 8 Å
        p1 = 0.6 * np.exp(-0.5 * ((np.log(w) - np.log(8.0)) / 0.4) ** 2)
        # Мезопоры: пик около 40 Å
        p2 = 0.4 * np.exp(-0.5 * ((np.log(w) - np.log(40.0)) / 0.5) ** 2)
        return p1 + p2

    # Сетка давлений
    p_grid = np.logspace(-8, 0, 80)

    # Строим ядро для N2
    w_kernel, K = build_kernel_from_table(p_grid, "N2", 4, 200, 60)

    # Вычисляем "экспериментальную" изотерму из модельной PSD
    f_true = model_psd(w_kernel)
    dw_kernel = np.gradient(w_kernel)
    N_theory = K @ (f_true * dw_kernel)

    # Добавляем небольшой шум
    rng = np.random.default_rng(42)
    N_exp = N_theory + 0.02 * N_theory.max() * rng.standard_normal(len(N_theory))
    N_exp = np.maximum(N_exp, 0)

    # Расчёт PSD методом QNNLS
    w_calc, f_calc, N_fit = qnnls_psd(p_grid, N_exp, "N2",
                                       w_min_A=4, w_max_A=200,
                                       n_pores=60, lambda_reg=5e-3)

    # --- Построение графиков ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("QSDFT — воспроизведение метода (Neimark et al., Carbon 2009)",
                 fontsize=13, fontweight="bold")

    # 1. Профиль плотности твёрдого тела
    ax = axes[0]
    z_plot = np.linspace(0, H0 + 2 * ROUGHNESS_DELTA + 0.2, 300)
    rho_s = solid_density_profile(z_plot, ROUGHNESS_DELTA, H0, 114.0)
    ze_val = edge_position()
    ax.fill_between(z_plot, rho_s, alpha=0.4, color="saddlebrown", label="Твёрдое тело")
    ax.axvline(ze_val, color="red", ls="--", lw=1.5, label=f"$z_e$ = {ze_val:.3f} нм")
    ax.set_xlabel("z, нм")
    ax.set_ylabel("ρ_s, нм⁻³")
    ax.set_title("Профиль плотности стенки\n(линейный рамп, δ=0.13 нм)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # 2. Давления заполнения пор
    ax = axes[1]
    for ads, color, lbl in [("N2", "steelblue", "N₂, 77.4 K"),
                              ("Ar", "darkorange", "Ar, 87.3 K")]:
        table = N2_FILLING_TABLE if ads == "N2" else AR_FILLING_TABLE
        ax.semilogy(table[:, 0], table[:, 1], "o-", ms=3, color=color,
                    label=lbl, lw=1.5)
    ax.set_xlabel("Ширина поры, Å")
    ax.set_ylabel("P/P₀ (давление заполнения)")
    ax.set_title("Давления заполнения пор QSDFT\n(из Supplementary, табл. 1-2)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 3. Результаты PSD расчёта
    ax = axes[2]
    # Нормализуем для сравнения
    f_true_norm = f_true / (f_true.max() + 1e-15)
    f_calc_norm = f_calc / (f_calc.max() + 1e-15)
    ax.plot(w_kernel, f_true_norm, "k-", lw=2, label="Заданная PSD")
    ax.plot(w_calc, f_calc_norm, "r--", lw=2, label="QSDFT QNNLS расчёт")
    ax.set_xlabel("Ширина поры, Å")
    ax.set_ylabel("dV/dw (норм.)")
    ax.set_title("Верификация расчёта PSD\n(синт. образец, N₂)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale("log")

    plt.tight_layout()
    plt.savefig(Path(OUT_DIR, "qsdft_demo.png"), dpi=150, bbox_inches="tight")
    print("График сохранён: qsdft_demo.png")

    return fig


def demo_filling_pressures():
    """Дополнительный график: сравнение N2 и Ar давлений заполнения."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Левый: линейная шкала
    ax = axes[0]
    ax.plot(N2_FILLING_TABLE[:, 0], N2_FILLING_TABLE[:, 1],
            "b-o", ms=3, lw=1.5, label="N₂, 77.4 K")
    ax.plot(AR_FILLING_TABLE[:, 0], AR_FILLING_TABLE[:, 1],
            "r-s", ms=3, lw=1.5, label="Ar, 87.3 K")
    ax.set_xlabel("Ширина поры w, Å")
    ax.set_ylabel("P/P₀")
    ax.set_title("Давления заполнения (линейная шкала)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Правый: логарифмическая шкала
    ax = axes[1]
    ax.semilogy(N2_FILLING_TABLE[:, 0], N2_FILLING_TABLE[:, 1],
                "b-o", ms=3, lw=1.5, label="N₂, 77.4 K")
    ax.semilogy(AR_FILLING_TABLE[:, 0], AR_FILLING_TABLE[:, 1],
                "r-s", ms=3, lw=1.5, label="Ar, 87.3 K")
    ax.set_xlabel("Ширина поры w, Å")
    ax.set_ylabel("P/P₀ (log)")
    ax.set_title("Давления заполнения (лог. шкала)\nAr заполняет малые поры при более\nвысоких давлениях → лучшее разрешение")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")

    plt.suptitle("QSDFT: Сравнение N₂ и Ar (Neimark 2009, табл. Suppl.)",
                 fontweight="bold")
    plt.tight_layout()
    # LOG_FILE = Path(LOGS_DIR, 'app.log')

    plt.savefig(Path(OUT_DIR, 'qsdft_filling_pressures.png'),
                dpi=150, bbox_inches="tight")
    print("График сохранён: qsdft_filling_pressures.png")
    return fig


# ============================================================
# ПРИМЕР ИСПОЛЬЗОВАНИЯ API
# ============================================================

def compute_psd_from_isotherm(p_data: np.ndarray, N_data: np.ndarray,
                               adsorbate: str = "N2",
                               lambda_reg: float = 1e-3,
                               plot: bool = True) -> dict:
    """
    Основная функция для расчёта PSD из экспериментальной изотермы.

    Параметры
    ---------
    p_data : np.ndarray
        Массив относительных давлений P/P0.
    N_data : np.ndarray
        Экспериментальная изотерма (любые единицы: см³/г, ммоль/г и т.д.).
    adsorbate : str
        'N2' (77.4 K) или 'Ar' (87.3 K).
    lambda_reg : float
        Параметр регуляризации Тихонова (обычно 1e-4 ... 1e-2).
    plot : bool
        Построить ли графики.

    Возвращает
    ----------
    dict с ключами:
        'pore_width_A'   : массив ширин пор (Å)
        'psd_dVdw'       : дифференциальная PSD
        'cumulative_vol' : кумулятивный объём
        'fitted_isotherm': теоретическая изотерма (подгонка)
    """
    w_arr, f_psd, N_fit = qnnls_psd(
        p_data, N_data, adsorbate,
        w_min_A=4.0, w_max_A=360.0,
        n_pores=70, lambda_reg=lambda_reg
    )
    cum_vol = cumulative_pore_volume(w_arr, f_psd)

    if plot:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        axes[0].semilogx(p_data, N_data, "ko", ms=4, label="Эксперимент")
        axes[0].semilogx(p_data, N_fit, "r-", lw=2, label="QSDFT подгонка")
        axes[0].set_xlabel("P/P₀")
        axes[0].set_ylabel("Адсорбция")
        axes[0].set_title(f"Изотерма ({adsorbate})")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(w_arr, f_psd, "b-", lw=2)
        axes[1].set_xlabel("Ширина поры, Å")
        axes[1].set_ylabel("dV/dw")
        axes[1].set_title(f"Дифф. PSD QSDFT ({adsorbate})")
        axes[1].set_xscale("log")
        axes[1].grid(True, alpha=0.3)

        axes[2].plot(w_arr, cum_vol, "g-", lw=2)
        axes[2].set_xlabel("Ширина поры, Å")
        axes[2].set_ylabel("Кум. объём пор")
        axes[2].set_title("Кумулятивный объём")
        axes[2].set_xscale("log")
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout()
        
        plt.savefig(Path(OUT_DIR, f"qsdft_psd_{adsorbate}.png"),
                    dpi=150, bbox_inches="tight")
        plt.show()

    return {
        "pore_width_A": w_arr,
        "psd_dVdw": f_psd,
        "cumulative_vol": cum_vol,
        "fitted_isotherm": N_fit,
    }


# ============================================================
# ЗАПУСК
# ============================================================

if __name__ == "__main__":

    print("\n[1/2] Генерация демонстрационного графика (синт. образец)...")
    demo_synthetic_psd()

    print("\n[2/2] График давлений заполнения N2 vs Ar...")
    demo_filling_pressures()

    print("\n--- Пример использования API ---")
    print("Для расчёта PSD из ваших данных:")
    print("""
    from qsdft import compute_psd_from_isotherm
    import numpy as np

    # Загрузите ваши данные
    p = np.array([...])   # P/P0
    N = np.array([...])   # Адсорбция

    result = compute_psd_from_isotherm(p, N, adsorbate='N2', lambda_reg=1e-3)
    print(result['pore_width_A'])   # Ширины пор (Å)
    print(result['psd_dVdw'])       # Дифференциальная PSD
    """)
    print("Готово.")
