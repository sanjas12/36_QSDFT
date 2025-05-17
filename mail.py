import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

# Физические константы
kB = 1.380649e-23       # Дж/К
T = 87.3                # Температура в K
NA = 6.022e23           # число Авогадро

# Lennard-Jones параметры для Ar–Ar и C–Ar
eps_ff = 111.95 * kB         # Дж
sigma_ff = 0.3358e-9         # м
eps_sf = 162.18 * kB         # Дж
sigma_sf = 0.2595e-9         # м

# Параметры пористой модели
roughness = 0.13e-9          # шероховатость поверхности, м
h0 = 0.68e-9                 # толщина стенки, м
z_max = 5e-9                 # глубина расчёта по z, м
n_points = 1000             # число точек дискретизации

z = np.linspace(0, z_max, n_points)

# Плотность твёрдой стенки
def rho_solid(z):
    q0 = 0.114e30  # плотность атомов углерода, м^-3
    d = roughness
    profile = np.zeros_like(z)
    for i, zi in enumerate(z):
        if zi < h0:
            profile[i] = q0
        elif h0 <= zi < h0 + 2*d:
            profile[i] = q0 * (1 - (zi - h0) / (2*d))
        else:
            profile[i] = 0
    return profile

rho_s = rho_solid(z)

# Lennard-Jones потенциал взаимодействия газ–твёрдое тело
def u_sf(z, rho_s, sigma, epsilon):
    r_cut = 1.5e-9
    dz = z[1] - z[0]
    u = np.zeros_like(z)
    for i in range(len(z)):
        r = np.abs(z[i] - z)
        mask = (r > 0) & (r < r_cut)
        r_valid = r[mask]
        rho_valid = rho_s[mask]
        lj = 4 * epsilon * ((sigma / r_valid)**12 - (sigma / r_valid)**6)
        integrand = lj * rho_valid
        u[i] = simpson(integrand, z[mask])
    return u

U_ext = u_sf(z, rho_s, sigma_sf, eps_sf)

# Давление в логарифмическом масштабе (P/P0)
pressures = np.logspace(-4, 0, 100)  # от 10^-4 до 1

adsorption = []

# Расчёт изотермы
for p in pressures:
    # Расчёт химпотенциала и плотности жидкости
    rho_bulk = (25e3 / 40e-3) * NA  # молекул/м³ (примерно)
    mu = kB * T * np.log(max(p, 1e-12))  # защита от log(0)

    # Расчёт экспоненты с обрезкой
    exponent = -(U_ext - mu) / (kB * T)
    exponent = np.clip(exponent, -100, 100)
    rho_f = rho_bulk * np.exp(exponent)

    # Интеграл по плотности жидкости — количество адсорбированного вещества
    N_ads = simpson(rho_f, z)  # молекул/м²
    adsorption.append(N_ads / NA * 1e3)  # ммоль/м²


print("Adsorption data:", adsorption[:])  # первые 5 значений
print("Pressures:", pressures[:])

with open("adsorption_data.txt", "w") as f:
    for p, N_ads in zip(pressures, adsorption):
        f.write(f"{p:.4e} {N_ads:.4e}\n")


# Построение графика изотермы
plt.figure(figsize=(7, 5))
plt.semilogx(pressures, adsorption, label='QSDFT-подобная модель')
plt.xlabel("Относительное давление P/P₀")
plt.ylabel("Адсорбция (ммоль/м²)")
plt.title("Изотерма адсорбции аргона при 87.3 K (упрощённая QSDFT)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
