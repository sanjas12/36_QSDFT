# import pygaps
# import pygaps.characterisation


# pygaps.characterisation.psd_kernel.psd_dft()
import numpy as np
import matplotlib.pyplot as plt
# from scipy.integrate import simps
from scipy.integrate import simpson  # Новый аналог

# Параметры модели
kB = 1.380649e-23  # Дж/К
T = 87.3           # Температура (K)
NA = 6.022e23      # число Авогадро

# LJ параметры для Ar–Ar
eps_ff = 111.95 * kB     # Дж
sigma_ff = 0.3358e-9     # м

# LJ параметры для C–Ar
eps_sf = 162.18 * kB     # Дж
sigma_sf = 0.2595e-9     # м

# Параметры модели пор
roughness = 0.13e-9      # шероховатость, м
h0 = 0.68e-9             # толщина стенки поры, м
z_max = 5e-9             # предел по оси z
n_points = 1000          # шаг сетки

z = np.linspace(0, z_max, n_points)

# Профиль плотности твёрдой фазы (стенка)
def rho_solid(z):
    d = roughness
    q0 = 0.114e30  # число/м³
    profile = np.zeros_like(z)
    for i, zi in enumerate(z):
        if zi < h0:
            profile[i] = q0
        elif h0 <= zi <= h0 + 2*d:
            profile[i] = q0 * (1 - (zi - h0) / (2*d))
        else:
            profile[i] = 0
    return profile

rho_s = rho_solid(z)

# Потенциал взаимодействия Lennard-Jones между газом и стенкой
def u_sf(z, rho_s, sigma, epsilon):
    r_cut = 1.5e-9  # отсечка
    dz = z[1] - z[0]
    u = np.zeros_like(z)
    for i in range(len(z)):
        integrand = np.zeros_like(z)
        for j in range(len(z)):
            r = np.abs(z[i] - z[j])
            if r > 0 and r < r_cut:
                lj = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
                integrand[j] = lj * rho_s[j]
        u[i] = simpson(integrand, z)
    return u

U_ext = u_sf(z, rho_s, sigma_sf, eps_sf)

# Расчёт плотности жидкости через Больцмана
def fluid_density(P_over_P0):
    rho_bulk = 25e3 / 40e-3 * NA  # ≈ молекул/м³ (при 87.3K и 1 атм)
    mu = kB * T * np.log(P_over_P0)
    rho_f = rho_bulk * np.exp(-(U_ext - mu) / (kB * T))
    return rho_f

# Расчёт изотермы
pressures = np.logspace(-6, 0, 100)
adsorption = []

for p in pressures:
    rho_f = fluid_density(p)
    N_ads = simpson(rho_f, z)  # интеграл по профилю плотности
    adsorption.append(N_ads / NA * 1e3)  # ммоль/м²

# Построение изотермы
plt.figure(figsize=(7,5))
plt.semilogx(pressures, adsorption)
plt.xlabel("P/P₀")
plt.ylabel("Адсорбция (ммоль/м²)")
plt.title("QSDFT-подобная изотерма аргона при 87.3K")
plt.grid(True)
plt.tight_layout()
plt.show()
