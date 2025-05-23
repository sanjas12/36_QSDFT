# mail.py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import simpson
from typing import Tuple, List, Optional, Dict, Any
from AdsorptionDataPlotter import AdsorptionDataPlotter


class QSDFT:
    # Physical constants (class level)
    KB: float = 1.380649e-23  # Boltzmann constant [J/K]
    NA: float = 6.022e23  # Avogadro's number

    def __init__(
        self,
        temperature: float = 87.3,
        eps_ff_k: float = 111.95,
        sigma_ff: float = 0.3358e-9,
        eps_sf_k: float = 162.18,
        sigma_sf: float = 0.2595e-9,
        roughness: float = 0.13e-9,
        h0: float = 0.68e-9,
        z_max: float = 5e-9,
        q0: float = 0.114e30,
        r_cut: float = 1.5e-9,
    ) -> None:
        """
        Initialize QSDFT adsorption calculator with physical parameters.

        Parameters:
            temperature: System temperature in Kelvin [K]
            eps_ff_k: Fluid-fluid interaction LJ ε/k_B parameter [K]
            sigma_ff: Fluid-fluid interaction LJ σ parameter [m]
            eps_sf_k: Solid-fluid interaction LJ ε/k_B parameter [K]
            sigma_sf: Solid-fluid interaction LJ σ parameter [m]
            roughness: Surface roughness [m]
            h0: Wall thickness [m]
            z_max: Maximum calculation depth in z-direction [m]
            q0: Carbon atom density [m^-3]
            r_cut: LJ potential cutoff radius [m]
        """
        # System parameters
        self.T: float = temperature
        self.eps_ff: float = eps_ff_k * self.KB
        self.sigma_ff: float = sigma_ff
        self.eps_sf: float = eps_sf_k * self.KB
        self.sigma_sf: float = sigma_sf
        self.roughness: float = roughness
        self.h0: float = h0
        self.z_max: float = z_max
        self.q0: float = q0
        self.r_cut: float = r_cut

        # Initialize state variables
        self.z: np.ndarray = np.array([])
        self.rho_s: np.ndarray = np.array([])
        self.U_ext: np.ndarray = np.array([])
        self.adsorption: List[float] = []

        # Load experimental data
        self._load_experimental_data()

    def _load_experimental_data(self) -> None:
        """Load experimental pressure data from AdsorptionDataPlotter."""
        try:
            plotter = AdsorptionDataPlotter()
            self.pressures: np.ndarray = plotter.get_all_data()["Pressure"].to_numpy()
            self.n_points: int = len(self.pressures)
        except Exception as e:
            raise RuntimeError(f"Failed to load experimental data: {str(e)}")

    def initialize_arrays(self) -> None:
        """Initialize spatial grid and solid density profile."""
        self.z = np.linspace(0, self.z_max, self.n_points)
        self.rho_s = self._calculate_solid_density(self.z)

    def _calculate_solid_density(self, z: np.ndarray) -> np.ndarray:
        """
        Calculate solid wall density profile using step function with roughness.

        Args:
            z: Spatial coordinates array [m]

        Returns:
            Density profile of solid wall [m^-3]
        """
        profile = np.zeros_like(z)
        d = self.roughness

        # Vectorized operations for better performance
        mask1 = z < self.h0
        mask2 = (self.h0 <= z) & (z < self.h0 + 2 * d)

        profile[mask1] = self.q0
        profile[mask2] = self.q0 * (1 - (z[mask2] - self.h0) / (2 * d))

        return profile

    def _calculate_external_potential(
        self, z: np.ndarray, rho_s: np.ndarray
    ) -> np.ndarray:
        """
        Calculate Lennard-Jones potential for gas-solid interaction.

        Args:
            z: Spatial coordinates [m]
            rho_s: Solid density profile [m^-3]

        Returns:
            External potential profile [J]
        """
        dz = z[1] - z[0]
        u = np.zeros_like(z)

        for i, zi in enumerate(z):
            r = np.abs(zi - z)
            mask = (r > 0) & (r < self.r_cut)

            if not np.any(mask):
                continue

            r_valid = r[mask]
            rho_valid = rho_s[mask]

            # Vectorized LJ potential calculation
            sigma_r = self.sigma_sf / r_valid
            lj = 4 * self.eps_sf * (sigma_r**12 - sigma_r**6)
            u[i] = simpson(lj * rho_valid, z[mask])

        return u

    def calculate_adsorption_isotherm(self) -> None:
        """Calculate adsorption isotherm over pressure range."""
        self.initialize_arrays()
        self.U_ext = self._calculate_external_potential(self.z, self.rho_s)

        # Bulk density calculation
        rho_bulk = (25e3 / 40e-3) * self.NA  # molecules/m³

        # Calculate adsorption for each pressure
        self.adsorption = []
        for p in self.pressures:
            mu = self.KB * self.T * np.log(max(p, 1e-12))  # Chemical potential

            # Calculate density profile
            exponent = np.clip(-(self.U_ext - mu) / (self.KB * self.T), -100, 100)
            rho_f = rho_bulk * np.exp(exponent)

            # Calculate adsorbed amount
            N_ads = simpson(rho_f, self.z)  # molecules/m²
            self.adsorption.append(N_ads / self.NA * 1e3)  # type: ignore # mmol/m²

    def save_results(self, filename: str = "data_out.csv") -> None:
        """
        Save calculated adsorption data to CSV file.

        Args:
            filename: Output file path

        Raises:
            ValueError: If no adsorption data is available
            IOError: If file writing fails
        """
        if not self.adsorption:
            raise ValueError(
                "No adsorption data to save. Run calculate_adsorption_isotherm() first."
            )

        data = np.column_stack((self.pressures, self.adsorption))
        try:
            with open(filename, "w", encoding="utf-8") as f:
                np.savetxt(
                    f,
                    data,
                    header="Pressure [P/P0], Adsorption [mmol/m2]",  # Removed ² character
                    fmt="%.4e",
                    delimiter=",",
                )
        except IOError as e:
            raise IOError(f"Failed to save results: {str(e)}")

    def plot_isotherm(self, **kwargs: Any) -> plt.Figure:  # type: ignore
        """
        Plot calculated adsorption isotherm.

        Args:
            **kwargs: Additional plotting parameters

        Returns:
            matplotlib Figure object
        """
        if not self.adsorption:
            raise ValueError(
                "No adsorption data to plot. Run calculate_adsorption_isotherm() first."
            )

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.semilogx(self.pressures, self.adsorption, **kwargs)

        ax.set_xlabel("Relative pressure P/P₀")
        ax.set_ylabel("Adsorption (mmol/m²)")
        ax.set_title(f"Adsorption isotherm at {self.T} K (QSDFT model)")
        ax.grid(True)

        plt.tight_layout()
        return fig


if __name__ == "__main__":
    try:
        calculator = QSDFT()
        calculator.calculate_adsorption_isotherm()
        calculator.save_results()

        fig = calculator.plot_isotherm(label="QSDFT model", linewidth=2)
        plt.show()

    except Exception as e:
        print(f"Error: {str(e)}")
