# mail.py
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from typing import Tuple, List, Optional

class QSDFT:
    def __init__(self) -> None:
        """
        Initialize the adsorption calculator with physical constants and parameters.
        """
        # Physical constants
        self.kB: float = 1.380649e-23       # J/K
        self.T: float = 87.3                # Temperature in K
        self.NA: float = 6.022e23           # Avogadro's number
        
        # Lennard-Jones parameters for Ar-Ar and C-Ar
        self.eps_ff: float = 111.95 * self.kB   # J
        self.sigma_ff: float = 0.3358e-9        # m
        self.eps_sf: float = 162.18 * self.kB   # J
        self.sigma_sf: float = 0.2595e-9        # m
        
        # Porous model parameters
        self.roughness: float = 0.13e-9         # surface roughness, m
        self.h0: float = 0.68e-9                # wall thickness, m
        self.z_max: float = 5e-9                # calculation depth in z, m
        self.n_points: int = 1000               # number of discretization points
        
        # Initialize arrays
        self.z: Optional[np.ndarray] = None
        self.rho_s: Optional[np.ndarray] = None
        self.U_ext: Optional[np.ndarray] = None
        self.pressures: Optional[np.ndarray] = None
        self.adsorption: Optional[List[float]] = None
        
    def initialize_arrays(self) -> None:
        """
        Initialize the spatial grid and solid density profile.
        """
        self.z = np.linspace(0, self.z_max, self.n_points)
        self.rho_s = self._calculate_solid_density(self.z)
        
    def _calculate_solid_density(self, z: np.ndarray) -> np.ndarray:
        """
        Calculate the solid wall density profile.
        
        Args:
            z: Spatial coordinates array
            
        Returns:
            Density profile of the solid wall
        """
        q0: float = 0.114e30  # carbon atom density, m^-3
        d: float = self.roughness
        profile = np.zeros_like(z)
        for i, zi in enumerate(z):
            if zi < self.h0:
                profile[i] = q0
            elif self.h0 <= zi < self.h0 + 2*d:
                profile[i] = q0 * (1 - (zi - self.h0) / (2*d))
            else:
                profile[i] = 0
        return profile
    
    def _calculate_external_potential(self, 
                                   z: np.ndarray, 
                                   rho_s: np.ndarray, 
                                   sigma: float, 
                                   epsilon: float) -> np.ndarray:
        """
        Calculate Lennard-Jones potential for gas-solid interaction.
        
        Args:
            z: Spatial coordinates
            rho_s: Solid density profile
            sigma: LJ sigma parameter
            epsilon: LJ epsilon parameter
            
        Returns:
            External potential profile
        """
        r_cut: float = 1.5e-9
        dz: float = z[1] - z[0]
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
    
    def calculate_adsorption_isotherm(self, 
                                    pressure_range: Tuple[float, float] = (-4, 0), 
                                    n_points: int = 100) -> None:
        """
        Calculate the adsorption isotherm over a pressure range.
        
        Args:
            pressure_range: Log10 range of relative pressures (P/P0)
            n_points: Number of pressure points to calculate
        """
        self.initialize_arrays()
        self.U_ext = self._calculate_external_potential(
            self.z, self.rho_s, self.sigma_sf, self.eps_sf # type: ignore
        )
        
        self.pressures = np.logspace(*pressure_range, n_points)
        self.adsorption = []
        
        rho_bulk: float = (25e3 / 40e-3) * self.NA  # molecules/m³
        
        for p in self.pressures:
            mu: float = self.kB * self.T * np.log(max(p, 1e-12))  # protection against log(0)
            
            exponent = -(self.U_ext - mu) / (self.kB * self.T)
            exponent = np.clip(exponent, -100, 100)
            rho_f = rho_bulk * np.exp(exponent)
            
            N_ads: float = simpson(rho_f, self.z)  # type: ignore # molecules/m²
            self.adsorption.append(N_ads / self.NA * 1e3)  # mmol/m²
            
    def save_results(self, filename: str = "data_out.csv") -> None:
        """
        Save the calculated adsorption data to a file.
        
        Args:
            filename: Output file name
        """
        with open(filename, "w") as f:
            for p, N_ads in zip(self.pressures, self.adsorption): # type: ignore
                f.write(f"{p:.4e} {N_ads:.4e}\n")
                
    def plot_isotherm(self) -> None:
        """
        Plot the calculated adsorption isotherm.
        """
        if self.adsorption is None:
            raise ValueError("Calculate adsorption isotherm first using calculate_adsorption_isotherm()")
            
        plt.figure(figsize=(7, 5))
        plt.semilogx(self.pressures, self.adsorption, label='QSDFT-like model')
        plt.xlabel("Relative pressure P/P₀")
        plt.ylabel("Adsorption (mmol/m²)")
        plt.title(f"Adsorption isotherm of argon at {self.T} K (simplified QSDFT)")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()


# Example usage:
if __name__ == "__main__":
    calculator = QSDFT()
    calculator.calculate_adsorption_isotherm()
    calculator.save_results()
    calculator.plot_isotherm()