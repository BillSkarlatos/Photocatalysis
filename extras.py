import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.constants import Avogadro, h, c

class PhotocatalysisModel:
    """
    Kinetic model for exciton dynamics in photocatalysis,
    including first-order decay, bimolecular annihilation,
    and optional surface adsorption/desorption.
    """
    def __init__(self,
                 k_rec,          # Radiative recombination rate (s^-1)
                 k_nr,           # Non-radiative decay rate (s^-1)
                 k_chem,         # Chemical (productive) reaction rate constant (s^-1 for homogeneous, M^-1 s^-1 for surface)
                 k_xx,           # Bimolecular exciton-exciton annihilation rate constant (cm^3/s)
                 k_ads=None,     # Adsorption rate constant (M^-1 s^-1)
                 k_des=None,     # Desorption rate constant (s^-1)
                 S_bulk=None,    # Bulk substrate concentration (M)
                 S_tot=None):    # Total surface site concentration (M)
        self.k_rec = k_rec
        self.k_nr  = k_nr
        self.k_chem = k_chem
        self.k_xx  = k_xx
        # surface parameters (None => homogeneous reaction)
        self.k_ads   = k_ads
        self.k_des   = k_des
        self.S_bulk  = S_bulk
        self.S_tot   = S_tot

    def kinetics(self, y, t, G):
        """
        ODE system:
          y[0] = n_exc  (exciton concentration)
          y[1] = S_ads  (adsorbed substrate)
          y[2] = P      (product)
        G = exciton generation rate (cm^-3 s^-1)
        """
        n, S_ads, P = y
        # Reaction consumption term
        if self.k_ads is not None and self.k_des is not None:
            reaction = self.k_chem * n * S_ads
            # surface adsorption/desorption kinetics
            dSads_dt = (self.k_ads * self.S_bulk * (self.S_tot - S_ads)
                        - self.k_des * S_ads
                        - reaction)
        else:
            reaction = self.k_chem * n
            dSads_dt = 0
        # exciton dynamics: generation minus losses
        dn_dt = (G
                 - (self.k_rec + self.k_nr) * n
                 - self.k_xx * n**2
                 - reaction)
        # product formation
        dP_dt = reaction
        return [dn_dt, dSads_dt, dP_dt]

    def simulate(self, G, t):
        """
        Numerically integrate ODEs over time array t for a given generation rate G.
        Returns time traces for n_exc and product P.
        """
        # initial conditions: no excitons, no ads, no product
        y0 = [0.0,
              0.0 if self.k_ads is not None else 0.0,
              0.0]
        sol = odeint(self.kinetics, y0, t, args=(G,))
        n_trace = sol[:, 0]
        P_trace = sol[:, 2]
        return n_trace, P_trace

    def steady_state(self, G):
        """
        Analytical steady-state solution for exciton density and quantum yield.
        Solves:     G = (k_rec + k_nr + k_chem) n + k_xx n^2
        Then R_prod = k_chem * n,  phi = R_prod / G
        """
        k_tot = self.k_rec + self.k_nr + self.k_chem
        if self.k_xx == 0:
            n_ss = G / k_tot
        else:
            # quadratic: k_xx n^2 + k_tot n - G = 0
            disc = np.sqrt(k_tot**2 + 4 * self.k_xx * G)
            n_ss = (-k_tot + disc) / (2 * self.k_xx)
        R_prod = self.k_chem * n_ss
        phi = R_prod / G
        return n_ss, R_prod, phi

if __name__ == "__main__":
    # Example parameters (tune to your material)
    k_rec = 1e8      # s^-1
    k_nr  = 1e7      # s^-1
    k_chem = 1e6     # s^-1 (homogeneous) or M^-1 s^-1 (surface)
    k_xx  = 1e-10    # cm^3/s
    # Optional surface kinetics
    # k_ads = 1e5    # M^-1 s^-1
    # k_des = 1e3    # s^-1
    # S_bulk = 1e-3  # M
    # S_tot  = 1e-6  # M

    model = PhotocatalysisModel(k_rec, k_nr, k_chem, k_xx)

    # Sweep generation rates
    G_vals = np.logspace(15, 25, 200)  # cm^-3 s^-1
    phi_vals = [model.steady_state(G)[2] for G in G_vals]

    plt.figure(figsize=(6, 4))
    plt.loglog(G_vals, phi_vals, linewidth=2)
    plt.xlabel('Exciton generation rate $G$ (cm$^{-3}$ s$^{-1}$)')
    plt.ylabel('Quantum yield $\phi$')
    plt.title('Quantum Yield vs. Generation Rate')
    plt.grid(True, which='both', ls='--', lw=0.5)
    plt.tight_layout()
    plt.show()
