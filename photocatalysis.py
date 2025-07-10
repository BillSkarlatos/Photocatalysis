import numpy as np
import matplotlib.pyplot as plt

# Physical constants
AVOGADRO = 6.022e23        # molecules/mol
PLANCK   = 6.626e-34       # J·s
C         = 3e8            # m/s

# Laser & SPDC parameters
AVG_POWER  = 0.1            # Average laser power [W]
WAVELENGTH = 800e-9         # Laser wavelength [m]
REP_RATE   = 80e6           # 80 MHz repetition rate
TAU_PULSE  = 100e-15        # 100 fs pulse duration
hnu = PLANCK * C / WAVELENGTH # Photon energy (J)

class PhotonSource:
    """
    Models a photon source: independent photons or entangled photon pairs.
    """
    def __init__(self, rate, entangled=False):
        self.rate = rate
        self.entangled = entangled

    def emit(self, dt):
        # Poisson-distributed photon (or pair) count in time dt
        return np.random.poisson(self.rate * dt)

class Payload:
    """
    A realistic photocatalytic medium for pulsed two-photon absorption:
      - Gaussian beam focus (beam waist w0, focal volume V)
      - Pulsed excitation (tau_pulse, rep_rate)
      - Correct one- and two-photon cross-section units
    """
    def __init__(self,
                 conc,                 # mol/L
                 w0=2e-4,              # beam waist [cm] (2 μm)
                 path_length=0.01,     # cm (100 μm)
                 sigma_single=1e-17,   # cm^2
                 sigma_two=1e-45,      # cm^4 s / photon (≈ 1000 GM)
                 wavelength=800e-9,    # m
                 tau_pulse=TAU_PULSE,  # s
                 rep_rate=REP_RATE     # Hz
                ):
        self.conc        = conc
        self.w0          = w0
        self.path_length = path_length
        # molecule density [#/cm^3]
        self.number_density = conc * 1e-3 * AVOGADRO
        # cross-sections
        self.sigma1 = sigma_single
        self.sigma2 = sigma_two
        # photon energy
        self.hnu     = PLANCK * C / wavelength
        # pulse parameters
        self.tau_pulse = tau_pulse
        self.rep_rate  = rep_rate
        # beam geometry
        self.beam_area    = np.pi * w0**2
        self.focal_volume = self.beam_area * path_length
        # molecules in focal volume
        self.N_mol = self.number_density * self.focal_volume

    def interact(self, photons, dt, entangled=False):
        # use fractional pulses in dt (even if <1)
        pulses = self.rep_rate * dt
        if photons <= 0 or pulses <= 0:
            return 0
        # photons per pulse
        ph_per_pulse = photons / pulses
        # pulse energy [J]
        E_pulse = ph_per_pulse * self.hnu
        # peak intensity [W/cm^2]
        I_peak = E_pulse / (self.beam_area * self.tau_pulse)
        # one-photon probability per molecule per pulse
        P1 = self.sigma1 * (I_peak / self.hnu) * self.tau_pulse
        # two-photon probability per molecule per pulse
        P2 = self.sigma2 * (I_peak / self.hnu)**2 * self.tau_pulse
        # clamp
        P1 = min(P1, 1.0)
        P2 = min(P2, 1.0)
        # mean events per pulse
        if entangled:
            mean_per_pulse = self.N_mol * P2
        else:
            mean_per_pulse = self.N_mol * P1**2
        # total mean in dt
        mean_total = pulses * mean_per_pulse
        return np.random.poisson(mean_total)


def run_sweep(concs, source_factory, payload_factory,
              T=1e-9, dt=1e-12, label=""):
    """
    Sweep over concentrations, returning yields array.
    concs: array of concentrations [M]
    source_factory(conc) -> PhotonSource
    payload_factory(conc) -> Payload
    T: total simulation time [s]
    dt: time step [s]
    label: prefix for prints
    """
    yields = np.zeros_like(concs)
    for i, c in enumerate(concs):
        src = source_factory(c)
        pld = payload_factory(c)
        total = 0
        t = 0.0
        while t < T:
            n = src.emit(dt)
            total += pld.interact(n, dt, entangled=src.entangled)
            t += dt
        yields[i] = total
        print(f"[{label}] c={c:.1e} M → yield={total}")
    return yields

if __name__ == '__main__':
    # concentration sweep: 50 pts from 1e-4 to 1 M
    concentrations = np.logspace(-4, 0, 50)

    def make_independent(conc):
        photon_rate = AVG_POWER / hnu
        return PhotonSource(rate=photon_rate, entangled=False)

    def make_entangled(conc):
        photon_rate = AVG_POWER / hnu
        return PhotonSource(rate=photon_rate, entangled=True)

    def make_payload(conc):
        return Payload(conc=conc)

    # run sweeps
    sweep_configs = {
        'Classical': make_independent,
        'Entangled': make_entangled
    }

    results = {}
    for label, src_fac in sweep_configs.items():
        print(f"Running {label} sweep...")
        results[label] = run_sweep(
            concentrations,
            source_factory=src_fac,
            payload_factory=make_payload,
            T=1e-9,
            dt=1e-12,
            label=label
        )

    # unpack
    y_ind = results['Classical']
    y_ent = results['Entangled']

    # plot
    plt.figure(figsize=(8, 5))
    plt.loglog(concentrations, y_ind+1, label='Classical', linewidth=2)
    plt.loglog(concentrations, y_ent+1, label='Entangled', linewidth=2)
    plt.xlabel('Concentration (M)')
    plt.ylabel('Yield (events per 1 ns + 1)')
    plt.title('Photocatalytic Yield vs. Concentration (1 ns)')
    plt.legend()
    plt.grid(True, which='both', ls='--', lw=0.5)
    plt.tight_layout()
    plt.show()