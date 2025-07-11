import numpy as np
import matplotlib.pyplot as plt
import random

# ------------------------------
# Simulation Functions
# ------------------------------
def simulate_classical(N, photon_rate, time_duration, tau,
                       p_deactivation, p_side_reaction,
                       substrate_amount=None, seed=None):
    """
    Monte Carlo simulation of classical two-photon photocatalysis.

    Parameters:
      N               : number of catalyst molecules
      photon_rate     : incident photons per second (Φ)
      time_duration   : total simulation time (s)
      tau             : excited-state lifetime (s)
      p_deactivation  : probability of catalyst deactivation per absorption
      p_side_reaction : probability that PC** does not yield product
      substrate_amount: initial substrate count (None for infinite)
      seed            : random seed for reproducibility

    Returns:
      products        : total product count
    """
    # Seed randomness
    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    # State arrays
    active = np.ones(N, dtype=bool)
    excited = np.zeros(N, dtype=bool)
    last_excitation = np.zeros(N)
    products = 0
    substrate = substrate_amount if substrate_amount is not None else np.inf

    # Determine number of photons
    expected_photons = photon_rate * time_duration
    M = np.random.poisson(expected_photons)
    # Sample arrival times
    times = np.random.rand(M) * time_duration
    times.sort()

    for t in times:
        if substrate <= 0:
            break
        # pick a random catalyst index
        i = np.random.randint(N)
        if not active[i]:
            continue
        if not excited[i]:
            # ground to PC*
            excited[i] = True
            last_excitation[i] = t
        else:
            dt = t - last_excitation[i]
            if dt <= tau:
                # two-photon event: PC**
                if random.random() >= p_side_reaction:
                    products += 1
                    substrate -= 1
                excited[i] = False
            else:
                # PC* decayed; treat this photon as new excitation
                last_excitation[i] = t
                excited[i] = True
        # catalyst deactivation
        if active[i] and random.random() < p_deactivation:
            active[i] = False
            excited[i] = False

    return products


def simulate_entangled(N, pair_rate, time_duration,
                       p_deactivation, p_side_reaction,
                       substrate_amount=None, seed=None):
    """
    Monte Carlo simulation with entangled photon pairs.

    Parameters:
      N               : number of catalyst molecules
      pair_rate       : entangled pairs per second (Φ/2)
      time_duration   : total simulation time (s)
      p_deactivation  : deactivation prob per absorption
      p_side_reaction : prob PC** fails to yield product
      substrate_amount: initial substrate (None for infinite)
      seed            : random seed

    Returns:
      products        : total product count
    """
    if seed is not None:
        np.random.seed(seed)
        random.seed(seed)

    active = np.ones(N, dtype=bool)
    products = 0
    substrate = substrate_amount if substrate_amount is not None else np.inf

    expected_pairs = pair_rate * time_duration
    K = np.random.poisson(expected_pairs)
    times = np.random.rand(K) * time_duration
    times.sort()

    for t in times:
        if substrate <= 0:
            break
        i = np.random.randint(N)
        if not active[i]:
            continue
        # direct two-photon absorption => reaction
        if random.random() >= p_side_reaction:
            products += 1
            substrate -= 1
        # deactivation
        if random.random() < p_deactivation:
            active[i] = False

    return products

# ------------------------------
# Main Routine and Plots
# ------------------------------
if __name__ == '__main__':
    # Common parameters
    N = 1000                  # catalyst count
    tau = 1e-6                # excited-state lifetime (s)
    p_deactivation = 1e-4     # per absorption
    p_side = 0.0              # ideal product yield
    T = 1.0                   # simulation time (s)

    # (a) Reaction rate vs. photon flux
    fluxes = np.logspace(2, 5, 10)  # photons/s from 1e2 to 1e5
    rates_classical = []
    rates_entangled = []

    for phi in fluxes:
        rate_c = simulate_classical(N, phi, T, tau, p_deactivation, p_side)
        rate_q = simulate_entangled(N, phi/2, T, p_deactivation, p_side)
        rates_classical.append(rate_c / T)
        rates_entangled.append(rate_q / T)

    plt.figure(figsize=(8, 5))
    plt.loglog(fluxes, rates_classical, 'r--o', label='Classical')
    plt.loglog(fluxes, rates_entangled, 'b-o', label='Entangled')
    plt.xlabel('Photon Flux (photons/s)')
    plt.ylabel('Reaction Rate (products/s)')
    plt.title('Reaction Rate vs. Photon Flux')
    plt.legend()
    plt.grid(True, which='both', ls='--', lw=0.5)
    plt.tight_layout()
    plt.show()

    # (b) Yield vs. substrate concentration for fixed photon dose
    phi = 5e3                       # photons/s
    dose_time = 1.0                 # s
    pairs_per_s = phi / 2           # entangled pairs rate
    S_list = [100, 1000, 5000, 10000]
    yield_classical = []
    yield_entangled = []

    for S0 in S_list:
        y_c = simulate_classical(N, phi, dose_time, tau, p_deactivation, p_side, substrate_amount=S0)
        y_q = simulate_entangled(N, pairs_per_s, dose_time, p_deactivation, p_side, substrate_amount=S0)
        yield_classical.append(100 * y_c / S0)
        yield_entangled.append(100 * y_q / S0)

    plt.figure(figsize=(8, 5))
    plt.plot(S_list, yield_classical, 'r--o', label='Classical')
    plt.plot(S_list, yield_entangled, 'b-o', label='Entangled')
    plt.xscale('log')
    plt.xlabel('Initial Substrate Count')
    plt.ylabel('Conversion Yield (%)')
    plt.title('Substrate Conversion vs. Initial Amount')
    plt.legend()
    plt.grid(True, which='both', ls='--', lw=0.5)
    plt.tight_layout()
    plt.show()
