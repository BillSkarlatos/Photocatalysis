import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Tuple, Optional
from enum import Enum
import scipy.constants as const

class CatalystState(Enum):
    """Possible states of a photocatalyst molecule"""
    GROUND = 0
    EXCITED_S1 = 1
    EXCITED_S2 = 2
    TRIPLET_T1 = 3
    DEACTIVATED = 4

@dataclass
class PhotocatalystMolecule:
    """Individual photocatalyst molecule with quantum state tracking"""
    state: CatalystState = CatalystState.GROUND
    excitation_time: Optional[float] = None
    coherence_phase: float = 0.0
    vibrational_level: int = 0
    
class QuantumPhotocatalysisSimulator:
    """
    Optimized Monte Carlo simulator for quantum-enhanced photocatalysis
    """
    
    def __init__(self, 
                 n_catalysts: int = 100,  # Reduced for faster computation
                 temperature: float = 298.15,
                 oxygen_concentration: float = 2.5e-4,
                 substrate_concentration: float = 0.01):
        
        # System parameters
        self.n_catalysts = n_catalysts
        self.temperature = temperature
        
        # Concentrations
        self.oxygen_conc = oxygen_concentration
        self.substrate_conc = substrate_concentration
        
        # Photophysical parameters (simplified but realistic)
        self.wavelength = 800e-9  # m
        self.photon_energy = const.h * const.c / self.wavelength
        
        # Absorption cross-sections
        self.sigma_1pa = 3e-20  # m^2
        self.sigma_2pa_classical = 1e-49  # m^4·s (~100 GM)
        self.entanglement_time = 100e-15  # s
        
        # State lifetimes
        self.tau_s1 = 5e-9  # s
        self.tau_s2 = 100e-15  # s
        self.tau_t1 = 1e-6  # s
        
        # Quantum yields
        self.phi_fluorescence = 0.1
        self.phi_isc = 0.6
        self.phi_reaction_s2 = 0.85
        self.phi_reaction_t1 = 0.4
        
        # Deactivation
        self.p_photobleaching = 1e-5
        self.p_side_reaction = 0.05
        
        # Initialize catalysts
        self.catalysts = [PhotocatalystMolecule() for _ in range(n_catalysts)]
        self.active_catalysts = set(range(n_catalysts))
        
        # Tracking
        self.products_formed = 0
        self.substrate_consumed = 0
        self.total_absorptions = 0
        self.two_photon_events = 0
        
    def calculate_entangled_enhancement(self) -> float:
        """Calculate enhancement factor for entangled two-photon absorption"""
        # Simplified enhancement model
        coherence_time = 1e-15
        enhancement = (self.entanglement_time / coherence_time) ** 0.5
        spectral_overlap = 0.8
        return min(enhancement * spectral_overlap * 100, 1000)  # Cap at 1000x
    
    def simulate_classical_illumination(self, 
                                      flux: float,
                                      duration: float) -> dict:
        """Optimized classical light simulation"""
        
        self._reset_system()
        
        # Use event-based simulation instead of time steps
        beam_area = 1e-6  # m^2
        intensity = flux / beam_area
        
        # Calculate expected number of events
        tpa_rate_per_molecule = self.sigma_2pa_classical * (intensity ** 2)
        total_tpa_rate = tpa_rate_per_molecule * len(self.active_catalysts)
        
        # Sample number of TPA events
        n_tpa_events = np.random.poisson(total_tpa_rate * duration)
        
        # Process events
        for _ in range(n_tpa_events):
            if len(self.active_catalysts) == 0:
                break
                
            cat_idx = np.random.choice(list(self.active_catalysts))
            cat = self.catalysts[cat_idx]
            
            if cat.state == CatalystState.GROUND:
                # Two-photon excitation to S2
                cat.state = CatalystState.EXCITED_S2
                self.two_photon_events += 1
                self.total_absorptions += 2
                
                # Immediate decay from S2 (very short lifetime)
                if np.random.random() < self.phi_reaction_s2:
                    # Direct reaction from S2
                    if self._attempt_reaction():
                        cat.state = CatalystState.GROUND
                else:
                    # IC to S1
                    cat.state = CatalystState.EXCITED_S1
                    
                    # S1 decay
                    rand = np.random.random()
                    if rand < self.phi_isc:
                        # ISC to triplet
                        cat.state = CatalystState.TRIPLET_T1
                        # Triplet reaction
                        if np.random.random() < self.phi_reaction_t1:
                            if self._attempt_reaction():
                                cat.state = CatalystState.GROUND
                    else:
                        cat.state = CatalystState.GROUND
                
                # Photobleaching check
                if np.random.random() < self.p_photobleaching:
                    cat.state = CatalystState.DEACTIVATED
                    self.active_catalysts.discard(cat_idx)
        
        return self._compile_results(duration, flux, "classical")
    
    def simulate_entangled_illumination(self,
                                       flux: float,
                                       duration: float) -> dict:
        """Optimized entangled photon simulation"""
        
        self._reset_system()
        
        # Entangled pairs
        pair_flux = flux / 2
        enhancement = self.calculate_entangled_enhancement()
        
        # Enhanced absorption probability per encounter
        absorption_prob = min(1.0, enhancement * self.sigma_2pa_classical * 1e20)
        
        # Number of pair arrivals
        n_pairs = np.random.poisson(pair_flux * duration)
        
        # Process each pair
        for _ in range(n_pairs):
            if len(self.active_catalysts) == 0:
                break
                
            cat_idx = np.random.choice(list(self.active_catalysts))
            cat = self.catalysts[cat_idx]
            
            if cat.state == CatalystState.GROUND and np.random.random() < absorption_prob:
                # Entangled TPA to S2
                cat.state = CatalystState.EXCITED_S2
                self.two_photon_events += 1
                self.total_absorptions += 2
                
                # Process excited state (same as classical)
                if np.random.random() < self.phi_reaction_s2:
                    if self._attempt_reaction():
                        cat.state = CatalystState.GROUND
                else:
                    cat.state = CatalystState.EXCITED_S1
                    
                    rand = np.random.random()
                    if rand < self.phi_isc:
                        cat.state = CatalystState.TRIPLET_T1
                        if np.random.random() < self.phi_reaction_t1:
                            if self._attempt_reaction():
                                cat.state = CatalystState.GROUND
                    else:
                        cat.state = CatalystState.GROUND
                
                # Photobleaching
                if np.random.random() < self.p_photobleaching * 2:
                    cat.state = CatalystState.DEACTIVATED
                    self.active_catalysts.discard(cat_idx)
        
        return self._compile_results(duration, flux, "entangled")
    
    def _attempt_reaction(self) -> bool:
        """Attempt a chemical reaction"""
        if self.substrate_conc > 0 and np.random.random() > self.p_side_reaction:
            self.products_formed += 1
            self.substrate_consumed += 1
            # Simple substrate depletion
            self.substrate_conc *= (1 - 0.001)
            return True
        return False
    
    def _reset_system(self):
        """Reset simulation state"""
        self.catalysts = [PhotocatalystMolecule() for _ in range(self.n_catalysts)]
        self.active_catalysts = set(range(self.n_catalysts))
        self.products_formed = 0
        self.substrate_consumed = 0
        self.total_absorptions = 0
        self.two_photon_events = 0
        self.substrate_conc = 0.01  # Reset substrate
    
    def _compile_results(self, duration: float, flux: float, light_type: str) -> dict:
        """Compile simulation results"""
        return {
            'light_type': light_type,
            'duration': duration,
            'flux': flux,
            'products_formed': self.products_formed,
            'substrate_consumed': self.substrate_consumed,
            'total_absorptions': self.total_absorptions,
            'two_photon_events': self.two_photon_events,
            'quantum_yield': self.products_formed / max(1, self.total_absorptions),
            'reaction_rate': self.products_formed / duration,
            'active_catalysts_remaining': len(self.active_catalysts),
            'photobleaching_fraction': 1 - len(self.active_catalysts) / self.n_catalysts
        }
    
    def run_flux_comparison(self, 
                           flux_range: np.ndarray,
                           duration: float = 1.0,
                           n_trials: int = 5) -> Tuple[dict, dict]:
        """Compare classical vs entangled performance"""
        
        classical_results = {
            'flux': flux_range,
            'rate_mean': [],
            'rate_std': [],
            'yield_mean': []
        }
        
        entangled_results = {
            'flux': flux_range,
            'rate_mean': [],
            'rate_std': [],
            'yield_mean': []
        }
        
        for flux in flux_range:
            print(f"Flux: {flux:.2e} photons/s", end=" ")
            
            classical_rates = []
            classical_yields = []
            entangled_rates = []
            entangled_yields = []
            
            for trial in range(n_trials):
                # Classical
                c_res = self.simulate_classical_illumination(flux, duration)
                classical_rates.append(c_res['reaction_rate'])
                classical_yields.append(c_res['quantum_yield'])
                
                # Entangled
                e_res = self.simulate_entangled_illumination(flux, duration)
                entangled_rates.append(e_res['reaction_rate'])
                entangled_yields.append(e_res['quantum_yield'])
            
            # Store statistics
            classical_results['rate_mean'].append(np.mean(classical_rates))
            classical_results['rate_std'].append(np.std(classical_rates))
            classical_results['yield_mean'].append(np.mean(classical_yields))
            
            entangled_results['rate_mean'].append(np.mean(entangled_rates))
            entangled_results['rate_std'].append(np.std(entangled_rates))
            entangled_results['yield_mean'].append(np.mean(entangled_yields))
            
            print(f"✓")
        
        return classical_results, entangled_results
    
    def plot_results(self, classical_results: dict, entangled_results: dict):
        """Create plots of results"""
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Reaction rate vs flux
        ax = axes[0, 0]
        ax.plot(classical_results['flux'], classical_results['rate_mean'],
                'o-', color='red', label='Classical', linewidth=2, markersize=8)
        ax.plot(entangled_results['flux'], entangled_results['rate_mean'],
                's-', color='blue', label='Entangled', linewidth=2, markersize=8)
        ax.set_xlabel('Photon Flux (photons/s)', fontsize=12)
        ax.set_ylabel('Reaction Rate (products/s)', fontsize=12)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_title('Reaction Rate Comparison', fontsize=14, fontweight='bold')
        
        # Quantum yield vs flux
        ax = axes[0, 1]
        ax.plot(classical_results['flux'], classical_results['yield_mean'],
                'o-', color='red', label='Classical', linewidth=2, markersize=8)
        ax.plot(entangled_results['flux'], entangled_results['yield_mean'],
                's-', color='blue', label='Entangled', linewidth=2, markersize=8)
        ax.set_xlabel('Photon Flux (photons/s)', fontsize=12)
        ax.set_ylabel('Quantum Yield', fontsize=12)
        ax.set_xscale('log')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_title('Quantum Yield Comparison', fontsize=14, fontweight='bold')
        
        # Enhancement factor
        ax = axes[1, 0]
        enhancement = np.array(entangled_results['rate_mean']) / (np.array(classical_results['rate_mean']) + 1e-10)
        ax.plot(classical_results['flux'], enhancement, 'g^-', linewidth=2, markersize=8)
        ax.set_xlabel('Photon Flux (photons/s)', fontsize=12)
        ax.set_ylabel('Enhancement Factor', fontsize=12)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        ax.set_title('Quantum Enhancement Factor', fontsize=14, fontweight='bold')
        ax.axhline(y=1, color='k', linestyle='--', alpha=0.5)
        ax.fill_between(classical_results['flux'], 1, enhancement, 
                        where=(enhancement > 1), alpha=0.2, color='green',
                        label='Quantum Advantage Region')
        ax.legend(fontsize=11)
        
        # Efficiency plot
        ax = axes[1, 1]
        classical_efficiency = np.array(classical_results['rate_mean']) / np.array(classical_results['flux'])
        entangled_efficiency = np.array(entangled_results['rate_mean']) / np.array(entangled_results['flux'])
        
        ax.plot(classical_results['flux'], classical_efficiency * 1e6, 'o-', 
                color='red', label='Classical', linewidth=2, markersize=8)
        ax.plot(entangled_results['flux'], entangled_efficiency * 1e6, 's-',
                color='blue', label='Entangled', linewidth=2, markersize=8)
        ax.set_xlabel('Photon Flux (photons/s)', fontsize=12)
        ax.set_ylabel('Efficiency (products per 10⁶ photons)', fontsize=12)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_title('Photon Utilization Efficiency', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        return fig

# Run the simulation
def main():
    print("="*60)
    print("QUANTUM-ENHANCED PHOTOCATALYSIS SIMULATION")
    print("="*60)
    print("\nInitializing simulator...")
    
    # Create simulator
    sim = QuantumPhotocatalysisSimulator(
        n_catalysts=100,  # Reduced for speed
        temperature=298.15,
        oxygen_concentration=2.5e-4,
        substrate_concentration=0.01
    )
    
    # Test flux range
    flux_range = np.logspace(2, 6, 15)  # 100 to 1,000,000 photons/s
    
    print(f"Running simulations across {len(flux_range)} flux values...")
    print(f"Duration: 1.0 s per simulation")
    print(f"Trials: 5 per condition\n")
    
    # Run comparison
    classical, entangled = sim.run_flux_comparison(
        flux_range=flux_range,
        duration=1.0,
        n_trials=5
    )
    
    print("\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)
    
    # Calculate enhancement
    c_rates = np.array(classical['rate_mean'])
    e_rates = np.array(entangled['rate_mean'])
    enhancement = e_rates / (c_rates + 1e-10)
    
    # Find key metrics
    max_enhancement_idx = np.argmax(enhancement)
    print(f"\nMaximum Enhancement: {enhancement[max_enhancement_idx]:.1f}x")
    print(f"At flux: {flux_range[max_enhancement_idx]:.2e} photons/s")
    
    # Low flux performance
    low_flux_idx = 2  # ~10^2.5 photons/s
    if low_flux_idx < len(flux_range):
        print(f"\nLow flux regime ({flux_range[low_flux_idx]:.2e} photons/s):")
        print(f"  Classical: {classical['rate_mean'][low_flux_idx]:.3e} products/s")
        print(f"  Entangled: {entangled['rate_mean'][low_flux_idx]:.3e} products/s")
        print(f"  Enhancement: {enhancement[low_flux_idx]:.1f}x")
    
    # High flux performance
    high_flux_idx = -1
    print(f"\nHigh flux regime ({flux_range[high_flux_idx]:.2e} photons/s):")
    print(f"  Classical: {classical['rate_mean'][high_flux_idx]:.3e} products/s")
    print(f"  Entangled: {entangled['rate_mean'][high_flux_idx]:.3e} products/s")
    print(f"  Enhancement: {enhancement[high_flux_idx]:.1f}x")
    
    # Create plots
    print("\nGenerating plots...")
    fig = sim.plot_results(classical, entangled)
    plt.show()
    
    print("\nSimulation complete!")
    return sim, classical, entangled, fig

# Execute
if __name__ == "__main__":
    sim, classical_results, entangled_results, figure = main()