import numpy as np

def crossSection(X, Transfo, E_neutron):
    """
    Approximate cross section σ(E) [barn] for a given nuclide and reaction.
    Based on typical nuclear data from IAEA (simplified interpolation).
    """

    # ==================================  Check arguments  ==================================
    valid_nuclides = ['U235', 'U238', 'Pu239', 'Pu240', 'Xe135']
    valid_transfo = ['Fission', 'Capture']

    if X not in valid_nuclides:
        raise ValueError(f"No database entry for element {X}")
    if Transfo not in valid_transfo:
        raise ValueError(f"Unknown transformation type: {Transfo}")

    E = np.array(E_neutron, dtype=float)
    sigma = np.zeros_like(E)

    # --- Define reference values [barn] ---
    data = {
        'U235': {'Fission': (585, 1.2), 'Capture': (100, 0.1)},
        'U238': {'Fission': (2, 0.3), 'Capture': (300, 0.3)},
        'Pu239': {'Fission': (742, 1.8), 'Capture': (270, 0.1)},
        'Pu240': {'Fission': (0.1, 0.01), 'Capture': (300, 0.1)},
        'Xe135': {'Capture': (2.6e6, 10)}
    }

    sigma_th, sigma_fast = data[X][Transfo]
    
    # --- Interpolation law ---
    E_th = 0.025  # eV
    E_fast = 1e6  # eV

    sigma = np.where(
        E < E_fast,
        sigma_th * (E / E_th) ** (-0.5),  # decrease with E
        sigma_fast
    )

    # Limit physically
    sigma = np.clip(sigma, sigma_fast, sigma_th)

    return sigma.tolist()
