# crossSection.py

import numpy as np


def crossSection(X, Transfo, E_neutron):
    """
    Computes the cross section for a specific energy level
    ENDF database: https://www-nds.iaea.org/exfor/endf.htm

    ----------------
    :param X: string
        nuclide or nucleon (for a neutron), follows the atomic notation of the element,
        e.g.: for Uranium-235, X = 'U235'
    :param Transfo: string
        name of the considered transformation, should be one of:
        ['Fission', 'Capture']
    :param E_neutron: array-like
        energy(ies) of the incident neutron in [eV]
        expected range [1e-5; 2e7] [eV]
    :return: numpy array of same shape as E_neutron
        cross section in [barn]
    """

    E = np.array(E_neutron, dtype=float)

    # output
    sigma = np.zeros_like(E)

    # define energy groups
    E_th = 1.0     # eV
    E_fast = 1e6   # eV

    # thermal values
    sig = {
        # U235
        ("U235", "Fission", "th"): 580.0,
        ("U235", "Capture", "th"): 100.0,

        # U238
        ("U238", "Fission", "th"): 0.02,
        ("U238", "Capture", "th"): 2.7,

        # U236 (peu fissile, plutôt capture)
        ("U236", "Fission", "th"): 0.1,
        ("U236", "Capture", "th"): 5.0,

        # U237 (intermédiaire, capture modérée)
        ("U237", "Fission", "th"): 0.05,
        ("U237", "Capture", "th"): 3.0,

        # U239 (intermédiaire vers Np239)
        ("U239", "Fission", "th"): 1.0,
        ("U239", "Capture", "th"): 5.0,

        # Np239
        ("Np239", "Fission", "th"): 200.0,
        ("Np239", "Capture", "th"): 50.0,

        # Pu239
        ("Pu239", "Fission", "th"): 750.0,
        ("Pu239", "Capture", "th"): 270.0,

        # Pu240 (plutôt poison)
        ("Pu240", "Fission", "th"): 0.5,
        ("Pu240", "Capture", "th"): 290.0,

        # Pu241 (bien fissile)
        ("Pu241", "Fission", "th"): 900.0,
        ("Pu241", "Capture", "th"): 280.0,

        # Thorium + cycle U233
        ("Th232", "Fission", "th"): 0.01,
        ("Th232", "Capture", "th"): 7.0,

        ("Th233", "Fission", "th"): 150.0,
        ("Th233", "Capture", "th"): 30.0,

        ("Pa233", "Fission", "th"): 200.0,
        ("Pa233", "Capture", "th"): 50.0,

        ("U233", "Fission", "th"): 530.0,
        ("U233", "Capture", "th"): 45.0,

        # Xénon-135 (énorme poison thermique)
        ("Xe135", "Capture", "th"): 2.5e6,
    }

    # fast values (smaller fission XS, small captures)
    sig_fast = {
        # U235
        ("U235", "Fission"): 1.0,
        ("U235", "Capture"): 1.0,

        # U238
        ("U238", "Fission"): 0.3,
        ("U238", "Capture"): 0.5,

        # U236
        ("U236", "Fission"): 0.1,
        ("U236", "Capture"): 0.5,

        # U237
        ("U237", "Fission"): 0.1,
        ("U237", "Capture"): 0.5,

        # U239
        ("U239", "Fission"): 0.5,
        ("U239", "Capture"): 0.5,

        # Np239
        ("Np239", "Fission"): 1.0,
        ("Np239", "Capture"): 0.5,

        # Pu239
        ("Pu239", "Fission"): 1.5,
        ("Pu239", "Capture"): 1.0,

        # Pu240
        ("Pu240", "Fission"): 0.1,
        ("Pu240", "Capture"): 1.0,

        # Pu241
        ("Pu241", "Fission"): 1.5,
        ("Pu241", "Capture"): 1.0,

        # Thorium + cycle U233
        ("Th232", "Fission"): 0.01,
        ("Th232", "Capture"): 0.3,

        ("Th233", "Fission"): 0.8,
        ("Th233", "Capture"): 0.3,

        ("Pa233", "Fission"): 1.0,
        ("Pa233", "Capture"): 0.5,

        ("U233", "Fission"): 1.2,
        ("U233", "Capture"): 0.5,

        # Xe135 (poison rapide beaucoup moins fort)
        ("Xe135", "Capture"): 1e3,
    }


    # masks
    mask_th = E <= E_th
    mask_fast = E >= E_fast
    mask_mid = (~mask_th) & (~mask_fast)

    key_th = (X, Transfo, "th")
    key_fast = (X, Transfo)

    if key_th not in sig and key_fast not in sig_fast:
        print("\n WARNING : No cross section data for", X, "/", Transfo)
        return sigma

    # thermal region
    if key_th in sig:
        sigma[mask_th] = sig[key_th]

    # fast region
    if key_fast in sig_fast:
        sigma[mask_fast] = sig_fast[key_fast]

    # intermediate region: log-linear interpolation between thermal and fast values
    if np.any(mask_mid):
        # if one of the endpoints is missing, just use the one we have
        sigma_mid = np.zeros_like(E[mask_mid])

        if key_th in sig and key_fast in sig_fast:
            s_th = sig[key_th]
            s_fast = sig_fast[key_fast]
            logE = np.log(E[mask_mid])
            logE_th = np.log(E_th)
            logE_fast = np.log(E_fast)
            w = (logE - logE_th) / (logE_fast - logE_th)
            sigma_mid = s_th + w * (s_fast - s_th)
        elif key_th in sig:
            sigma_mid[:] = sig[key_th]
        elif key_fast in sig_fast:
            sigma_mid[:] = sig_fast[key_fast]

        sigma[mask_mid] = sigma_mid

    return sigma


# Simple test
if __name__ == "__main__":
    E_test = np.logspace(-5, 7, 5)
    cs = crossSection("U235", "Fission", E_test)
    print("E [eV] :", E_test)
    print("sigma  :", cs, "barn")
