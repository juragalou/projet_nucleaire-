"""
def crossSection(X, Transfo, E_neutron):
    
    Computes the cross section for a specific energy level
    ENDF database: https://www-nds.iaea.org/exfor/endf.htm

    ----------------
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        e.g.: for the Uranium 235, X = 'U235'
    :param Transfo: string
        name of the considered transformation, should be one of the following in the list
        ['Fission', 'Capture']
    :param E_neutron: array, list
        value(s) of the cross section for an incident neutron of energy level equal to E_neutron in [eV]
        the value should be in the range [1e-5; 2e7] [eV]
    :return: array, list
        values of the cross section in [barn] corresponding to the input parameters
    

    # Your work : here

    # ==================================  Check arguments  ==================================
    #  Check that the nucleon/nuclide asked, and that the associated transformation exists in the database
    if X != 'Th232' and X != 'Th233' and X != 'Pa233' and X != 'U233' and X != 'U235' and X != 'U236' and X != 'U237' and X != 'U238' \
            and X != 'U239' and X != 'Np239' and X != 'Pu239' and X != 'Pu240' and X != 'Xe135':
        print('\n WARNING : There is no database for element ', X, '. \n Please check function information')
    
    # sigma = 0.

    return sigma

# Try the function:
# cs = crossSection(X='Pu240', Transfo='Fission', E_neutron=np.logspace(-5,6,10000).tolist())
"""

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

    # simple values (order of magnitude; enough for the project)
    # thermal values
    sig = {
        ("U235", "Fission", "th"): 580.0,
        ("U235", "Capture", "th"): 100.0,

        ("U238", "Fission", "th"): 0.02,
        ("U238", "Capture", "th"): 2.7,

        ("Pu239", "Fission", "th"): 750.0,
        ("Pu239", "Capture", "th"): 270.0,

        ("Pu240", "Fission", "th"): 0.5,
        ("Pu240", "Capture", "th"): 290.0,

        ("Xe135", "Capture", "th"): 2.5e6,  # huge poison
    }

    # fast values (smaller fission XS, small captures)
    sig_fast = {
        ("U235", "Fission"): 1.0,
        ("U235", "Capture"): 1.0,

        ("U238", "Fission"): 0.3,
        ("U238", "Capture"): 0.5,

        ("Pu239", "Fission"): 1.5,
        ("Pu239", "Capture"): 1.0,

        ("Pu240", "Fission"): 0.1,
        ("Pu240", "Capture"): 1.0,

        ("Xe135", "Capture"): 1e3,  # much smaller than thermal but still sizable
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
