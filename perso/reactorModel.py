import numpy as np
import matplotlib.pyplot as plt

import molarMass as mM
import halfLife as hL
import crossSection as cS

# ------------------- CONSTANTES PHYSIQUES -------------------

NA = 6.022e23           # [1/mol]
V_CORE = 10.0           # [m^3] volume du coeur
Q_FISSION = 200e6 * 1.602e-19   # [J/fission] ~200 MeV
NU = 2.4                # neutrons par fission (moyenne PWR)
BETA = 0.0065           # fraction de neutrons retardés

# Groupes d'énergie
E_TH = 0.025      # [eV] thermique
E_FAST = 1e6      # [eV] rapide
EV_TO_J = 1.602e-19
M_NEUTRON = 1.6749e-27  # [kg]

# Ralentissement fast -> thermal
T_SLOW = 5e-4
LAMBDA_SLOW = np.log(2.0) / T_SLOW

# Pertes (barres + fuites) [1/s] (à ajuster)
SIGMA_CTRL_FAST = 5.0
SIGMA_CTRL_TH   = 5.0


def decay_constant(species, transfo):
    T = hL.halfLife(species, transfo)
    if T <= 0:
        return 0.0
    return np.log(2.0) / T


# Décroissances
LAMBDA_FP    = decay_constant("FP",    "BetaMinus")
LAMBDA_XE    = decay_constant("Xe135", "BetaMinus")
LAMBDA_U239  = decay_constant("U239",  "BetaMinus")
LAMBDA_NP239 = decay_constant("Np239", "BetaMinus")
LAMBDA_TH233 = decay_constant("Th233", "BetaMinus")
LAMBDA_PA233 = decay_constant("Pa233", "BetaMinus")
LAMBDA_PU241 = decay_constant("Pu241", "BetaMinus")  # Pu-241 -> Am-241


def reactorModel(fuelCompo, FPCompo, t_final, n_th_init, n_fa_init, mTot):
    """
    Modèle cinétique avec toutes les filières possibles du fuel (U, Pu, Th).
    """

    # --------- 1. COMPOSITION INITIALE DU FUEL ---------

    # masses [kg]
    m_U235  = mTot * fuelCompo.U235  / 100.0
    m_U238  = mTot * fuelCompo.U238  / 100.0
    m_Pu239 = mTot * fuelCompo.Pu239 / 100.0
    m_Th232 = mTot * fuelCompo.Th232 / 100.0

    # nombres de noyaux
    N_U235  = m_U235  / mM.molarMass("U235")  * NA
    N_U238  = m_U238  / mM.molarMass("U238")  * NA
    N_Pu239 = m_Pu239 / mM.molarMass("Pu239") * NA
    N_Th232 = m_Th232 / mM.molarMass("Th232") * NA

    # noyaux créés pendant l'irradiation
    N_U236  = 0.0
    N_U237  = 0.0
    N_U239  = 0.0
    N_Np239 = 0.0
    N_Pu240 = 0.0
    N_Pu241 = 0.0
    N_Th233 = 0.0
    N_Pa233 = 0.0
    N_U233  = 0.0

    # produits de fission
    N_FP = 0.0
    N_Xe = 0.0

    # neutrons
    n_fast = n_fa_init
    n_th   = n_th_init

    # fraction de PF allant dans Xe135
    y_XE = FPCompo.Xe135 / 100.0

    # --------- 2. VITESSES + SECTIONS EFFICACES ---------

    v_th   = np.sqrt(2.0 * E_TH   * EV_TO_J / M_NEUTRON)
    v_fast = np.sqrt(2.0 * E_FAST * EV_TO_J / M_NEUTRON)

    def xs_th_fast(nuclide, transfo):
        vals = cS.crossSection(nuclide, transfo, [E_TH, E_FAST])
        return vals[0] * 1e-28, vals[1] * 1e-28

    # Uranium
    sigma_U235_fis_th, sigma_U235_fis_fast = xs_th_fast("U235", "Fission")
    sigma_U235_cap_th, sigma_U235_cap_fast = xs_th_fast("U235", "Capture")
    sigma_U236_fis_th, sigma_U236_fis_fast = xs_th_fast("U236", "Fission")
    sigma_U236_cap_th, sigma_U236_cap_fast = xs_th_fast("U236", "Capture")
    sigma_U237_fis_th, sigma_U237_fis_fast = xs_th_fast("U237", "Fission")
    sigma_U237_cap_th, sigma_U237_cap_fast = xs_th_fast("U237", "Capture")
    sigma_U238_fis_th, sigma_U238_fis_fast = xs_th_fast("U238", "Fission")
    sigma_U238_cap_th, sigma_U238_cap_fast = xs_th_fast("U238", "Capture")
    sigma_U239_fis_th, sigma_U239_fis_fast = xs_th_fast("U239", "Fission")
    sigma_U239_cap_th, sigma_U239_cap_fast = xs_th_fast("U239", "Capture")

    # Np / Pu
    sigma_Np239_fis_th, sigma_Np239_fis_fast = xs_th_fast("Np239", "Fission")
    sigma_Np239_cap_th, sigma_Np239_cap_fast = xs_th_fast("Np239", "Capture")
    sigma_Pu239_fis_th, sigma_Pu239_fis_fast = xs_th_fast("Pu239", "Fission")
    sigma_Pu239_cap_th, sigma_Pu239_cap_fast = xs_th_fast("Pu239", "Capture")
    sigma_Pu240_fis_th, sigma_Pu240_fis_fast = xs_th_fast("Pu240", "Fission")
    sigma_Pu240_cap_th, sigma_Pu240_cap_fast = xs_th_fast("Pu240", "Capture")
    sigma_Pu241_fis_th, sigma_Pu241_fis_fast = xs_th_fast("Pu241", "Fission")
    sigma_Pu241_cap_th, sigma_Pu241_cap_fast = xs_th_fast("Pu241", "Capture")

    # Filière thorium
    sigma_Th232_cap_th, sigma_Th232_cap_fast = xs_th_fast("Th232", "Capture")
    sigma_Th233_fis_th, sigma_Th233_fis_fast = xs_th_fast("Th233", "Fission")
    sigma_Th233_cap_th, sigma_Th233_cap_fast = xs_th_fast("Th233", "Capture")
    sigma_Pa233_fis_th, sigma_Pa233_fis_fast = xs_th_fast("Pa233", "Fission")
    sigma_Pa233_cap_th, sigma_Pa233_cap_fast = xs_th_fast("Pa233", "Capture")
    sigma_U233_fis_th,  sigma_U233_fis_fast  = xs_th_fast("U233",  "Fission")
    sigma_U233_cap_th,  sigma_U233_cap_fast  = xs_th_fast("U233",  "Capture")

    # Xe
    sigma_Xe_cap_th, sigma_Xe_cap_fast = xs_th_fast("Xe135", "Capture")
    kappa_xe = sigma_Xe_cap_th * v_th / V_CORE

    # --------- 3. TEMPS + TABLEAUX ---------

    dt      = 1e-4
    n_steps = int(t_final / dt)

    t_arr      = np.zeros(n_steps)
    P_arr      = np.zeros(n_steps)
    n_fast_arr = np.zeros(n_steps)
    n_th_arr   = np.zeros(n_steps)

    N_U235_arr  = np.zeros(n_steps)
    N_U236_arr  = np.zeros(n_steps)
    N_U237_arr  = np.zeros(n_steps)
    N_U238_arr  = np.zeros(n_steps)
    N_U239_arr  = np.zeros(n_steps)
    N_Np239_arr = np.zeros(n_steps)
    N_Pu239_arr = np.zeros(n_steps)
    N_Pu240_arr = np.zeros(n_steps)
    N_Pu241_arr = np.zeros(n_steps)
    N_Th232_arr = np.zeros(n_steps)
    N_Th233_arr = np.zeros(n_steps)
    N_Pa233_arr = np.zeros(n_steps)
    N_U233_arr  = np.zeros(n_steps)
    N_FP_arr    = np.zeros(n_steps)
    N_Xe_arr    = np.zeros(n_steps)

    # --------- 4. BOUCLE TEMPORELLE ---------

    for k in range(n_steps):
        t = k * dt
        t_arr[k] = t

        flux_th   = n_th   * v_th   / V_CORE
        flux_fast = n_fast * v_fast / V_CORE

        # --- Uranium ---

        R_U235_fis_th   = sigma_U235_fis_th   * flux_th   * N_U235
        R_U235_fis_fast = sigma_U235_fis_fast * flux_fast * N_U235
        R_U235_cap_th   = sigma_U235_cap_th   * flux_th   * N_U235
        R_U235_cap_fast = sigma_U235_cap_fast * flux_fast * N_U235
        R_U235_fis = R_U235_fis_th + R_U235_fis_fast
        R_U235_cap = R_U235_cap_th + R_U235_cap_fast

        R_U236_fis_th   = sigma_U236_fis_th   * flux_th   * N_U236
        R_U236_fis_fast = sigma_U236_fis_fast * flux_fast * N_U236
        R_U236_cap_th   = sigma_U236_cap_th   * flux_th   * N_U236
        R_U236_cap_fast = sigma_U236_cap_fast * flux_fast * N_U236
        R_U236_fis = R_U236_fis_th + R_U236_fis_fast
        R_U236_cap = R_U236_cap_th + R_U236_cap_fast

        R_U237_fis_th   = sigma_U237_fis_th   * flux_th   * N_U237
        R_U237_fis_fast = sigma_U237_fis_fast * flux_fast * N_U237
        R_U237_cap_th   = sigma_U237_cap_th   * flux_th   * N_U237
        R_U237_cap_fast = sigma_U237_cap_fast * flux_fast * N_U237
        R_U237_fis = R_U237_fis_th + R_U237_fis_fast
        R_U237_cap = R_U237_cap_th + R_U237_cap_fast

        R_U238_fis_th   = sigma_U238_fis_th   * flux_th   * N_U238
        R_U238_fis_fast = sigma_U238_fis_fast * flux_fast * N_U238
        R_U238_cap_th   = sigma_U238_cap_th   * flux_th   * N_U238
        R_U238_cap_fast = sigma_U238_cap_fast * flux_fast * N_U238
        R_U238_fis = R_U238_fis_th + R_U238_fis_fast
        R_U238_cap = R_U238_cap_th + R_U238_cap_fast

        R_U239_fis_th   = sigma_U239_fis_th   * flux_th   * N_U239
        R_U239_fis_fast = sigma_U239_fis_fast * flux_fast * N_U239
        R_U239_cap_th   = sigma_U239_cap_th   * flux_th   * N_U239
        R_U239_cap_fast = sigma_U239_cap_fast * flux_fast * N_U239
        R_U239_fis = R_U239_fis_th + R_U239_fis_fast
        R_U239_cap = R_U239_cap_th + R_U239_cap_fast
        R_U239_beta = LAMBDA_U239 * N_U239

        # --- Np239 ---

        R_Np239_fis_th   = sigma_Np239_fis_th   * flux_th   * N_Np239
        R_Np239_fis_fast = sigma_Np239_fis_fast * flux_fast * N_Np239
        R_Np239_cap_th   = sigma_Np239_cap_th   * flux_th   * N_Np239
        R_Np239_cap_fast = sigma_Np239_cap_fast * flux_fast * N_Np239
        R_Np239_fis = R_Np239_fis_th + R_Np239_fis_fast
        R_Np239_cap = R_Np239_cap_th + R_Np239_cap_fast
        R_Np239_beta = LAMBDA_NP239 * N_Np239

        # --- Pu239 ---

        R_Pu239_fis_th   = sigma_Pu239_fis_th   * flux_th   * N_Pu239
        R_Pu239_fis_fast = sigma_Pu239_fis_fast * flux_fast * N_Pu239
        R_Pu239_cap_th   = sigma_Pu239_cap_th   * flux_th   * N_Pu239
        R_Pu239_cap_fast = sigma_Pu239_cap_fast * flux_fast * N_Pu239
        R_Pu239_fis = R_Pu239_fis_th + R_Pu239_fis_fast
        R_Pu239_cap = R_Pu239_cap_th + R_Pu239_cap_fast

        # --- Pu240 ---

        R_Pu240_fis_th   = sigma_Pu240_fis_th   * flux_th   * N_Pu240
        R_Pu240_fis_fast = sigma_Pu240_fis_fast * flux_fast * N_Pu240
        R_Pu240_cap_th   = sigma_Pu240_cap_th   * flux_th   * N_Pu240
        R_Pu240_cap_fast = sigma_Pu240_cap_fast * flux_fast * N_Pu240
        R_Pu240_fis = R_Pu240_fis_th + R_Pu240_fis_fast
        R_Pu240_cap = R_Pu240_cap_th + R_Pu240_cap_fast

        # --- Pu241 ---

        R_Pu241_fis_th   = sigma_Pu241_fis_th   * flux_th   * N_Pu241
        R_Pu241_fis_fast = sigma_Pu241_fis_fast * flux_fast * N_Pu241
        R_Pu241_cap_th   = sigma_Pu241_cap_th   * flux_th   * N_Pu241
        R_Pu241_cap_fast = sigma_Pu241_cap_fast * flux_fast * N_Pu241
        R_Pu241_fis = R_Pu241_fis_th + R_Pu241_fis_fast
        R_Pu241_cap = R_Pu241_cap_th + R_Pu241_cap_fast
        R_Pu241_beta = LAMBDA_PU241 * N_Pu241   # Pu-241 -> Am-241 (non suivi, considéré comme perte)

        # --- Filière thorium ---

        R_Th232_cap_th   = sigma_Th232_cap_th   * flux_th   * N_Th232
        R_Th232_cap_fast = sigma_Th232_cap_fast * flux_fast * N_Th232
        R_Th232_cap = R_Th232_cap_th + R_Th232_cap_fast

        R_Th233_fis_th   = sigma_Th233_fis_th   * flux_th   * N_Th233
        R_Th233_fis_fast = sigma_Th233_fis_fast * flux_fast * N_Th233
        R_Th233_cap_th   = sigma_Th233_cap_th   * flux_th   * N_Th233
        R_Th233_cap_fast = sigma_Th233_cap_fast * flux_fast * N_Th233
        R_Th233_fis = R_Th233_fis_th + R_Th233_fis_fast
        R_Th233_cap = R_Th233_cap_th + R_Th233_cap_fast
        R_Th233_beta = LAMBDA_TH233 * N_Th233

        R_Pa233_fis_th   = sigma_Pa233_fis_th   * flux_th   * N_Pa233
        R_Pa233_fis_fast = sigma_Pa233_fis_fast * flux_fast * N_Pa233
        R_Pa233_cap_th   = sigma_Pa233_cap_th   * flux_th   * N_Pa233
        R_Pa233_cap_fast = sigma_Pa233_cap_fast * flux_fast * N_Pa233
        R_Pa233_fis = R_Pa233_fis_th + R_Pa233_fis_fast
        R_Pa233_cap = R_Pa233_cap_th + R_Pa233_cap_fast
        R_Pa233_beta = LAMBDA_PA233 * N_Pa233

        R_U233_fis_th   = sigma_U233_fis_th   * flux_th   * N_U233
        R_U233_fis_fast = sigma_U233_fis_fast * flux_fast * N_U233
        R_U233_cap_th   = sigma_U233_cap_th   * flux_th   * N_U233
        R_U233_cap_fast = sigma_U233_cap_fast * flux_fast * N_U233
        R_U233_fis = R_U233_fis_th + R_U233_fis_fast
        R_U233_cap = R_U233_cap_th + R_U233_cap_fast

        # --- Xénon ---

        R_Xe_cap_th = sigma_Xe_cap_th * flux_th * N_Xe

        # Fissions totales
        F_tot = (
            R_U235_fis + R_U236_fis + R_U237_fis + R_U238_fis + R_U239_fis +
            R_Np239_fis +
            R_Pu239_fis + R_Pu240_fis + R_Pu241_fis +
            R_Th233_fis + R_Pa233_fis + R_U233_fis
        )

        # Fissions / captures par groupe d'énergie
        R_fis_fast = (
            R_U235_fis_fast + R_U236_fis_fast + R_U237_fis_fast +
            R_U238_fis_fast + R_U239_fis_fast + R_Np239_fis_fast +
            R_Pu239_fis_fast + R_Pu240_fis_fast + R_Pu241_fis_fast +
            R_Th233_fis_fast + R_Pa233_fis_fast + R_U233_fis_fast
        )
        R_fis_th = (
            R_U235_fis_th + R_U236_fis_th + R_U237_fis_th +
            R_U238_fis_th + R_U239_fis_th + R_Np239_fis_th +
            R_Pu239_fis_th + R_Pu240_fis_th + R_Pu241_fis_th +
            R_Th233_fis_th + R_Pa233_fis_th + R_U233_fis_th
        )

        R_cap_fast = (
            R_U235_cap_fast + R_U236_cap_fast + R_U237_cap_fast +
            R_U238_cap_fast + R_U239_cap_fast +
            R_Np239_cap_fast + R_Pu239_cap_fast + R_Pu240_cap_fast + R_Pu241_cap_fast +
            R_Th232_cap_fast + R_Th233_cap_fast +
            R_Pa233_cap_fast + R_U233_cap_fast
        )
        R_cap_th = (
            R_U235_cap_th + R_U236_cap_th + R_U237_cap_th +
            R_U238_cap_th + R_U239_cap_th +
            R_Np239_cap_th + R_Pu239_cap_th + R_Pu240_cap_th + R_Pu241_cap_th +
            R_Th232_cap_th + R_Th233_cap_th +
            R_Pa233_cap_th + R_U233_cap_th +
            R_Xe_cap_th
        )

        # Source retardée
        S_d = BETA * LAMBDA_FP * N_FP

        # --- Neutrons ---

        dn_fast_dt = (
            NU * F_tot
            - R_fis_fast
            - R_cap_fast
            - SIGMA_CTRL_FAST * n_fast
            - LAMBDA_SLOW * n_fast
            + S_d
        )

        dn_th_dt = (
            LAMBDA_SLOW * n_fast
            - R_fis_th
            - R_cap_th
            - SIGMA_CTRL_TH * n_th
        )

        # --- Actinides U / Np / Pu ---

        dN_U235_dt = -R_U235_fis - R_U235_cap
        dN_U236_dt = -R_U236_fis - R_U236_cap + R_U235_cap
        dN_U237_dt = -R_U237_fis - R_U237_cap + R_U236_cap
        dN_U238_dt = -R_U238_fis - R_U238_cap + R_U237_cap
        dN_U239_dt = -R_U239_fis - R_U239_cap - R_U239_beta + R_U238_cap

        dN_Np239_dt = -R_Np239_fis - R_Np239_cap - R_Np239_beta + R_U239_beta

        dN_Pu239_dt = -R_Pu239_fis - R_Pu239_cap + R_Np239_beta
        dN_Pu240_dt = -R_Pu240_fis - R_Pu240_cap + R_Pu239_cap
        dN_Pu241_dt = -R_Pu241_fis - R_Pu241_cap - R_Pu241_beta + R_Pu240_cap

        # --- Filière thorium ---

        dN_Th232_dt = -R_Th232_cap
        dN_Th233_dt = -R_Th233_fis - R_Th233_cap - R_Th233_beta + R_Th232_cap
        dN_Pa233_dt = -R_Pa233_fis - R_Pa233_cap - R_Pa233_beta + R_Th233_beta
        dN_U233_dt  = -R_U233_fis  - R_U233_cap  + R_Pa233_beta

        # --- FP / Xe ---

        dN_FP_dt = (1.0 - y_XE) * 2.0 * F_tot - LAMBDA_FP * N_FP
        dN_Xe_dt = (
            y_XE * 2.0 * F_tot
            - LAMBDA_XE * N_Xe
            - R_Xe_cap_th
        )

        # --- Mise à jour Euler ---

        n_fast += dn_fast_dt * dt
        n_th   += dn_th_dt * dt

        N_U235  += dN_U235_dt * dt
        N_U236  += dN_U236_dt * dt
        N_U237  += dN_U237_dt * dt
        N_U238  += dN_U238_dt * dt
        N_U239  += dN_U239_dt * dt
        N_Np239 += dN_Np239_dt * dt
        N_Pu239 += dN_Pu239_dt * dt
        N_Pu240 += dN_Pu240_dt * dt
        N_Pu241 += dN_Pu241_dt * dt
        N_Th232 += dN_Th232_dt * dt
        N_Th233 += dN_Th233_dt * dt
        N_Pa233 += dN_Pa233_dt * dt
        N_U233  += dN_U233_dt * dt
        N_FP    += dN_FP_dt * dt
        N_Xe    += dN_Xe_dt * dt

        # clamp
        n_fast = max(n_fast, 0.0)
        n_th   = max(n_th,   0.0)
        for N in [N_U235, N_U236, N_U237, N_U238, N_U239,
                  N_Np239, N_Pu239, N_Pu240, N_Pu241,
                  N_Th232, N_Th233, N_Pa233, N_U233, N_FP, N_Xe]:
            if N < 0:
                N = 0.0  # pas ultra propre mais évite les valeurs < 0 numériques

        # Puissance
        P = Q_FISSION * F_tot

        # stockage
        P_arr[k]       = P
        n_fast_arr[k]  = n_fast
        n_th_arr[k]    = n_th
        N_U235_arr[k]  = N_U235
        N_U236_arr[k]  = N_U236
        N_U237_arr[k]  = N_U237
        N_U238_arr[k]  = N_U238
        N_U239_arr[k]  = N_U239
        N_Np239_arr[k] = N_Np239
        N_Pu239_arr[k] = N_Pu239
        N_Pu240_arr[k] = N_Pu240
        N_Pu241_arr[k] = N_Pu241
        N_Th232_arr[k] = N_Th232
        N_Th233_arr[k] = N_Th233
        N_Pa233_arr[k] = N_Pa233
        N_U233_arr[k]  = N_U233
        N_FP_arr[k]    = N_FP
        N_Xe_arr[k]    = N_Xe

    E_tot = np.trapz(P_arr, t_arr)
    burnup = E_tot / mTot

    return {
        "time": t_arr,
        "power": P_arr,
        "n_fast": n_fast_arr,
        "n_thermal": n_th_arr,
        "N_U235": N_U235_arr,
        "N_U236": N_U236_arr,
        "N_U237": N_U237_arr,
        "N_U238": N_U238_arr,
        "N_U239": N_U239_arr,
        "N_Np239": N_Np239_arr,
        "N_Pu239": N_Pu239_arr,
        "N_Pu240": N_Pu240_arr,
        "N_Pu241": N_Pu241_arr,
        "N_Th232": N_Th232_arr,
        "N_Th233": N_Th233_arr,
        "N_Pa233": N_Pa233_arr,
        "N_U233": N_U233_arr,
        "N_FP": N_FP_arr,
        "N_Xe": N_Xe_arr,
        "burnup": burnup,
    }


# ------------------- CLASSES & TEST -------------------

class Fuel:
    def __init__(self):
        self.U235 = 3.0
        self.U238 = 97.0
        self.Pu239 = 0.0
        self.Th232 = 0.0


class FP:
    def __init__(self):
        self.Xe135 = 5.0
        self.FP    = 95.0


if __name__ == "__main__":
    fuel = Fuel()
    fp = FP()

    t_final    = 10.0
    n_th_init  = 1e10
    n_fa_init  = 0.0
    mTot       = 25.0

    res = reactorModel(
        fuelCompo=fuel,
        FPCompo=fp,
        t_final=t_final,
        n_th_init=n_th_init,
        n_fa_init=n_fa_init,
        mTot=mTot,
    )

    t = res["time"]

    # puissance
    plt.figure()
    plt.plot(t, res["power"])
    plt.xlabel("Time [s]")
    plt.ylabel("Power [W]")
    plt.grid(True)
    plt.title("Reactor power")
    plt.tight_layout()

    # neutrons
    plt.figure()
    plt.loglog(t, res["n_thermal"], label="n_thermal")
    plt.loglog(t, res["n_fast"], label="n_fast")
    plt.xlabel("Time [s]")
    plt.ylabel("Neutrons")
    plt.legend()
    plt.grid(True, which="both")
    plt.title("Neutron populations")
    plt.tight_layout()

    # quelques actinides
    plt.figure()
    plt.loglog(t, res["N_U233"], label="U-233")
    plt.loglog(t, res["N_U235"], label="U-235")
    plt.loglog(t, res["N_U236"], label="U-236")
    plt.loglog(t, res["N_U237"], label="U-237")
    plt.loglog(t, res["N_U238"], label="U-238")

    plt.loglog(t, res["N_Pu239"], label="Pu-239")
    plt.loglog(t, res["N_Pu240"], label="Pu-240")
    plt.loglog(t, res["N_Pu241"], label="Pu-241")


    plt.loglog(t, res["N_Xe"], label="Xe-135")
    plt.loglog(t, res["N_FP"], label="FP")

    plt.xlabel("Time [s]")
    plt.ylabel("Number of nuclei")
    plt.legend()
    plt.grid(True, which="both")
    plt.title("Main nuclides")
    plt.tight_layout()

    plt.show()
