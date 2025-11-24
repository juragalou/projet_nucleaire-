import numpy as np
import matplotlib.pyplot as plt

import molarMass as mM
import halfLife as hL
import crossSection as cS

# ------------------- CONSTANTES GLOBALES -------------------

NA = 6.022e23           # [1/mol]
V_CORE = 10.0           # m^3 (volume du coeur, donné dans l'énoncé)

# Energie par fission ~ 200 MeV
Q_FISSION = 200e6 * 1.602e-19   # [J/fission]

# Nombre moyen de neutrons par fission (ordre de grandeur PWR)
NU = 2.4

# Fraction de neutrons retardés (beta effectif)
BETA = 0.0065

# Ralentissement fast -> thermal : demi-vie = 5e-4 s
T_SLOW = 5e-4
LAMBDA_SLOW = np.log(2.0) / T_SLOW

# Paramètres de contrôle (barres + fuites) [1/s] (à choisir dans [0;20])
SIGMA_CTRL_FAST = 5.0
SIGMA_CTRL_TH = 5.0

# Energies typiques des neutrons
E_TH = 0.025      # eV (thermal)
E_FAST = 1e6      # eV (fast)
EV_TO_J = 1.602e-19
M_NEUTRON = 1.6749e-27  # kg

# ------------------- FONCTIONS UTILES -------------------

def decay_constant(species, transfo):
    """Retourne lambda = ln(2)/T1/2 pour un noyau donné."""
    T = hL.halfLife(species, transfo)
    if T <= 0:
        return 0.0
    return np.log(2.0) / T


# Constantes de décroissance utilisées
LAMBDA_FP    = decay_constant("FP",    "BetaMinus")
LAMBDA_XE    = decay_constant("Xe135", "BetaMinus")
LAMBDA_U239  = decay_constant("U239",  "BetaMinus")
LAMBDA_NP239 = decay_constant("Np239", "BetaMinus")
LAMBDA_TH233 = decay_constant("Th233", "BetaMinus")
LAMBDA_PA233 = decay_constant("Pa233", "BetaMinus")

# ------------------- FONCTION PRINCIPALE -------------------

def reactorModel(fuelCompo, FPCompo, t_final, n_th_init, n_fa_init, mTot):
    """
    Computes the evolution of the species in the reactor, and compute the power evolution in time
    :param fuelCompo: class
    :param FPCompo: class
            FPCompo.Xe135: percentage of FP that are considered as poison (in this project, we assume only
                Xe135 is produced as a poison)
            FPCompo.FP: percentage of FP that are considered "others fission products" (delayed neutron source)
    :param t_final: double
            final time of the simulation in seconds [s]
    :param n_th_init: double
            initial number of thermal neutrons
    :param n_fa_init: double
            initial number of fast neutrons
    :param mTot: double
            total mass of fuel at the initial state in [kg]
    :return results: dict
            time series for all main species + power + burnup
    """

    # -------- 1. Initialisation des quantités de combustible --------

    # masses des isotopes (en kg)
    m_U235 = mTot * fuelCompo.U235 / 100.0
    m_U238 = mTot * fuelCompo.U238 / 100.0
    m_Pu239 = mTot * fuelCompo.Pu239 / 100.0
    m_Th232 = mTot * fuelCompo.Th232 / 100.0

    # nombres de noyaux
    N_U235 = m_U235 / mM.molarMass("U235") * NA
    N_U238 = m_U238 / mM.molarMass("U238") * NA
    N_Pu239 = m_Pu239 / mM.molarMass("Pu239") * NA
    N_Th232 = m_Th232 / mM.molarMass("Th232") * NA  

    # Chaîne de conversion U238 -> U239 -> Np239 -> Pu239
    N_U239 = 0.0
    N_U237  = 0.0
    N_U236  = 0.0
    N_U233 = 0.0

    N_Np239 = 0.0

    N_Th233 = 0.0

    N_Pa233 = 0.0

    N_Pu240 = 0.0


    # Produits de fission initiaux
    N_FP = 0.0
    N_Xe = 0.0

    # Neutrons
    n_fast = n_fa_init
    n_th = n_th_init

    # Fraction de fission products qui sont du Xe-135 (donnée par FPCompo)
    y_XE = FPCompo.Xe135 / 100.0

    # -------- 2. Vitesses des neutrons et sections efficaces --------

    # vitesses (issues de E = 1/2 m v^2)
    v_th = np.sqrt(2.0 * E_TH * EV_TO_J / M_NEUTRON)
    v_fast = np.sqrt(2.0 * E_FAST * EV_TO_J / M_NEUTRON)

    
    def xs_th_fast(nuclide, transfo):
        """Retourne (sigma_th, sigma_fast) en m^2 à partir de crossSection (en barn)."""
        vals = cS.crossSection(nuclide, transfo, [E_TH, E_FAST])
        return vals[0] * 1e-28, vals[1] * 1e-28

    # U235 fission & capture (thermal + fast)
    sigma_U235_fis = cS.crossSection("U235", "Fission", [E_TH, E_FAST])
    sigma_U235_cap = cS.crossSection("U235", "Capture", [E_TH, E_FAST])
    sigma_U235_fis_th = sigma_U235_fis[0] * 1e-28
    sigma_U235_fis_fast = sigma_U235_fis[1] * 1e-28
    sigma_U235_cap_th = sigma_U235_cap[0] * 1e-28
    sigma_U235_cap_fast = sigma_U235_cap[1] * 1e-28

    # U238 fission (fast) & capture (th + fast)
    sigma_U238_fis = cS.crossSection("U238", "Fission", [E_TH, E_FAST])
    sigma_U238_cap = cS.crossSection("U238", "Capture", [E_TH, E_FAST])
    sigma_U238_fis_fast = sigma_U238_fis[1] * 1e-28
    sigma_U238_cap_th = sigma_U238_cap[0] * 1e-28
    sigma_U238_cap_fast = sigma_U238_cap[1] * 1e-28

    # Pu239 fission & capture (au début N_Pu239 = 0, mais va monter)
    sigma_Pu239_fis = cS.crossSection("Pu239", "Fission", [E_TH, E_FAST])
    sigma_Pu239_cap = cS.crossSection("Pu239", "Capture", [E_TH, E_FAST])
    sigma_Pu239_fis_th = sigma_Pu239_fis[0] * 1e-28
    sigma_Pu239_fis_fast = sigma_Pu239_fis[1] * 1e-28
    sigma_Pu239_cap_th = sigma_Pu239_cap[0] * 1e-28
    sigma_Pu239_cap_fast = sigma_Pu239_cap[1] * 1e-28

    # Xe135 capture thermique [m^2]
    sigma_Xe_cap_vals = cS.crossSection("Xe135", "Capture", [E_TH, E_FAST])
    sigma_Xe_cap_th = sigma_Xe_cap_vals[0] * 1e-28

    # coefficient effectif pour le poison au xénon :
    # R_cap_Xe = sigma_Xe * flux_th * N_Xe
    # flux_th = n_th * v_th / V_CORE
    # => R_cap_Xe = (sigma_Xe * v_th / V_CORE) * N_Xe * n_th
    kappa_xe = sigma_Xe_cap_th * v_th / V_CORE

    # -------- 3. Temps et stockage --------

    dt = 1e-4
    n_steps = int(t_final / dt)

    t_arr = np.zeros(n_steps)
    P_arr = np.zeros(n_steps)
    n_fast_arr = np.zeros(n_steps)
    n_th_arr = np.zeros(n_steps)

    N_U235_arr = np.zeros(n_steps)
    N_U238_arr = np.zeros(n_steps)
    N_U239_arr = np.zeros(n_steps)
    N_Np239_arr = np.zeros(n_steps)
    N_Pu239_arr = np.zeros(n_steps)
    N_FP_arr = np.zeros(n_steps)
    N_Xe_arr = np.zeros(n_steps)

    # -------- 4. Boucle temporelle (Euler explicite) --------

    for k in range(n_steps):

        t = k * dt
        t_arr[k] = t

        # --- 4.1 Flux des neutrons ---
        # flux [neutrons / (m^2 s)] = n * v / V
        flux_th = n_th * v_th / V_CORE
        flux_fast = n_fast * v_fast / V_CORE

        # --- 4.2 Taux de fission et captures physiques : σ * flux * N ---

        # --- 4.2 Taux de fission et captures physiques : σ * flux * N ---

        # U235 : fission + capture, séparées en th / fast
        R_U235_fis_th = sigma_U235_fis_th * flux_th * N_U235
        R_U235_fis_fast = sigma_U235_fis_fast * flux_fast * N_U235
        R_U235_fis = R_U235_fis_th + R_U235_fis_fast

        R_U235_cap_th = sigma_U235_cap_th * flux_th * N_U235
        R_U235_cap_fast = sigma_U235_cap_fast * flux_fast * N_U235
        R_U235_cap = R_U235_cap_th + R_U235_cap_fast

        # U238 : fission (fast seulement) + capture
        R_U238_fis_fast = sigma_U238_fis_fast * flux_fast * N_U238
        R_U238_fis = R_U238_fis_fast

        R_U238_cap_th = sigma_U238_cap_th * flux_th * N_U238
        R_U238_cap_fast = sigma_U238_cap_fast * flux_fast * N_U238
        R_U238_cap = R_U238_cap_th + R_U238_cap_fast

        # U239 -> Np239 (beta)
        R_U239_beta = LAMBDA_U239 * N_U239

        # Np239 -> Pu239 (beta)
        R_Np239_beta = LAMBDA_NP239 * N_Np239

        # Pu239 : fission + capture
        R_Pu239_fis_th = sigma_Pu239_fis_th * flux_th * N_Pu239
        R_Pu239_fis_fast = sigma_Pu239_fis_fast * flux_fast * N_Pu239
        R_Pu239_fis = R_Pu239_fis_th + R_Pu239_fis_fast

        R_Pu239_cap_th = sigma_Pu239_cap_th * flux_th * N_Pu239
        R_Pu239_cap_fast = sigma_Pu239_cap_fast * flux_fast * N_Pu239
        R_Pu239_cap = R_Pu239_cap_th + R_Pu239_cap_fast

        # Fissions totales (toutes énergies confondues)
        F_tot = R_U235_fis + R_U238_fis + R_Pu239_fis

        # Fissions / captures par groupe d'énergie (pour enlever les neutrons)
        R_fis_fast = R_U235_fis_fast + R_U238_fis_fast + R_Pu239_fis_fast
        R_fis_th   = R_U235_fis_th   + R_Pu239_fis_th

        R_cap_fast = R_U235_cap_fast + R_U238_cap_fast + R_Pu239_cap_fast
        R_cap_th   = R_U235_cap_th   + R_U238_cap_th   + R_Pu239_cap_th

        # Fissions totales
        F_tot = R_U235_fis + R_U238_fis + R_Pu239_fis

        # --- 4.3 Source de neutrons retardés via FP ---
        S_d = BETA * LAMBDA_FP * N_FP

        # --- 4.4 Équations différentielles ---

        # Neutrons rapides :
        # + ν F_tot (tous les neutrons de fission naissent rapides)
        # - fissions rapides (un neutron utilisé par fission)
        # - captures rapides
        # - pertes contrôle / fuites
        # - ralentissement vers thermiques
        # + neutrons retardés
        dn_fast_dt = (
            NU * F_tot
            - R_fis_fast
            - R_cap_fast
            - SIGMA_CTRL_FAST * n_fast
            - LAMBDA_SLOW * n_fast
            + S_d
        )

        # Neutrons thermiques :
        # + ralentissement des rapides
        # - fissions thermiques
        # - captures thermiques
        # - pertes contrôle
        # - absorption sur xénon
        dn_th_dt = (
            LAMBDA_SLOW * n_fast
            - R_fis_th
            - R_cap_th
            - SIGMA_CTRL_TH * n_th
            - kappa_xe * N_Xe * n_th
        )


        # Actinides :
        dN_U235_dt = -R_U235_fis - R_U235_cap
        dN_U238_dt = -R_U238_fis - R_U238_cap
        dN_U239_dt = +R_U238_cap - R_U239_beta
        dN_Np239_dt = +R_U239_beta - R_Np239_beta
        dN_Pu239_dt = +R_Np239_beta - R_Pu239_fis - R_Pu239_cap

        # FP : produits de fission (hors Xe) + décroissance
        dN_FP_dt = (1.0 - y_XE) * 2.0 * F_tot - LAMBDA_FP * N_FP

        # Xe135 : produit par fission + décroissance + capture
        dN_Xe_dt = (
            y_XE * 2.0 * F_tot
            - LAMBDA_XE * N_Xe
            - kappa_xe * N_Xe * n_th
        )

        # --- 4.5 Mise à jour (Euler explicite) ---

        n_fast += dn_fast_dt * dt
        n_th += dn_th_dt * dt
        N_U235 += dN_U235_dt * dt
        N_U238 += dN_U238_dt * dt
        N_U239 += dN_U239_dt * dt
        N_Np239 += dN_Np239_dt * dt
        N_Pu239 += dN_Pu239_dt * dt
        N_FP += dN_FP_dt * dt
        N_Xe += dN_Xe_dt * dt

        # éviter les valeurs négatives numériquement
        n_fast = max(n_fast, 0.0)
        n_th = max(n_th, 0.0)
        N_U235 = max(N_U235, 0.0)
        N_U238 = max(N_U238, 0.0)
        N_U239 = max(N_U239, 0.0)
        N_Np239 = max(N_Np239, 0.0)
        N_Pu239 = max(N_Pu239, 0.0)
        N_FP = max(N_FP, 0.0)
        N_Xe = max(N_Xe, 0.0)

        # --- 4.6 Puissance instantanée ---

        P = Q_FISSION * F_tot

        # stockage
        P_arr[k] = P
        n_fast_arr[k] = n_fast
        n_th_arr[k] = n_th
        N_U235_arr[k] = N_U235
        N_U238_arr[k] = N_U238
        N_U239_arr[k] = N_U239
        N_Np239_arr[k] = N_Np239
        N_Pu239_arr[k] = N_Pu239
        N_FP_arr[k] = N_FP
        N_Xe_arr[k] = N_Xe

    # -------- 5. Burnup simple --------

    E_tot = np.trapz(P_arr, t_arr)     # [J]
    burnup = E_tot / mTot              # [J/kg]

    # -------- 6. Résultats --------

    results = {
        "time": t_arr,
        "power": P_arr,
        "n_fast": n_fast_arr,
        "n_thermal": n_th_arr,
        "N_U235": N_U235_arr,
        "N_U238": N_U238_arr,
        "N_U239": N_U239_arr,
        "N_Np239": N_Np239_arr,
        "N_Pu239": N_Pu239_arr,
        "N_FP": N_FP_arr,
        "N_Xe": N_Xe_arr,
        "burnup": burnup,
    }

    return results


# ------------------- EXÉCUTION / TEST + GRAPHS -------------------

class Fuel:
    def __init__(self):
        self.U235 = 3
        self.U238 = 97
        self.Pu239 = 0
        self.Th232 = 0

class FP:
    def __init__(self):
        self.Xe135 = 5   # [%]
        self.FP = 95     # [%]


if __name__ == "__main__":
    fuel = Fuel()
    fp = FP()

    t_final = 100.0
    n_th_init = 1e10
    n_fa_init = 0.0
    mTot = 25.0

    res = reactorModel(
        fuelCompo=fuel,
        FPCompo=fp,
        t_final=t_final,
        n_th_init=n_th_init,
        n_fa_init=n_fa_init,
        mTot=mTot
    )

    t = res["time"]
    P = res["power"]
    n_fast = res["n_fast"]
    n_th = res["n_thermal"]
    N_U235 = res["N_U235"]
    N_U238 = res["N_U238"]
    N_Pu239 = res["N_Pu239"]
    N_Xe = res["N_Xe"]
    N_FP = res["N_FP"]

    # Nuclides (actinides + FP + Xe)
    plt.figure()
    plt.loglog(t, N_U235, label="U-235")
    plt.loglog(t, N_U238, label="U-238")
    plt.loglog(t, N_Pu239, label="Pu-239")
    plt.loglog(t, N_FP, label="FP (delayed)")
    plt.loglog(t, N_Xe, label="Xe-135")
    plt.xlabel("Time [s]")
    plt.ylabel("Number of nuclei")
    plt.legend()
    plt.grid(True, which="both")
    plt.title("Nuclide populations")
    plt.tight_layout()

    # Neutrons
    plt.figure()
    plt.loglog(t, n_th, label="n_thermal")
    plt.loglog(t, n_fast, label="n_fast")
    plt.xlabel("Time [s]")
    plt.ylabel("Number of neutrons")
    plt.legend()
    plt.grid(True, which="both")
    plt.title("Neutron populations")
    plt.tight_layout()

    # Power
    plt.figure()
    plt.plot(t, P)
    plt.xlabel("Time [s]")
    plt.ylabel("Power [W]")
    plt.grid(True)
    plt.title("Reactor power")
    plt.tight_layout()

    plt.show()
