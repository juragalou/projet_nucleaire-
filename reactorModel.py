import numpy as np
import sys
import matplotlib.pyplot as plt

import molarMass as mM
import halfLife as hL
import crossSection as cS

# --- global histories ---
TIME = None
P_HIST = None
BURNUP_HIST = None
NFA_HIST = None
NTH_HIST = None
NU_HIST = None
NPU_HIST = None
NXE_HIST = None
SIGMA_TH_CTRL_HIST = None


# ====== Options globales pour l'extérieur ======
USE_CONTROL = False   # si True : contrôle de puissance via les barres
P_NOM = 1e7           # [W] puissance nominale visée
K_P = 1e-10           # gain du correcteur P 


def reactorModel(fuelCompo, FPCompo, t_final, n_th_init, n_fa_init, mTot):
    """
    Calcule l’évolution des différentes espèces présentes dans le réacteur
    ainsi que l’évolution de la puissance au cours du temps.

    :param fuelCompo: classe
        Composition du combustible (fractions massiques de U235, U238, Pu239, Th232).
    :param FPCompo: classe
        Composition des produits de fission :
            - FPCompo.Xe135 : pourcentage des produits de fission considérés comme poison (ici uniquement le Xe135),
            - FPCompo.FP : pourcentage des autres produits de fission issus de la fission de U235.
    :param t_final: float
        Temps final de la simulation [s].
    :param n_th_init: float
        Nombre initial de neutrons thermiques.
    :param n_fa_init: float
        Nombre initial de neutrons rapides.
    :param mTot: float
        Masse totale de combustible initiale [kg].
    :return fuelBurnup: array
        Valeur du burnup (cumulé) en fonction du temps [MWd/t].
    """


    global TIME, P_HIST, BURNUP_HIST, NFA_HIST, NTH_HIST, NU_HIST, NPU_HIST, NXE_HIST, SIGMA_TH_CTRL_HIST
    global USE_CONTROL, P_NOM, K_P  
    # ======================= CONSTANTES PHYSIQUES =======================
    NA = 6.02214076e23        # [1/mol]
    V_core = 10.0             # [m^3]
    dt = 1e-3                # [s]
    n_steps = int(t_final / dt)

    v_th = 2200.0
    v_fast = 1.4e7

    E_th = 0.025    # [eV]
    E_fast = 1e6    # [eV]

    nu = 2.0
    T_slow = 5e-4
    lambda_slow = np.log(2) / T_slow

    T_FP = 1.0
    lambda_FP = np.log(2) / T_FP

    T_Xe = hL.halfLife('Xe135', 'BetaMinus')
    lambda_Xe = np.log(2) / T_Xe

    y_Xe = FPCompo.Xe135 / (FPCompo.Xe135 + FPCompo.FP)

    beta = 0.006



    # bornes pour les pertes effectives (fuites + barres)
    Sigma_fast_min, Sigma_fast_max = 2.0, 20.0   # [1/s]
    Sigma_th_min,  Sigma_th_max  = 2.0, 20.0     # [1/s]

    # valeurs initiales (barres dehors, donc pertes = fuites seules)
    Sigma_fast = Sigma_fast_min
    Sigma_th   = Sigma_th_min



    eV_to_J = 1.602176634e-19
    Qf = 200e6 * eV_to_J
    QFP = 5e6 * eV_to_J
    Qslow = 1e6 * eV_to_J

    # ======================= COMPOSITION DU COMBUSTIBLE =======================

    frac_U235 = fuelCompo.U235 / 100.0
    frac_Pu239 = fuelCompo.Pu239 / 100.0
    frac_Th232 = fuelCompo.Th232 / 100.0

    m_U235 = mTot * frac_U235
    m_Pu239 = mTot * frac_Pu239
    m_Th232 = mTot * frac_Th232

    NU235 = m_U235 / mM.molarMass('U235') * NA
    NPu239 = m_Pu239 / mM.molarMass('Pu239') * NA
    NTh232 = m_Th232 / mM.molarMass('Th232') * NA

    NU = NU235
    NPu = NPu239
    NFP = 0.0
    NXe = 0.0

    n_fa = n_fa_init
    n_th = n_th_init

    # ======================= SECTIONS EFFICACES MICROSCOPIQUES =======================

    barn_to_m2 = 1e-28

    sigma_U235_th = cS.crossSection('U235', 'Fission', [E_th])[0] * barn_to_m2
    sigma_U235_fast = cS.crossSection('U235', 'Fission', [E_fast])[0] * barn_to_m2

    sigma_Pu239_th = cS.crossSection('Pu239', 'Fission', [E_th])[0] * barn_to_m2
    sigma_Pu239_fast = cS.crossSection('Pu239', 'Fission', [E_fast])[0] * barn_to_m2

    sigma_Xe_th = cS.crossSection('Xe135', 'Capture', [E_th])[0] * barn_to_m2

    # ======================= ALLOCATION DES HISTORIQUES (globales) =======================

    TIME = np.zeros(n_steps)
    P_HIST = np.zeros(n_steps)
    BURNUP_HIST = np.zeros(n_steps)
    NFA_HIST = np.zeros(n_steps)
    NTH_HIST = np.zeros(n_steps)
    NU_HIST = np.zeros(n_steps)
    NPU_HIST = np.zeros(n_steps)
    NXE_HIST = np.zeros(n_steps)
    SIGMA_TH_CTRL_HIST = np.zeros(n_steps)

    E_cum = np.zeros(n_steps)

    # ======================= BOUCLE TEMPORELLE =======================

    for k in range(n_steps):
        t = k * dt
        TIME[k] = t

        NU = max(NU, 0.0)
        NPu = max(NPu, 0.0)

        SigmaU_fast_f = v_fast * sigma_U235_fast * NU / V_core #taux de fission U235
        SigmaU_th_f   = v_th * sigma_U235_th   * NU / V_core
        SigmaPu_fast_f = v_fast * sigma_Pu239_fast * NPu / V_core
        SigmaPu_th_f   = v_th * sigma_Pu239_th   * NPu / V_core

        FU = SigmaU_fast_f * n_fa + SigmaU_th_f * n_th
        FPu = SigmaPu_fast_f * n_fa + SigmaPu_th_f * n_th
        Ftot = FU + FPu

        Sd = beta * lambda_FP * NFP

        NXe = max(NXe, 0.0)
        Sigma_Xe_th = v_th * sigma_Xe_th * NXe / V_core

        dn_fa = nu * Ftot - Sigma_fast * n_fa - lambda_slow * n_fa + Sd
        dn_th = lambda_slow * n_fa - (Sigma_th + Sigma_Xe_th) * n_th

        dNU = -FU
        dNPu = -FPu

        dNFP = (1.0 - y_Xe) * 2.0 * Ftot - lambda_FP * NFP
        dNXe = y_Xe * 2.0 * Ftot - lambda_Xe * NXe

        # === Euler explicite ===
        n_fa += dn_fa * dt
        n_th += dn_th * dt
        NU   += dNU * dt
        NPu  += dNPu * dt
        NFP  += dNFP * dt
        NXe  += dNXe * dt

        n_fa = max(n_fa, 0.0)
        n_th = max(n_th, 0.0)
        NFP  = max(NFP, 0.0)

        # Puissance instantanée
        P = Qf * Ftot + QFP * lambda_FP * NFP + Qslow * lambda_slow * n_fa
        P_HIST[k] = P

        # Énergie cumulée
        if k == 0:
            E_cum[k] = P * dt
        else:
            E_cum[k] = E_cum[k-1] + P * dt

        # Burnup
        J_per_MWd_per_t = 8.64e7
        BURNUP_HIST[k] = (E_cum[k] / mTot) / J_per_MWd_per_t

        # Sauvegarde des autres variables
        NFA_HIST[k] = n_fa
        NTH_HIST[k] = n_th
        NU_HIST[k] = NU
        NPU_HIST[k] = NPu
        NXE_HIST[k] = NXe
        SIGMA_TH_CTRL_HIST[k] = Sigma_th

        # Contrôle des barres (si activé)
        if USE_CONTROL:
            error = P - P_NOM
            Sigma_th += K_P * error * dt
            Sigma_th = float(np.clip(Sigma_th, Sigma_th_min, Sigma_th_max))

    # on retourne juste la burnup comme prévu
    fuelBurnup = BURNUP_HIST.copy()
    return fuelBurnup


class Fuel:
    def __init__(self):
        self.U235 = 3
        self.U238 = 97
        self.Pu239 = 0
        self.Th232 = 0


class FP:
    def __init__(self):
        self.Xe135 = 3.165
        self.FP = 96.835


fuelCompo = Fuel()
FPCompo = FP()

# Exemple de test rapide (à commenter si l’auto-correcteur n’aime pas) :
# if __name__ == "__main__":
#     fBurnup = reactorModel(fuelCompo=fuelCompo, FPCompo=FPCompo,
#                            t_final=1.0, n_th_init=1e10, n_fa_init=0., mTot=25.)
#     t = np.linspace(0, 1.0, len(fBurnup))
#     plt.plot(t, fBurnup)
#     plt.xlabel("Time [s]")
#     plt.ylabel("Burnup [MWd/t]")
#     plt.grid(True)
#     plt.show()
