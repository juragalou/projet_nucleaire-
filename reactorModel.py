import numpy as np
import sys
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

# Paramètres de contrôle (barres + fuites) [1/s] (a choisir dans un arange [0;20])
SIGMA_CTRL_FAST = 5.0
SIGMA_CTRL_TH = 5.0

# Fraction de produits de fission qui sont du Xe-135
#y_XE = 0.065    # typique, tu peux ajuster

# ------------------- FONCTIONS UTILES -------------------

def decay_constant(species, transfo):
    """Retourne lambda = ln(2)/T1/2 pour un noyau donné."""
    T = hL.halfLife(species, transfo)
    if T <= 0:
        return 0.0
    return np.log(2.0) / T


# Constantes de décroissance pour FP et Xe
LAMBDA_FP = decay_constant("FP", "BetaMinus")     # ~1 s^-1
LAMBDA_XE = decay_constant("Xe135", "BetaMinus")  # ~9.14 h

# ------------------- FONCTION PRINCIPALE -------------------

def reactorModel(fuelCompo, FPCompo, t_final, n_th_init, n_fa_init, mTot):

    """
    Computes the evolution of the species in the reactor, and compute the power evolution in time
    :param fuelCompo: class
    :param FPCompo: class
            FPCompo.Xe135: percentage of FP that are considered as poison (in this project, we assume only
                Xe135 is produced as a poison)
            FPCompo.FP: percentage of FP that are considered "others fission products" of U235 (assume only one isotope
                different than Xe135 is produced by the fission)
    :param t_final: double
            final time of the simulation in seconds [s]
    :param n_th_init: double
            initial number of thermal neutrons
    :param n_fa_init: double
            initial number of fast neutrons
    :param mTot: double
            total number of kg of fuel at the initial state in [kg]
    :return fuelBurnup: array
            depending on the fuel, returns the burnup
    """
    # -------- 1. Initialisation des quantités de combustible --------

    # masses molaires [kg/mol]
    M_U235 = mM.molarMass("U235")
    M_U238 = mM.molarMass("U238")

    # masses des isotopes (en kg)
    m_U235 = mTot * fuelCompo.U235 / 100.0
    m_U238 = mTot * fuelCompo.U238 / 100.0

    # nombres de noyaux
    N_U235 = m_U235 / M_U235 * NA
    N_U238 = m_U238 / M_U238 * NA

    # pour info / extension : Pu, Th etc. si tu veux les ajouter plus tard
    N_Pu239 = 0.0

    # Produits de fission initiaux (on peut commencer à 0 pour simplifier)
    N_FP = 0.0
    N_Xe = 0.0

    # Neutrons
    n_fast = n_fa_init
    n_th = n_th_init
    
    # -------- 2. Temps et stockage --------

    dt = 1e-4
    n_steps = int(t_final / dt)


    return fuelBurnup




# ------------------- EXÉCUTION DU MODÈLE -------------------
class Fuel:
    def __init__(self):
        self.U235 = 3
        self.U238 = 97
        self.Pu239 = 0
        self.Th232 = 0

class FP:
    def __init__(self):
        self.Xe135 = 5
        self.FP = 95

fuelCompo = Fuel()
FPCompo = FP()


# Try the function

# Mm = mM.molarMass('U235')
# hl = hL.halfLife(X='Pa233', Transfo='BetaMinus')
# cs = cS.crossSection(X='Pu240', Transfo='Fission', E_neutron=np.logspace(-5,6,10000).tolist())
#
# print(Mm, hl, cs)

# fBurnup = reactorModel(fuelCompo=fuelCompo, FPCompo=FPCompo, t_final=200., n_th_init=1e+10, n_fa_init=0., mTot=25.)