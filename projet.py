import numpy as np
import matplotlib.pyplot as plt

import reactorModel as rm


def run_case(label, use_xe=True, use_control=False, t_final= 100.0):
    # Configurer la composition du combustible
    fuel = rm.Fuel()
    FPc = rm.FP()

    # Cas sans Xe : on met tout en "autres FP"
    if not use_xe:
        FPc.Xe135 = 0
        FPc.FP = 100

    # Contrôle de puissance
    rm.USE_CONTROL = use_control
    rm.P_NOM = 1e7
    rm.K_P = 1e-10

    # Lancer la simulation
    _ = rm.reactorModel(
        fuelCompo=fuel,
        FPCompo=FPc,
        t_final=t_final,
        n_th_init=1e10,
        n_fa_init=0.,
        mTot=25.
    )

    # Récupération des historiques globaux
    t = rm.TIME
    P = rm.P_HIST
    burnup = rm.BURNUP_HIST
    NXe = rm.NXE_HIST
    n_fa = rm.NFA_HIST
    n_th = rm.NTH_HIST
    NU = rm.NU_HIST
    NPu = rm.NPU_HIST
    Sigma_th_ctrl = rm.SIGMA_TH_CTRL_HIST

    return t, P, burnup, NXe, n_fa, n_th, NU, NPu, Sigma_th_ctrl, label


if __name__ == "__main__":
    t_final = 100.0

    # Cas 1 : sans Xe, sans contrôle
    res1 = run_case("Sans Xe, sans contrôle", use_xe=False, use_control=False, t_final=t_final)
    # Cas 2 : avec Xe, sans contrôle
    res2 = run_case("Avec Xe, sans contrôle", use_xe=True, use_control=False, t_final=t_final)
    # Cas 3 : avec Xe, avec contrôle
    res3 = run_case("Avec Xe + contrôle", use_xe=True, use_control=True, t_final=t_final)

    # Déballage
    t1, P1, B1, Xe1, nf1, nt1, NU1, NPu1, Sctrl1, lab1 = res1
    t2, P2, B2, Xe2, nf2, nt2, NU2, NPu2, Sctrl2, lab2 = res2
    t3, P3, B3, Xe3, nf3, nt3, NU3, NPu3, Sctrl3, lab3 = res3

    # ===== 1. Puissance P(t) =====
    plt.figure()
    plt.plot(t1, P1, label=lab1)
    plt.plot(t2, P2, label=lab2)
    plt.plot(t3, P3, label=lab3)
    plt.xlabel("Temps [s]")
    plt.ylabel("Puissance [W]")
    plt.grid(True)
    plt.legend()
    plt.title("Évolution de la puissance du réacteur")
    plt.tight_layout()

    # ===== 2. Xe-135(t) pour les cas avec Xe =====
    plt.figure()
    plt.plot(t2, Xe2, label="Xe135 – sans contrôle")
    plt.plot(t3, Xe3, label="Xe135 – avec contrôle")
    plt.xlabel("Temps [s]")
    plt.ylabel("N_Xe [atomes]")
    plt.grid(True)
    plt.legend()
    plt.title("Accumulation de Xe-135")
    plt.tight_layout()

    # ===== 3. Neutrons rapides / thermiques (cas avec Xe + contrôle) =====
    plt.figure()
    plt.plot(t3, nf3, label="Neutrons rapides")
    plt.plot(t3, nt3, label="Neutrons thermiques")
    plt.xlabel("Temps [s]")
    plt.ylabel("Population de neutrons")
    plt.grid(True)
    plt.legend()
    plt.title("Évolution des populations de neutrons (cas Xe + contrôle)")
    plt.tight_layout()

    # ===== 4. U235 et Pu239 (cas Xe + contrôle) =====
    plt.figure()
    plt.plot(t3, NU3, label="U235")
    plt.plot(t3, NPu3, label="Pu239")
    plt.xlabel("Temps [s]")
    plt.ylabel("Nombre d'atomes")
    plt.grid(True)
    plt.legend()
    plt.title("Consommation du combustible fissile")
    plt.tight_layout()

    # ===== 5. Burnup cumulée =====
    plt.figure()
    plt.plot(t1, B1, label=lab1)
    plt.plot(t2, B2, label=lab2)
    plt.plot(t3, B3, label=lab3)
    plt.xlabel("Temps [s]")
    plt.ylabel("Burnup [MWd/t]")
    plt.grid(True)
    plt.legend()
    plt.title("Burnup cumulée du combustible")
    plt.tight_layout()

    # ===== 6. Sigma_th_ctrl(t) (cas contrôle) =====
    plt.figure()
    plt.plot(t3, Sctrl3)
    plt.xlabel("Temps [s]")
    plt.ylabel("Σ_th_ctrl [1/s]")
    plt.grid(True)
    plt.title("Action des barres de contrôle (cas Xe + contrôle)")
    plt.tight_layout()

    plt.show()
