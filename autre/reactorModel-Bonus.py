
import numpy as np
import matplotlib.pyplot as plt

import molarMass as mM
import halfLife as hL
import crossSection as cS

N_a = 6.02214076e23
m_n = 1.6749e-27 
V_reactor = 10 
E_fa = 1e6 
E_th = 0.025 
barn_number = 1e-24 

sigma_rod_min = 20
sigma_rod_max = 100

E_decay = 4e6
E_fission = 2e8
eV_to_Joules = 1.60218e-19
Power_reactor = 3e9

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
    
    time_step = 1e-4
    time = np.arange(0,t_final,time_step)
    sigma_rod = sigma_rod_max #no speed of control rod
    dsigma_rod = 5e-1 * sigma_rod_max
    sigma_rod_list = np.zeros(len(time))
    sigma_rod_list[0] = sigma_rod 
    fuelBurnup = 0
    Ntot = mTot / ((fuelCompo.U235/100 * mM.molarMass('U235')) + (fuelCompo.U238/100 * mM.molarMass('U238')))
    v_th = np.sqrt( 2*(E_th * eV_to_Joules) / (mM.molarMass("neutron")/N_a) )
    v_fa = np.sqrt( 2*(E_fa * eV_to_Joules) / (mM.molarMass("neutron")/N_a) )

    n_th = np.zeros(len(time))
    n_fa = np.zeros(len(time))
    flux_n_th = np.zeros(len(time))
    flux_n_fa = np.zeros(len(time))
    U235 = np.zeros(len(time))
    U236 = np.zeros(len(time))
    U237 = np.zeros(len(time))
    U238 = np.zeros(len(time))
    U239 = np.zeros(len(time))
    Np239 = np.zeros(len(time))
    Pu239 = np.zeros(len(time))
    Pu240 = np.zeros(len(time))
    Xe135 = np.zeros(len(time))
    Pu239 = np.zeros(len(time))
    Pu241 = np.zeros(len(time))
    Th232 = np.zeros(len(time))
    Th233 = np.zeros(len(time))
    Pa233 = np.zeros(len(time))
    U233 = np.zeros(len(time))
    FP_star = np.zeros(len(time))
    Power = np.zeros(len(time))

    U235[0] = fuelCompo.U235/100 * Ntot * N_a
    U238[0] = fuelCompo.U238/100 * Ntot * N_a
    Pu239[0] = fuelCompo.Pu239/100 * Ntot * N_a
    Th232[0] = fuelCompo.Th232/100 * Ntot * N_a
    n_th[0] = n_th_init 
    n_fa[0] = n_fa_init
    flux_n_th[0] = (n_th_init * v_th) / V_reactor
    flux_n_fa[0] = (n_fa_init * v_fa) / V_reactor

    sigma_U235_fis_th = cS.crossSection('U235','Fission', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_U236_fis_th = cS.crossSection('U236','Fission', [E_th, E_fa])[0] * barn_number * 0.0001 
    sigma_U237_fis_th = cS.crossSection('U237','Fission', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_U238_fis_th = cS.crossSection('U238','Fission', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_U239_fis_th = cS.crossSection('U239','Fission', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Pu239_fis_th = cS.crossSection('Pu239','Fission', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Pu240_fis_th = cS.crossSection('Pu240', 'Fission', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Pu241_fis_th = cS.crossSection('Pu241','Fission', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Np_239_fis_th = cS.crossSection('Np239', 'Fission', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Th233_fis_th = cS.crossSection('Th233', 'Fission', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_U233_fis_th = cS.crossSection('U233', 'Fission', [E_th, E_fa])[0] * barn_number * 0.0001
    
    sigma_U235_fis_fa = cS.crossSection('U235','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_U236_fis_fa = cS.crossSection('U236','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_U237_fis_fa = cS.crossSection('U237','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_U238_fis_fa = cS.crossSection('U238','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_U239_fis_fa = cS.crossSection('U239','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Pu239_fis_fa = cS.crossSection('Pu239','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Np239_fis_fa = cS.crossSection('Np239','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Pu240_fis_fa = cS.crossSection('Pu240','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Pu241_fis_fa = cS.crossSection('Pu241','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Th232_fis_fa = cS.crossSection('Th232','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Th233_fis_fa = cS.crossSection('Th233','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Pa233_fis_fa = cS.crossSection('Pa233','Fission', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_U233_fis_fa = cS.crossSection('U233','Fission', [E_th, E_fa])[1] * barn_number * 0.0001

    sigma_U235_cap_th = cS.crossSection('U235','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_U236_cap_th = cS.crossSection('U236','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_U237_cap_th = cS.crossSection('U237','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_U238_cap_th = cS.crossSection('U238','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_U239_cap_th = cS.crossSection('U239','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Pu239_cap_th = cS.crossSection('Pu239','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Pu240_cap_th = cS.crossSection('Pu240','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Pu241_cap_th = cS.crossSection('Pu241','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Np239_cap_th = cS.crossSection('Np239','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Xe135_cap_th = cS.crossSection('Xe135','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Th232_cap_th = cS.crossSection('Th232','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Th233_cap_th = cS.crossSection('Th233','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_Pa233_cap_th = cS.crossSection('Pa233','Capture', [E_th, E_fa])[0] * barn_number * 0.0001
    sigma_U233_cap_th = cS.crossSection('U233','Capture', [E_th, E_fa])[0] * barn_number * 0.0001

    sigma_U235_cap_fa = cS.crossSection('U235','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_U236_cap_fa = cS.crossSection('U236','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_U237_cap_fa = cS.crossSection('U237','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_U238_cap_fa = cS.crossSection('U238','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_U239_cap_fa = cS.crossSection('U239','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Pu239_cap_fa = cS.crossSection('Pu239','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Pu240_cap_fa = cS.crossSection('Pu240','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Pu241_cap_fa = cS.crossSection('Pu241', 'Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Np239_cap_fa = cS.crossSection('Np239','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Xe135_cap_fa = cS.crossSection('Xe135','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Th232_cap_fa = cS.crossSection('Th232','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Th233_cap_fa = cS.crossSection('Th233','Capture', [E_th, E_fa])[1] * barn_number * 0.0001
    sigma_Pa233_cap_fa = cS.crossSection('Pa233','Capture', [E_th, E_fa])[1] * barn_number * 0.0001

    half_life_U237 = hL.halfLife('U237','BetaMinus')
    half_life_U239 = hL.halfLife('U239','BetaMinus')
    half_life_Np239 = hL.halfLife('Np239','BetaMinus')
    half_life_Xe135 = hL.halfLife('Xe135','BetaMinus')
    half_life_Th233 = hL.halfLife('Th233','BetaMinus')
    half_life_Pa233 = hL.halfLife('Pa233','BetaMinus')
    lambda_U237 = np.log(2)/half_life_U237
    lambda_U239 = np.log(2)/half_life_U239
    lambda_Np239 = np.log(2)/half_life_Np239
    lambda_Xe135 = np.log(2)/half_life_Xe135
    lambda_Th233 = np.log(2)/half_life_Th233
    lambda_Pa233 = np.log(2)/half_life_Pa233
    lambda_n_fa = np.log(2)/(5e-4)
    lambda_FP_star = np.log(2)
      
    for i in range(len(time)-1) :
        
        U235_fis_th = sigma_U235_fis_th * flux_n_th[i]*U235[i] 
        U235_fis_fa = sigma_U235_fis_fa * flux_n_fa[i]*U235[i] 
        U235_cap_fa = sigma_U235_cap_fa * flux_n_fa[i]*U235[i]
        U235_cap_th = sigma_U235_cap_th * flux_n_th[i]*U235[i]      

        U236_cap_th = sigma_U236_cap_th * flux_n_th[i]*U236[i]       
        U236_cap_fa = sigma_U236_cap_fa * flux_n_fa[i]*U236[i]
        U236_fis_fa = sigma_U236_fis_fa * flux_n_fa[i]*U236[i]
        U236_fis_th = sigma_U236_fis_th * flux_n_th[i]*U236[i]

        U237_fis_th = sigma_U237_fis_th * flux_n_th[i]*U237[i] 
        U237_cap_th = sigma_U237_cap_th * flux_n_th[i]*U237[i]      
        U237_fis_fa = sigma_U237_fis_fa * flux_n_fa[i]*U237[i] 
        U237_cap_fa = sigma_U237_cap_fa * flux_n_fa[i]*U237[i]
        U237_beta = lambda_U237 * U237[i] 

        U238_fis_th = sigma_U238_fis_th * flux_n_th[i]*U238[i] 
        U238_cap_th = sigma_U238_cap_th * flux_n_th[i]*U238[i]      
        U238_fis_fa = sigma_U238_fis_fa * flux_n_fa[i]*U238[i] 
        U238_cap_fa = sigma_U238_cap_fa * flux_n_fa[i]*U238[i]
    
        U239_fis_th = sigma_U239_fis_th * flux_n_th[i]*U239[i]
        U239_cap_th = sigma_U239_cap_th * flux_n_th[i]*U239[i]
        U239_fis_fa = sigma_U239_fis_fa * flux_n_fa[i]*U239[i] 
        U239_cap_fa = sigma_U239_cap_fa * flux_n_fa[i]*U239[i]
        U239_beta = lambda_U239 * U239[i]
  
        Pu239_fis_th = sigma_Pu239_fis_th * flux_n_th[i]*Pu239[i] 
        Pu239_cap_th = sigma_Pu239_cap_th * flux_n_th[i]*Pu239[i]     
        Pu239_fis_fa = sigma_Pu239_fis_fa * flux_n_fa[i]*Pu239[i] 
        Pu239_cap_fa = sigma_Pu239_cap_fa * flux_n_fa[i]*Pu239[i]
         
        Pu240_cap_th = sigma_Pu240_cap_th * flux_n_th[i]*Pu240[i]    
        Pu240_cap_fa = sigma_Pu240_cap_fa * flux_n_fa[i]*Pu240[i]
        Pu240_fis_fa = sigma_Pu240_fis_fa * flux_n_fa[i]*Pu240[i]
        Pu240_fis_th = sigma_Pu240_fis_th * flux_n_th[i]*Pu240[i]
        
        Pu241_fis_th = sigma_Pu241_fis_th * flux_n_th[i]*Pu241[i] 
        Pu241_fis_fa = sigma_Pu241_fis_fa * flux_n_fa[i]*Pu241[i]
        Pu241_cap_th = sigma_Pu241_cap_th * flux_n_th[i]*Pu241[i]
        Pu241_cap_fa = sigma_Pu241_cap_fa * flux_n_fa[i]*Pu241[i]
        
        Np239_cap_th = sigma_Np239_cap_th * flux_n_th[i]*Np239[i]      
        Np239_cap_fa = sigma_Np239_cap_fa * flux_n_fa[i]*Np239[i]
        Np239_fis_fa = sigma_Np239_fis_fa * flux_n_fa[i]*Np239[i]
        Np239_fis_th = sigma_Np_239_fis_th * flux_n_th[i]*Np239[i]
        Np239_beta = lambda_Np239 * Np239[i]

        Th232_fis_fa = sigma_Th232_fis_fa * flux_n_fa[i]*Th232[i]
        Th232_cap_th = sigma_Th232_cap_th * flux_n_th[i]*Th232[i]
        Th232_cap_fa = sigma_Th232_cap_fa * flux_n_fa[i]*Th232[i]

        Th233_fis_th = sigma_Th233_fis_th * flux_n_th[i]*Th233[i]
        Th233_cap_th = sigma_Th233_cap_th * flux_n_th[i]*Th233[i]
        Th233_fis_fa = sigma_Th233_fis_fa * flux_n_fa[i]*Th233[i]
        Th233_cap_fa = sigma_Th233_cap_fa * flux_n_fa[i]*Th233[i]
        Th233_beta = lambda_Th233 * Th233[i]

        Pa233_fis_fa = sigma_Pa233_fis_fa * flux_n_fa[i]*Pa233[i]
        Pa233_cap_fa = sigma_Pa233_cap_fa * flux_n_fa[i]*Pa233[i]
        Pa233_cap_th = sigma_Pa233_cap_th * flux_n_th[i]*Pa233[i]
        Pa233_beta = lambda_Pa233 * Pa233[i]

        U233_fis_fa = sigma_U233_fis_fa * flux_n_fa[i]*U233[i]
        U233_fis_th = sigma_U233_fis_th * flux_n_th[i]*U233[i]
        U233_cap_th = sigma_U233_cap_th * flux_n_th[i]*U233[i]

        Xe135_cap_th = sigma_Xe135_cap_th * flux_n_th[i]*Xe135[i]
        Xe135_cap_fa = sigma_Xe135_cap_fa * flux_n_fa[i]*Xe135[i]
        Xe135_beta = lambda_Xe135 * Xe135[i]
        
        dU235 = -U235_fis_th - U235_cap_th - U235_fis_fa
        dU236 = -U236_cap_th - U236_cap_fa - U236_fis_fa + U235_cap_th
        dU237 = -U237_cap_th - U237_fis_fa - U237_cap_fa - U237_beta + U236_cap_th + U236_cap_fa
        dU238 = -U238_cap_th - U238_fis_fa - U238_cap_fa + U237_cap_th + U237_cap_fa
        dU239 = -U239_fis_th - U239_cap_th - U239_fis_fa - U239_cap_fa - U239_beta + U238_cap_th + U238_cap_fa
        dNp239 = -Np239_cap_th - Np239_cap_fa - Np239_fis_fa - Np239_beta + U239_beta
        dPu239 = -Pu239_fis_th - Pu239_cap_th - Pu239_fis_fa + Np239_beta
        dPu240 = -Pu240_cap_th - Pu240_fis_fa + Pu239_cap_th
        dPu241 = -Pu241_fis_fa - Pu241_fis_th - Pu241_cap_th + Pu240_cap_th
        dTh232 = -Th232_fis_fa - Th232_cap_fa - Th232_cap_th
        dTh233 = -Th233_beta - Th233_fis_fa - Th233_fis_th - Th233_cap_th - Th233_cap_fa + Th232_cap_fa + Th232_cap_fa
        dPa233 = -Pa233_beta - Pa233_fis_fa - Pa233_cap_fa - Pa233_cap_th + Th233_beta
        dU233 = -U233_fis_th - U233_fis_fa - U233_cap_th + Pa233_beta

        n_fis_th = U235_fis_th + U239_fis_th + Pu239_fis_th + Pu241_fis_th + Th233_fis_th + U233_fis_th
        n_cap_th = U235_cap_th + U236_cap_th + U237_cap_th + U238_cap_th + U239_cap_th + Np239_cap_th + Pu239_cap_th + Pu240_cap_th + Pu241_cap_th + Xe135_cap_th + Th232_cap_th + Th233_cap_th + Pa233_cap_th
        n_fis_fa = U235_fis_fa + U236_fis_fa + U237_fis_fa + U238_fis_fa + U239_fis_fa + Np239_fis_fa + Pu239_fis_fa + Pu240_fis_fa + Pu241_fis_fa + Th232_fis_fa + Th233_fis_fa + Pa233_fis_fa + U233_fis_fa
        n_cap_fa = U236_cap_fa + U237_cap_fa + U238_cap_fa + U239_cap_fa + Np239_cap_fa + Xe135_cap_fa + Th233_cap_fa + Pa233_cap_fa

        n_th_rod = sigma_rod_list[i] * n_th[i] 
        n_fa_rod = sigma_rod_list[i] * n_fa[i] 

        dn_th = lambda_n_fa * n_fa[i] - n_cap_th - n_fis_th - n_th_rod 
        dn_fa = -n_fa[i] * lambda_n_fa + n_fis_fa - n_cap_fa + 2 * n_fis_th + lambda_FP_star * FP_star[i] - n_fa_rod

        dXe135 = -Xe135_beta - Xe135_cap_th - Xe135_cap_fa + 2 * n_fis_th * FPCompo.Xe135/100 + 2 * n_fis_fa * FPCompo.Xe135/100
        dFP_star = -lambda_FP_star * FP_star[i] + 2 * n_fis_th * FPCompo.FP/100 + 2 * n_fis_fa * FPCompo.FP/100

        U235[i+1] = U235[i] + time_step*dU235
        U236[i+1] = U236[i] + time_step*dU236
        U237[i+1] = U237[i] + time_step*dU237
        U238[i+1] = U238[i] + time_step*dU238
        U239[i+1] = U239[i] + time_step*dU239
        Np239[i+1] = Np239[i] + time_step*dNp239
        Pu239[i+1] = Pu239[i] + time_step*dPu239
        Pu240[i+1] = Pu240[i] + time_step*dPu240
        Pu241[i+1] = Pu241[i] + time_step*dPu241
        Xe135[i+1] = Xe135[i] + time_step*dXe135
        Th232[i+1] = Th232[i] + time_step*dTh232
        Th233[i+1] = Th233[i] + time_step*dTh233
        Pa233[i+1] = Pa233[i] + time_step*dPa233
        U233[i+1] = U233[i] + time_step*dU233
        n_th[i+1] = n_th[i] + time_step*dn_th
        n_fa[i+1] = n_fa[i] + time_step*dn_fa
        FP_star[i+1] = FP_star[i] + time_step*dFP_star 

        flux_n_th[i+1] = (n_th[i] * v_th) / V_reactor
        flux_n_fa[i+1] = (n_fa[i] * v_fa) / V_reactor

        Power_fission = E_fission * eV_to_Joules * (n_fis_fa + n_fis_th) / time_step 
        Power_decay_FP = E_decay * eV_to_Joules * lambda_FP_star * FP_star[i] / time_step 
        Power_decay_n = (E_fa - E_th) * eV_to_Joules * lambda_n_fa* n_fa[i] / time_step
        Power[i+1] = (Power_fission + Power_decay_n + Power_decay_FP) / (1e6)

        K = 0.01
        deriv_sigma = K * (Power[i+1] - Power_reactor) / time_step
        if deriv_sigma > dsigma_rod:
            deriv_sigma = dsigma_rod
        elif deriv_sigma < -dsigma_rod:
            deriv_sigma = -dsigma_rod

        sigma_rod_next = sigma_rod_list[i] + deriv_sigma*time_step

        if sigma_rod_next > sigma_rod_max:
            sigma_rod_next = sigma_rod_max
        elif sigma_rod_next < sigma_rod_min:
            sigma_rod_next = sigma_rod_min

        sigma_rod_list[i+1] = sigma_rod_next

        fuelBurnup += (U235_fis_th + U235_fis_fa) / U235[0] + (U238_fis_th + U238_fis_fa) / U238[0]
        
    print(U233)
    # plt.loglog(time,U235,label="U235")
    # plt.loglog(time,U236,label="U236")
    # plt.loglog(time,U237,label="U237")
    # plt.loglog(time,U238,label="U238")
    # plt.loglog(time,U239,label="U239")
    # plt.loglog(time,Np239,label="Np239")
    # plt.loglog(time,Pu239,label="Pu239")
    # plt.loglog(time,Pu240,label="Pu240")
    # plt.loglog(time,Pu241,label="Pu241")
    plt.loglog(time,FP_star,label='FP*')
    plt.loglog(time,Xe135,label="Xe135")
    plt.loglog(time,Th232,label="Th232")
    plt.loglog(time,Th233,label="Th233")
    plt.loglog(time,Pa233,label='Pa233')
    plt.loglog(time,U233,label="U233")
    plt.ylabel("Number of atoms")
    plt.xlabel("Time [s]")
    plt.legend()
    plt.grid()
    plt.show()

    plt.loglog(time,n_th,label='n_th')
    plt.loglog(time,n_fa,label='n_fa')
    plt.xlabel("Time [s]")
    plt.ylabel("Number of neutrons")
    plt.legend()
    plt.grid()
    plt.show()
    
    plt.plot(time, Power)
    plt.xlabel("Time [s]")
    plt.ylabel("Power [W]")
    plt.grid()
    plt.show()

    plt.plot(time, sigma_rod_list)
    plt.xlabel('Time [s]')
    plt.ylabel('Cross Section of Control Bars [barn]')
    plt.grid()
    plt.show()
    
    return fuelBurnup

class Fuel:
    def __init__(self):
        self.U235 = 3
        self.U238 = 93
        self.Pu239 = 2
        self.Th232 = 2

class FP:
    def __init__(self):
        self.Xe135 = 5
        self.FP = 95

fuelCompo = Fuel()
FPCompo = FP()

fBurnup = reactorModel(fuelCompo=fuelCompo, FPCompo=FPCompo, t_final=100, n_th_init=1e10, n_fa_init=0, mTot=25.)

