def molarMass(X):
    """
    Database: http://wwwndc.jaea.go.jp/NuC/
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        i.e.: for the Uranium 235, X = 'U235'; for a neutron X = 'n'
    :return: double
        molar mass of the nucleon, or nuclide X in [kg/mol]
    """

    # Your work : add all species

    if X == 'Th232':
        M = 232.038055760/1000
    elif X == 'Th233':
        M = 233.041586541/1000
    elif X == 'Np239':
        M = 239.052940487/1000
    elif X == 'Pa233':
        M = 233.040248815/1000
    elif X == 'Pu239':
        M = 239.052164844/1000
    elif X == 'Pu240':
        M = 240.053815008/1000
    elif X == 'Pu241':
        M = 241.056852919/1000
    elif X == 'U233':
        M = 233.039636574/1000
    elif X == 'U235':
        M = 235.043931368/1000
    elif X == 'U236':
        M = 236.045569468/1000
    elif X == 'U237':
        M = 237.048731636/1000
    elif X == 'U238':
        M = 238.050789466/1000
    elif X == 'U239':
        M = 239.054294518/1000
    elif X == 'neutron':
        M = 1.00866/1000
    else:
        print('\n ---------------- WARNING ----------------- \n No molar mass for species (', X, ')')
        M = 0

    return M

