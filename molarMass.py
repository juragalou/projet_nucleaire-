'''
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
        M = 0
    elif X == 'Th233':
        M = 0
    else:
        print('\n ---------------- WARNING ----------------- \n No molar mass for species (', X, ')')
        M = 0

    return M

# Try your function:
Mm = molarMass('U235')
print(Mm)# molarMass.py
'''


################################################################################ Database of molar masses in [kg/mol]
MOLAR_MASS_DB = {
    # neutron
    "n": 1.008665e-3,  # kg/mol

    # Thorium
    "Th232": 232.0381e-3,
    "Th233": 233.0403e-3,

    # Protactinium
    "Pa233": 233.0407e-3,

    # Uranium
    "U233": 233.0396e-3,
    "U235": 235.0439e-3,
    "U236": 236.0456e-3,
    "U237": 237.0482e-3,
    "U238": 238.0508e-3,
    "U239": 239.0543e-3,

    # Neptunium
    "Np239": 239.0529e-3,

    # Plutonium
    "Pu239": 239.0522e-3,
    "Pu240": 240.0538e-3,

    # Xenon
    "Xe135": 134.9072e-3,

    # Effective fictitious fission product (for delayed neutrons)
    "FP": 100.0e-3,  # arbitrary ~100 g/mol
}


def molarMass(X):
    """
    Database: http://wwwndc.jaea.go.jp/NuC/
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        i.e.: for Uranium-235, X = 'U235'; for a neutron X = 'n'
    :return: double
        molar mass of the nucleon or nuclide X in [kg/mol]
    """

    if X in MOLAR_MASS_DB:
        return MOLAR_MASS_DB[X]
    else:
        print("\n ---------------- WARNING -----------------")
        print(" No molar mass for species (", X, ")")
        print(" ------------------------------------------\n")
        return 0.0


# Small test
if __name__ == "__main__":
    Mm = molarMass("U235")
    print("M(U235) =", Mm, "kg/mol")