############################################################################### database of half-lives in seconds
import math

# half-lives database in seconds
# (only what we really need; others are approximations)
HL_DB = {
    # Thorium
    ("Th232", "Alpha"): 4.4e17,      # ~1.4e10 yr
    ("Th233", "BetaMinus"): 5.0e3,   # ~1.4 h

    # Protactinium
    ("Pa233", "BetaMinus"): 9.9e4,   # ~27.4 h

    # Uranium (mostly irrelevant on short times, but put something)
    ("U233", "Alpha"): 5.2e12,       # ~1.6e5 yr
    ("U235", "Alpha"): 2.2e16,       # ~7e8 yr
    ("U236", "Alpha"): 7.4e14,       # ~2.3e7 yr
    ("U237", "BetaMinus"): 2.2e3,    # ~38 min
    ("U238", "Alpha"): 1.4e17,       # ~4.5e9 yr
    ("U239", "BetaMinus"): 2.3e2,    # ~23 min

    # Neptunium
    ("Np239", "BetaMinus"): 8.4e4,   # ~2.35 d

    # Plutonium
    ("Pu239", "Alpha"): 7.6e15,      # ~2.4e4 yr
    ("Pu240", "Alpha"): 2.0e15,      # ~6.6e3 yr

    # Xenon-135
    ("Xe135", "BetaMinus"): 3.29e4,  # ~9.14 h

    # Fictitious fission product (delayed neutrons source)
    ("FP", "BetaMinus"): 1.0,        # given: T1/2 ~ 1 s
}


def halfLife(X, Transfo):
    """
    Return the half life of a nuclide or nucleon for a specific transformation
    Data available on http://wwwndc.jaea.go.jp/NuC/

    ----------------
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        i.e.: for Uranium-235, X = 'U235'
    :param Transfo: string
        name of the considered transformation, should be one of:
        ['Alpha', 'BetaMinus', 'BetaPlus', 'Gamma']
    :return: the half life of X for the transformation Transfo in seconds
    """

    if (X, Transfo) in HL_DB:
        hl = HL_DB[(X, Transfo)]
    else:
        print("\n WARNING : No half-life in database for (", X, ",", Transfo, ")")
        hl = 0.0

    return hl


# Small test
if __name__ == "__main__":
    hl_fp = halfLife("FP", "BetaMinus")
    hl_xe = halfLife("Xe135", "BetaMinus")
    print("T1/2(FP)   =", hl_fp, "s")
    print("T1/2(Xe135) =", hl_xe, "s")