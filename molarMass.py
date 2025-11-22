def molarMass(X):
    """
    Return the molar mass of a given nuclide or nucleon [kg/mol].
    Data reference: http://wwwndc.jaea.go.jp/NuC/

    Parameters
    ----------
    X : str
        Nuclide or nucleon symbol (e.g. 'U235', 'n')

    Returns
    -------
    M : float
        Molar mass [kg/mol]
    """

    # ========================== CONSTANT ==========================
    g_per_mol_to_kg_per_mol = 1e-3

    # ========================== DATABASE ==========================
    M_database = {
        'n': 1.008665 * g_per_mol_to_kg_per_mol,   # neutron
        'Th232': 232.038 * g_per_mol_to_kg_per_mol,
        'Th233': 233.040 * g_per_mol_to_kg_per_mol,
        'Pa233': 233.040 * g_per_mol_to_kg_per_mol,
        'U233': 233.039 * g_per_mol_to_kg_per_mol,
        'U235': 235.044 * g_per_mol_to_kg_per_mol,
        'U236': 236.046 * g_per_mol_to_kg_per_mol,
        'U237': 237.048 * g_per_mol_to_kg_per_mol,
        'U238': 238.051 * g_per_mol_to_kg_per_mol,
        'U239': 239.054 * g_per_mol_to_kg_per_mol,
        'Np239': 239.052 * g_per_mol_to_kg_per_mol,
        'Pu239': 239.052 * g_per_mol_to_kg_per_mol,
        'Pu240': 240.054 * g_per_mol_to_kg_per_mol,
        'Xe135': 134.905 * g_per_mol_to_kg_per_mol
    }

    # ========================== RETURN ==========================
    if X in M_database:
        M = M_database[X]
    else:
        print(f"\n ---------------- WARNING ----------------- \n No molar mass for species ({X})")
        M = 0.0

    return M
