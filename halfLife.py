import numpy as np

def halfLife(X, Transfo):
    """
    Return the half life of a nuclide or nucleon for a specific transformation
    Data available on http://wwwndc.jaea.go.jp/NuC/ (more precisely: http://wwwndc.jaea.go.jp/CN14/index.html)

    ----------------
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        i.e.: for the Uranium 235, X = 'U235'
    :param Transfo: string
        name of the considered transformation, should be one of the following in the list
        ['Alpha', 'BetaMinus', 'BetaPlus', 'Gamma']
    :return: the half life of X for the transformation Transfo in seconds
    """

    # ==================================  Check arguments  ==================================
    #  Check that the nucleon/nuclide asked exists in our database
    valid_nuclides = [
        'Th232', 'Th233', 'Pa233', 'U233', 'U235', 'U236', 'U237', 'U238',
        'U239', 'Np239', 'Pu239', 'Pu240', 'Xe135'
    ]
    if X not in valid_nuclides:
        print('\n WARNING : There is no database for element ', X, '. \n Please check function information')

    # Check whether the transformation exists or not
    valid_transfo = ['Alpha', 'BetaMinus', 'BetaPlus', 'Gamma']
    if Transfo not in valid_transfo:  # Gamma very short, usually neglected
        print('\n WARNING : These transformations are not implemented :', Transfo, '.\n Please check function information')

    # =============================  Internal half-life database  =============================
    # Values are total half-lives for the dominant decay mode (in seconds).
    # (Orders of grandeur based on nuclear data tables; sufficient for this simplified model.)

    year = 365.25 * 24 * 3600
    day  = 24 * 3600
    hour = 3600
    minute = 60

    hl_database = {
        # Thorium
        'Th232': 1.4e9 * year,          # ~1.4×10^9 years, alpha :contentReference[oaicite:0]{index=0}
        'Th233': 21.8 * minute,          # ~21.8 min, beta- to Pa-233 :contentReference[oaicite:1]{index=1}

        # Protactinium
        'Pa233': 26.97 * day,            # ~26.97 days, beta- to U-233 :contentReference[oaicite:2]{index=2}

        # Uranium
        'U233': 1.592e5 * year,          # ~1.59×10^5 years, alpha :contentReference[oaicite:3]{index=3}
        'U235': 7.04e8 * year,           # ~7.04×10^8 years, alpha :contentReference[oaicite:4]{index=4}
        'U236': 2.342e7 * year,          # ~2.34×10^7 years :contentReference[oaicite:5]{index=5}
        'U237': 6.75 * day,              # ~6.75 days :contentReference[oaicite:6]{index=6}
        'U238': 4.468e9 * year,          # ~4.47×10^9 years, alpha :contentReference[oaicite:7]{index=7}
        'U239': 23.45 * minute,          # ~23.45 min, beta- to Np-239 :contentReference[oaicite:8]{index=8}

        # Neptunium
        'Np239': 2.3565 * day,           # ~2.3565 days, beta- to Pu-239 :contentReference[oaicite:9]{index=9}

        # Plutonium
        'Pu239': 2.411e4 * year,         # ~2.41×10^4 years, alpha :contentReference[oaicite:10]{index=10}
        'Pu240': 6561 * year,            # ~6561 years, alpha :contentReference[oaicite:11]{index=11}

        # Xenon
        'Xe135': 9.14 * hour             # ~9.14 h, beta- to Cs-135 :contentReference[oaicite:12]{index=12}
    }

    # =============================  Return value  =============================
    # Ici, on ne différencie pas la demi-vie selon Alpha/BetaMinus/etc.,
    # on renvoie la demi-vie totale de l’isotope (suffisant pour ton modèle).
    if X in hl_database:
        hl = hl_database[X]
    else:
        # valeur fallback très grande (quasi stable)
        hl = np.inf

    return hl
