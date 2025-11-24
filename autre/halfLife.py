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
    #  Check that the nucleon/nuclide asked, and that the associated transformation exists in the database
    if X != 'Th232' and X != 'Th233' and X != 'Pa233' and X != 'U233' and X != 'U235'and X != 'U236'and X != 'U237' and X != 'U238'\
            and X != 'U239' and X != 'Np239' and X != 'Pu239' and X != 'Pu240' and X != 'Xe135':
        print('\n WARNING : There is no database for element ', X, '. \n Please check function information')

     # Check whether the transformation exists or not
    if Transfo != 'Alpha' and Transfo != 'BetaMinus' and Transfo != 'BetaPlus' and Transfo != 'Gamma': # Transfos Gamma très courtes
        print('\n WARNING : These transformation are not implemented :', Transfo, '.\n Please check function information')

    hl = 0

    if X == 'U233' and 'Alpha':
        hl = 159.2e3 * 365 * 24 * 3600
    if X == 'U235' and 'Alpha':
        hl = 703.8e6 * 365 * 24 * 3600
    if X == 'U236' and 'Alpha':
        hl = 2.342e7 * 365 * 24 * 3600
    if X == 'U237' and 'BetaMinus':
        hl = 6.75 * 24 * 3600
    if X == 'U238' and 'Alpha':
        hl = 4.468e9 * 365 * 24 * 3600
    if X == 'U239' and 'BetaMinus':
        hl = 23.45 * 3600
    if X == 'Np239' and 'BetaMinus':
        hl = 2.356 * 24 * 3600
    if X == 'Pu239' and 'Alpha':
        hl = 24110 * 365 * 24 * 3600
    if X == 'Pu240' and 'Alpha':
        hl = 6561 * 365 * 24 * 3600
    if X == 'Pu241' and ('Alpha' or 'BetaMinus'):
        hl = 14.290 * 365 * 24 * 3600
    if X == 'Xe135' and 'BetaMinus':
        hl = 15.29 * 60
    if X == 'Th232' and 'Alpha':
        hl = 14.05e9 * 365 * 24 * 3600
    if X == 'Th233' and 'BetaMinus':
        hl = 22.3 * 3600
    if X == 'Pa233' and 'BetaMinus':
        hl = 26.975 * 24 * 3600

    return hl 
