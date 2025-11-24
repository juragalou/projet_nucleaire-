import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def crossSection(X, Transfo, E_neutron):
    """
    Computes the cross section for a specific energy level
    ENDF database: https://www-nds.iaea.org/exfor/endf.htm

    ----------------
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        e.g.: for the Uranium 235, X = 'U235'
    :param Transfo: string
        name of the considered transformation, should be one of the following in the list
        ['Fission', 'Capture']
    :param E_neutron: array, list
        value(s) of the cross section for an incident neutron of energy level equal to E_neutron in [eV]
        the value should be in the range [1e-5; 2e7] [eV]
    :return: array, list
        values of the cross section in [barn] corresponding to the input parameters
    """
    
    # ==================================  Check arguments  ==================================
    #  Check that the nucleon/nuclide asked, and that the associated transformation exists in the database
    if X != 'Th232' and X != 'Th233' and X != 'Pa233' and X != 'U233' and X != 'U235' and X != 'U236' and X != 'U237' and X != 'U238' \
            and X != 'U239' and X != 'Np239' and X != 'Pu239' and X != 'Pu240' and X != 'Pu241' and X != 'Xe135':
        print('\n WARNING : There is no database for element ', X, '. \n Please check function information')
        return
    
     
    df = pd.read_csv('Database/'+X+'_'+Transfo+'.txt', delim_whitespace=True, names=["Energy", "Cross", "Interpolation"])
    cross_section = []
    for x in E_neutron:
        resultat = df[df["Energy"] == x]
        if not resultat.empty:
            valeur_associee = resultat['Cross'].values[0]
            cross_section.append(valeur_associee)
        else:
            for i in range(1, df["Energy"].size):
                if x > df["Energy"][i-1] and x < df["Energy"][i]:
                    pente = (df["Cross"][i]-df["Cross"][i-1])/(df["Energy"][i]-df["Energy"][i-1])
                    chatte = pente * x + df["Cross"][i] - pente*df["Energy"][i]
                    cross_section.append(chatte)
    
    return cross_section

file_path_1 = "Database/U235_Fission.txt"  # Replace with your actual file path

# Read the file and extract the columns
x_values = []
y_values = []

with open(file_path_1, "r") as file:
    for line in file:
        if line.strip():  # Skip empty lines
            x, y = map(float, line.split()[:2])  # Assuming space-separated values
            x_values.append(x)
            y_values.append(y)

# Plot the data
plt.loglog(x_values, y_values)
plt.xlabel("Energy [eV]")
plt.ylabel("Cross section [barn]")
plt.grid()
plt.show()



