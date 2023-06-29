from proton_disorder.src.constants import *
import random
import numpy as np

def shake_bonds(ice, nshakes=20):
    """
    First step in monte carlo algorithm. This function will randomly choose a link,
    then swap the hydrogen bond from one oxygen to another. This 'mixes' the hydrogen
    configuration.
    """
    # Get links dictionary
    ice_links = ice[LINKS_INDEX]
    keys = list(ice_links.keys())

    for i in range(nshakes):
        # Generate random index, get random key
        z = random.random()
        index = int(len(keys) * z)
        if index >= len(keys):
            break
        
        link = keys[i]
        swap_bond(ice, link)

    return ice

def adjust_bonds(ice):
    """
    After initializing and shaking hydrogens at random, the each oxygen will have a range
    of hydrogens between 0-4. This will correct that by adjusting at random the bonds such that
    all oxygens will have two bonds to hydrogens. 

    This is the slowest function in the entire codebase. Because of the randomness of choosing
    a link, it previously would contiue to miss links that had the incorrect number of hydrogens,
    meaning that it would take forever with large systems. I corrected this by focusing the random
    selection of links on the indeces that do not have the right number of hydrogens. This is not 
    a purely random approach, but it saves significantly on computational time and complexity.
    """
    ice_links = ice[LINKS_INDEX]
    nbonds = ice[NBONDS_INDEX]

    tm_links= ice_links.copy()
    while True:

        success, tm_links = two_bonds(ice)
        if success:
            break
        
        keys = list(tm_links.keys())
        z = random.random()
        index = int(len(keys)*z)
        if index >= len(keys):
            continue
        link = keys[index]
        oxy1_index = link[0]
        oxy2_index = link[1]

        oxy1_nbonds = nbonds[oxy1_index]
        oxy2_nbonds = nbonds[oxy2_index]

        # Get the original difference between nbonds (should be 0 if they are both 2)
        diff_old = np.abs(oxy1_nbonds - oxy2_nbonds)
        # If the bond belongs to oxy1, swap it
        # Otherwise, it belongs to oxy2
        if ice_links[link] == 0:
            oxy1_nbonds -= 1
            oxy2_nbonds += 1
        else:
            oxy2_nbonds -= 1
            oxy1_nbonds += 1
        
        diff_new = np.abs(oxy1_nbonds - oxy2_nbonds)

        if diff_new <= diff_old:
            swap_bond(ice, link)

    return ice

def two_bonds(ice):
    """
    This function iterates through the ice and confirms that each oxygen has two
    hydrogens bonded to it. If it finds an oxygen with not two hydrogens bonded to it,
    it will flag the run as a failure, and record all the indeces that still do not have
    the correct number of hydrogens.
    """
    success = True
    nbonds = ice[NBONDS_INDEX]
    ice_links = ice[LINKS_INDEX]
    tm_links = {}

    keys = list(ice_links.keys())
    for link in keys:
        if nbonds[link[0]] == 2 and nbonds[link[1]] == 2:
            continue
        else:
            success = False
            tm_links[link] = ice_links[link]

    return success, tm_links

def count_bonds(ice, oxy_index):
    """
    This function will intially count all the bonds that each oxygen has. This
    is only run once, and then the nbonds parameter is stored in the ice array and 
    kept up to date.
    """
    ice_links = ice[LINKS_INDEX]
    oxy_links = ice[OXY_LINKS_INDEX]
    links = oxy_links[oxy_index]
    
    nbonds = 0
    for link in links:
        if link[0] == oxy_index and ice_links[link] == 0:
            nbonds += 1
        elif link[1] == oxy_index and ice_links[link] == 1:
            nbonds += 1

    return nbonds

def swap_bond(ice, link):
    """
    This function swaps the hydrogen bond between two oxygens in a link, and then 
    updates the number of bonds each of those oxygens has.
    """
    ice_links = ice[LINKS_INDEX]
    nbonds = ice[NBONDS_INDEX]

    if ice_links[link] == 0:
        ice_links[link] = 1
        nbonds[link[0]] -= 1
        nbonds[link[1]] += 1
    else:
        ice_links[link] = 0
        nbonds[link[0]] += 1
        nbonds[link[1]] -= 1

def get_bonds(ice, oxy_index):
    """
    Returns the index of the bonded oxygens
    """
    ice_links = ice[LINKS_INDEX]
    oxy_links = ice[OXY_LINKS_INDEX]
    links = oxy_links[oxy_index]
    
    bonds = []
    for link in links:
        if link[0] == oxy_index and ice_links[link] == 0:
            bonds.append(link[1])
        elif link[1] == oxy_index and ice_links[link] == 1:
            bonds.append(link[0])
    return bonds

def get_dipole(ice,h,hinv):
    """
    Calculates the dipole moment of the ice based on the coordinates of the oxygen atom and
    the coordinate of the oxygen atoms with which is is bonded with and owns the hydrogen in said
    bond. To calculate the dipole moment, only the direction of the hydrogen bonds is necessary, 
    not the actual position of the hydrogens. 
    """
    coordinates = ice[OXY_COORD_INDEX]

    dipole = 0.
    for oxy1_index in range(len(coordinates)):
        bonds = get_bonds(ice, oxy1_index)
        poles = []
        for oxy2_index in bonds:
            dist = coordinates[oxy2_index] - coordinates[oxy1_index]
            s = np.matmul(hinv, dist)
            s = s - np.rint(s)
            dist = np.matmul(h,s)
            poles.append(dist)

        u = poles[0] + poles[1]
        u = u / np.sqrt(np.sum(u**2))
        dipole = dipole + u

    return np.sqrt(np.sum(dipole**2))