from proton_disorder.src.constants import *
import random
import numpy as np

def shake_bonds(ice, nshakes=20):
    """
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
        # However, there is a possibility to have 3 and 3 or 1 and 1 - not perfect
        diff_old = np.abs(oxy1_nbonds - oxy2_nbonds)
        # if oxy1_nbonds == 2 and oxy2_nbonds == 2:
        #     pass
        # else:
        # If the bond belongs to oxy1, swap it
        # Otherwise, it belongs to oxy2
        if ice_links[link] == 0:
            oxy1_nbonds -= 1
            oxy2_nbonds += 1
        else:
            oxy2_nbonds -= 1
            oxy1_nbonds += 1
            # link.swap_bond()
            # ice[link.oxy1].set_nbonds(oxy1_nbonds)
            # ice[link.oxy2].set_nbonds(oxy2_nbonds)
        
        diff_new = np.abs(oxy1_nbonds - oxy2_nbonds)

        if diff_new <= diff_old:
            swap_bond(ice, link)

        # Check if all oxygens have two bonds
        # if two_bonds(ice):
        #     break
    return ice

def two_bonds(ice):
    """
    """
    success = True
    nbonds = ice[NBONDS_INDEX]
    ice_links = ice[LINKS_INDEX]
    oxy_links = ice[OXY_LINKS_INDEX]
    # new_tm_indeces = []
    tm_links = {}

    # for i in tm_indeces:
    #     n = nbonds[i]
    #     if n != 2:
    #         success = False
    #         new_tm_indeces.append(i)
    #         for link in oxy_links[i]:
    #             new_tm_links[link] = ice_links[link]

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

def gen_matrices(box_dim):
    """
    """
    # Create 3x3 matrix with system box coordinates
    h = np.zeros((3,3))
    for i in range(len(box_dim)):
        h[i,i] = box_dim[i]
    hinv = np.linalg.inv(h)

    return h, hinv
