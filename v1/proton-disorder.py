import numpy as np
import pandas as pd
from Oxygen import Oxygen
from Link import Link
import random


def parse_csv(csv_file):
    """
    Author: Gabe Miles
    Date: 6/14/23
    Description: Translates an csv with XYZ coordinates of an ice into an array
    of Oxygen objects. 

    Input
    - csv_file: an xyz file converted in csv format (comma deliminated). Should include
      columns with x, y, z coordinates, and atom name. 
    
    Output
    - ice: an array of Oxygen objects that represents the inputted structure.
    """
    
    ice = []

    # Parse XYZ file, get only oxygen atom coordinates
    coord=pd.read_csv(csv_file,index_col=0)
    oxy_coord=pd.DataFrame()
    oxy_coord[['X','Y','Z']]=coord[['OX','OY','OZ']]
    box_dim = oxy_coord.max().values - oxy_coord.min().values

    # Generate array of Oxygen objects
    for i in range(len(oxy_coord.index)):
        new_coord = oxy_coord.iloc[i].values
        oxy = Oxygen()
        oxy.set_coord(new_coord)
        ice.append(oxy)

    return ice, box_dim

def neighbour_list(ice, box_dim):
    """
    """
    rcut = 3.
    rcut2 = rcut**2

    links = []

    # Create 3x3 matrix with system box coordinates
    h = np.zeros((3,3))
    for i in range(len(box_dim)):
        h[i,i] = box_dim[i]
    hinv = np.linalg.inv(h)

    for i in range(len(ice)):
        for j in range(i+1,len(ice)):
            # Calculate the distance between the oxygens, taking into
            # account PBC
            dist = ice[i].coord - ice[j].coord
            s = np.matmul(hinv, dist)
            s = s - np.rint(s)
            dist_out = np.matmul(h, s)
            # Calculate if the corrected distance is within the cutoff radius
            if np.sum(dist_out**2) <= rcut2:
                if ice[i].nneighbours < 4 and ice[j].nneighbours < 4:
                    # Add link to list
                    new_link = Link(i,j)
                    links.append(new_link)
                    # Update the number of neighbours an oxygen has, record link index
                    ice[i].add_link(len(links) - 1)
                    ice[j].add_link(len(links) - 1)
                else:
                    print("Too many neighbours")
    
    return ice, links

def init_hydrogens(ice, links):
    """
    """
    for link in links:
        z = random.random()
        # Assign hydrogen to one oxygen of the other randomly
        # Does not matter is the oxygen already has 2 hydrogens
        if (z < 0.5):
            link.make_oxy1_bond()
            ice[link.oxy1].add_bond()
        else:
            link.make_oxy2_bond()
            ice[link.oxy2].add_bond()

    return ice, links

def shake_bonds(nshakes, ice, links):
    """
    """
    for i in range(nshakes):
        z = random.random()
        index = int(len(links) * z)
        if index >= len(links):
            break
        link = links[index]
        link.swap_bond()
        # If this oxygen just gained a bond, add a bond and remove one from the neightbour
        if link.oxy1_bond:
            ice[link.oxy1].add_bond()
            ice[link.oxy2].remove_bond()
        else:
            ice[link.oxy2].add_bond()
            ice[link.oxy1].remove_bond()

        oxy1_nbonds = ice[link.oxy1].nbonds
        oxy2_nbonds = ice[link.oxy2].nbonds
    
    return ice, links

def adjust_bonds(ice, links):
    """
    """
    while True:
        z = random.random()
        index = int(len(links)*z)
        if index >= len(links):
            continue

        link = links[index]
        oxy1_nbonds = ice[link.oxy1].nbonds
        oxy2_nbonds = ice[link.oxy2].nbonds
        # Get the original difference between nbonds (should be 0 if they are both 2)
        # However, there is a possibility to have 3 and 3 or 1 and 1 - not perfect
        diff_old = np.abs(oxy1_nbonds - oxy2_nbonds)
        if oxy1_nbonds == 2 and oxy2_nbonds == 2:
            pass
        else:
            if link.oxy1_bond:
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
            link.swap_bond()
            ice[link.oxy1].set_nbonds(oxy1_nbonds)
            ice[link.oxy2].set_nbonds(oxy2_nbonds)

        if two_bonds(ice):
            break

    return ice, links

def two_bonds(ice):
    """
    """
    val = True
    for oxygen in ice:
        if oxygen.nbonds > 2:
            val = False
            break
        
    return val

def get_dipole(ice, links):
    dipole = 0.0
    for i in range(len(ice)):
        pass

ice, box_dim = parse_csv('input/s2-hydrate.csv')
ice, links = neighbour_list(ice, box_dim)
ice, links = init_hydrogens(ice, links)
# for i in range(50):
#     print(ice[i].nbonds)
ice, links = shake_bonds(50, ice, links)
ice, links = adjust_bonds(ice, links)

print("Link example:\nOxy 1: {}\nOxy 2: {}\nOxy 1 bond: {}\nOxy 2 bond: {}".format(
    links[0].oxy1, links[0].oxy2, links[0].oxy1_bond, links[0].oxy2_bond
))

# Debug
# max_neighbours = 0
# for oxygen in ice:
#     if oxygen.nneighbours > max_neighbours:
#         max_neighbours = oxygen.nneighbours

# print(max_neighbours)
# def add_hydrogens():


