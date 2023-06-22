import numpy as np
import pandas as pd
import random

# Ice indeces for clarity
OXY_INDEX = 0
LINKS_INDEX = 1

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

    # Coordinate array to ice, get box dimension
    coordinates = oxy_coord.values
    box_dim = oxy_coord.max().values - oxy_coord.min().values

    # Instantiate tuple list
    # coord_tuple = []
    # for oxygen in coordinates:
    #     new_tuple = (oxygen[0], oxygen[1], oxygen[2])
    #     coord_tuple.append(new_tuple)
    ice.append(coordinates)

    return ice, box_dim

def neighbour_list(ice,box_dim):
    """
    """

    # Bond cutoff radius
    rcut = 3.
    rcut2 = rcut**2

    # Create a dictionary of links
    ice_links = []
    for i in range(len(ice[OXY_INDEX])):
        ice_links.append({})

    # Create 3x3 matrix with system box coordinates
    h = np.zeros((3,3))
    for i in range(len(box_dim)):
        h[i,i] = box_dim[i]
    hinv = np.linalg.inv(h)

    for i in range(len(ice[OXY_INDEX])):
        oxy_links = []
        for j in range(len(ice[OXY_INDEX])):
            # Check if the looking at the same oxygen atom
            if i == j:
                continue

            # Calculate the distance between the oxygens, taking into
            # account PBC
            dist = ice[OXY_INDEX][i] - ice[OXY_INDEX][j]
            s = np.matmul(hinv, dist)
            s = s - np.rint(s)
            dist_out = np.matmul(h, s)
            # Calculate if the corrected distance is within the cutoff radius
            if np.sum(dist_out**2) <= rcut2:
                # If the current ice we are checking has less than 4 bond
                # add a bond - for debugging
                if len(ice_links[i]) < 4:
                    # Add link to list
                    oxy_links.append((j,False))
                else:
                    print("Too many neighbours")
        
        ice_links[i] = oxy_links

    ice.append(ice_links)

    return ice

def init_hydrogens(ice):
    """
    """
    ice_links = ice[LINKS_INDEX]
    for i in range(2 * len(ice_links)):
        # Get a random oxygens links
        z1 = random.random()
        oxy1_index = int(z1 * 2 * len(ice_links))

        # Get that oxygens keys
        oxy1 = ice_links[oxy1_index]
        oxy1_links = list(oxy1.keys())
        
        rand_bond = random.randint(0,3)
        # Assign hydrogen to one oxygen of the other randomly
        # Does not matter is the oxygen already has 2 hydrogens
        oxy1 = ice_links[oxy1_index]
        oxy1_links = list(oxy1.keys())
        oxy2 = halflink1[0]

        link_oxygens = list(halflink1.keys())
        rand_index = random.randint(0,3)
        halflink


        halflink2 = ice_links[halflink1[0]][i]
        if (z1 < 0.5):

            link.make_oxy1_bond()
            ice[link.oxy1].add_bond()
        else:
            link.make_oxy2_bond()
            ice[link.oxy2].add_bond()

    return ice, links

ice, box_dim = parse_csv('input/s2-hydrate.csv')
ice = neighbour_list(ice, box_dim)

print(ice)
