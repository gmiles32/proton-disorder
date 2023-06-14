import numpy as np
import pandas as pd
from Oxygen import Oxygen
from Link import Link


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
                    new_link = Link(i,j,True)
                    links.append(new_link)
                    # Update the number of neighbours an oxygen has
                    ice[i].add_neighbour()
                    ice[j].add_neighbour()
                else:
                    print("Too many neighbours")
    
    return ice, links
                    

ice, box_dim = parse_csv('input/s2-hydrate.csv')
ice, links = neighbour_list(ice, box_dim)

max_neighbours = 0
for oxygen in ice:
    if oxygen.nneighbours > max_neighbours:
        max_neighbours = oxygen.nneighbours

print(max_neighbours)
# def add_hydrogens():


