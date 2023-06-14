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

    h = np.zeros((3,3))
    for i in range(len(box_dim)):
        h[i,i] = box_dim[i]

    hinv = np.linalg.inv(h)

    for i in range(len(ice)):
        for j in range(i+1,len(ice)):
            dist = ice[i].coord - ice[j].coord



ice, box_dim = parse_csv('input/s2-hydrate.csv')
neighbour_list(ice, box_dim)
# def add_hydrogens():


