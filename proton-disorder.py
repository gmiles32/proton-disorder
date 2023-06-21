import numpy as np
import pandas as pd
import random

# global variables
ice = []
# indexes for ice parts
oxy = 0
link = 1

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
    
    global ice
    global oxy, link

    # Parse XYZ file, get only oxygen atom coordinates
    coord=pd.read_csv(csv_file,index_col=0)
    oxy_coord=pd.DataFrame()
    oxy_coord[['X','Y','Z']]=coord[['OX','OY','OZ']]

    # Coordinate array to ice, get box dimension
    coordinates = oxy_coord.values
    box_dim = oxy_coord.max().values - oxy_coord.min().values
    ice.append(coordinates)

    return box_dim

def neighbour_list(box_dim):
    """
    """
    global ice
    global oxy_index, links_index

    rcut = 3.
    rcut2 = rcut**2

    # Create a dictionary of links
    links = {}
    for i in range(len(ice[oxy_index])):
        links[i] = []

    # Add links to ice
    ice[links_index] = links

    # Create 3x3 matrix with system box coordinates
    h = np.zeros((3,3))
    for i in range(len(box_dim)):
        h[i,i] = box_dim[i]
    hinv = np.linalg.inv(h)

    for i in range(len(ice[oxy_index])):
        oxy_links = []
        for j in range(i+1,ice[oxy]):
            # Calculate the distance between the oxygens, taking into
            # account PBC
            dist = ice[oxy_index][i] - ice[oxy_index][j]
            s = np.matmul(hinv, dist)
            s = s - np.rint(s)
            dist_out = np.matmul(h, s)
            # Calculate if the corrected distance is within the cutoff radius
            if np.sum(dist_out**2) <= rcut2:
                if len(links[i]) < 4 and len(links[j]) < 4:
                    # Add link to list
                    new_link = {j:False}
                    oxy_links.append(new_link)
                else:
                    print("Too many neighbours")
        
        links[i] = oxy_links

    ice[links_index] = links

ice, box_dim = parse_csv('input/s2-hydrate.csv')