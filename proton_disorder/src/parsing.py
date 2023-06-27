from proton_disorder.src.constants import *
import random
import pandas as pd
import numpy as np
import argparse
import os

def parser():
    """
    """
    parser = argparse.ArgumentParser(prog="Proton Disorder",
                                     description="Generates a proton disordered ice structure",
                                     epilog="nerd")
    
    parser.add_argument('-i','--input-file',type=str,help='Name of input xyz file')
    parser.add_argument('-o','--output-file',type=str,help='Name for output file')
    parser.add_argument('-d','--dipole-target',type=float,default=0.1,help='Target dipole for hydrate')
    parser.add_argument('-n','--nshakes',type=int,default=-1,help='Number of shakes to do per round')
    parser.add_argument('-t','--tip3p', action='store_true',default=False,help='Use TIP3P (Not recommended)')

    return parser

def parse_input(filename):
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

    filename_parts = os.path.splitext(filename)
    fileformat = filename_parts[len(filename_parts) - 1]
    
    ice = []
    if 'csv' in fileformat:
        coordinates=pd.read_csv(filename,header=None,names=['ATOM','X','Y','Z'],skiprows=[0,1])
    elif 'xyz' in fileformat:
        coordinates=pd.read_csv(filename,delimiter='\t',header=None,names=['ATOM','X','Y','Z'],skiprows=[0,1])
    else:
        raise Exception('Input file format "{}" not supported'.format(fileformat))
    
    oxy_coord = coordinates.query('`ATOM` == "O"')[['X','Y','Z']]

    # Coordinate array to ice, get box dimension
    # coordinates = oxy_coord.values
    box_dim = oxy_coord.max().values - oxy_coord.min().values

    # Instantiate tuple list
    # coord_tuple = []
    # for oxygen in coordinates:
    #     new_tuple = (oxygen[0], oxygen[1], oxygen[2])
    #     coord_tuple.append(new_tuple)
    ice.append(oxy_coord.values)

    return ice, box_dim

def neighbour_list(ice,h,hinv):
    """
    """

    # Bond cutoff radius
    rcut = 3.
    rcut2 = rcut**2

    # Create a dictionary of links
    ice_links = {}
    oxy_bonds = []
    nbonds = []

    for i in range(len(ice[OXY_COORD_INDEX])):
        oxy_bonds.append([])
        nbonds.append(0)

    for i in range(len(ice[OXY_COORD_INDEX])):
        # oxy_links = []
        for j in range(i+1,len(ice[OXY_COORD_INDEX])):
            # Check if the looking at the same oxygen atom
            # if i == j:
            #     continue

            # Calculate the distance between the oxygens, taking into
            # account PBC
            dist = ice[OXY_COORD_INDEX][i] - ice[OXY_COORD_INDEX][j]
            s = np.matmul(hinv, dist)
            s = s - np.rint(s)
            dist_out = np.matmul(h, s)
            # Calculate if the corrected distance is within the cutoff radius
            if np.sum(dist_out**2) <= rcut2:
                # If the current ice we are checking has less than 4 bond
                # add a bond - for debugging
                # if len(ice_links[i]) < 4:
                    # Add link to list
                    # oxy_links.append((j,False))
                # Order tuple (small,large) to test for duplications
                # if i < j:
                #     link_tuple = (i,j)
                # else:
                #     link_tuple = (j,i)

                link_tuple = (i,j)
                # assign to a nonesense value
                if link_tuple not in ice_links.keys():
                    ice_links[link_tuple] = -1
                    # Add tuple key to oxy bonds
                    oxy_bonds[i].append(link_tuple)
                    oxy_bonds[j].append(link_tuple)

                # else:
                #     print("Too many neighbours")
        
        # ice_links[i] = oxy_links

    ice.append(ice_links)
    ice.append(oxy_bonds)
    ice.append(nbonds)

    return ice

def init_hydrogens(ice):
    """
    """
    ice_links = ice[LINKS_INDEX]
    nbonds = ice[NBONDS_INDEX]
    for key in ice_links.keys():
        # Get a random oxygens links
        z = random.random()

        if (z < 0.5):
            ice_links[key] = 0
            nbonds[key[0]] += 1
        else:
            ice_links[key] = 1
            nbonds[key[1]] += 1

    return ice

def gen_matrices(box_dim):
    """
    """
    # Create 3x3 matrix with system box coordinates
    h = np.zeros((3,3))
    for i in range(len(box_dim)):
        h[i,i] = box_dim[i]
    hinv = np.linalg.inv(h)

    return h, hinv