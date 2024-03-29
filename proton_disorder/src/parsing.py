from proton_disorder.src.constants import *
import random
import pandas as pd
import numpy as np
import argparse
import os

def parser():
    """
    Command line argument parser for __main__ function.
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
    - filename: a csv or xyz file containing atom names and xyz coordinates
    
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
    box_dim = oxy_coord.max().values - oxy_coord.min().values
    ice.append(oxy_coord.values)

    return ice, box_dim

def neighbour_list(ice,h,hinv):
    """
    Finds and records the neighboring oxygens for each oxygen atom in the ice.
    By ice rules, there should be exactly 4 neighbours for each oxygen atom. This
    will return the completed ice array.
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
        for j in range(i+1,len(ice[OXY_COORD_INDEX])):
            # Calculate the distance between the oxygens, taking into
            # account PBC
            dist = ice[OXY_COORD_INDEX][i] - ice[OXY_COORD_INDEX][j]
            s = np.matmul(hinv, dist)
            s = s - np.rint(s)
            dist_out = np.matmul(h, s)
            # Calculate if the corrected distance is within the cutoff radius
            if np.sum(dist_out**2) <= rcut2:
                link_tuple = (i,j)
                # assign to a nonesense value
                if link_tuple not in ice_links.keys():
                    ice_links[link_tuple] = -1
                    # Add tuple key to oxy bonds
                    oxy_bonds[i].append(link_tuple)
                    oxy_bonds[j].append(link_tuple)

    # Add lists to overall ice array
    ice.append(ice_links)
    ice.append(oxy_bonds)
    ice.append(nbonds)

    return ice

def init_hydrogens(ice):
    """
    Adds a hydrogen to each link present in the ice structure. It will randomly assign
    ownership of the hydrogen (or the bond) to one of the oxygens, and then updates the 
    number of bonds those oxygens have. This is only used in initialization of the ice.
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
    Generates the box matrix in a 3x3 array for later use in PBC corrections.
    Also generates the inverse matrix for the same purpose.
    """
    # Create 3x3 matrix with system box coordinates
    h = np.zeros((3,3))
    for i in range(len(box_dim)):
        h[i,i] = box_dim[i]
    hinv = np.linalg.inv(h)

    return h, hinv