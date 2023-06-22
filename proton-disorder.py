import numpy as np
import pandas as pd
import random

# Ice indeces for clarity
OXY_COORD_INDEX = 0
LINKS_INDEX = 1
OXY_LINKS_INDEX = 2

# Def general water constants (for TIP4P)
R_OH = 0.9572 # Angstroms
ANG_HOH = 104.52 * (np.pi/180.) # in radians
R_OM = 0.15 # Angstroms

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

def neighbour_list(ice,h,hinv):
    """
    """

    # Bond cutoff radius
    rcut = 3.
    rcut2 = rcut**2

    # Create a dictionary of links
    ice_links = {}
    oxy_bonds = []
    for i in range(len(ice[OXY_COORD_INDEX])):
        oxy_bonds.append([])

    for i in range(len(ice[OXY_COORD_INDEX])):
        # oxy_links = []
        for j in range(len(ice[OXY_COORD_INDEX])):
            # Check if the looking at the same oxygen atom
            if i == j:
                continue

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
                if i < j:
                    link_tuple = (i,j)
                else:
                    link_tuple = (j,i)
                
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

    return ice

def init_hydrogens(ice):
    """
    """
    ice_links = ice[LINKS_INDEX]
    for key in ice_links.keys():
        # Get a random oxygens links
        z = random.random()

        if (z < 0.5):
            ice_links[key] = 0
        else:
            ice_links[key] = 1

    return ice

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
    oxy_links = ice[OXY_LINKS_INDEX]

    keys = list(ice_links.keys())
    while True:
        z = random.random()
        index = int(len(keys)*z)
        if index >= len(keys):
            continue

        link = keys[index]
        oxy1_index = link[0]
        oxy2_index = link[1]

        oxy1_nbonds = count_bonds(ice, oxy1_index)
        oxy2_nbonds = count_bonds(ice, oxy2_index)

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
        if two_bonds(ice):
            break

    return ice

def two_bonds(ice):
    """
    """
    val = True
    oxy_links = ice[OXY_LINKS_INDEX]

    for oxy_index in range(len(oxy_links)):
        nbonds = count_bonds(ice, oxy_index)
        if nbonds > 2:
            val = False
            break

    return val

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

    if ice_links[link] == 0:
        ice_links[link] = 1
    else:
        ice_links[link] = 0

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
            dist = coordinates[oxy1_index] - coordinates[oxy2_index]
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

def put_hydrogens(ice,h,hinv):
    """
    """
    cosal2 = np.cos(ANG_HOH/2.)
    sinal2 = np.sin(ANG_HOH/2.)

    coordinates = ice[OXY_COORD_INDEX]
    ice_links = ice[LINKS_INDEX]
    oxy_links = ice[OXY_LINKS_INDEX]

    atoms = []
    resnum = 1
    for oxy1_index in range(len(coordinates)):
        bonds = get_bonds(ice, oxy1_index)
        poles = []
        for oxy2_index in bonds:
            dist = coordinates[oxy1_index] - coordinates[oxy2_index]
            s = np.matmul(hinv,dist)
            s = s - np.rint(s)
            dist = np.matmul(h,s)
            poles.append(dist)
        
        u = poles[0] + poles[1]
        u = u / np.sqrt(np.sum(u**2))
        v = poles[0] - poles[1]
        v = v / np.sqrt(np.sum(v**2))
        coord_h1 = coordinates[oxy1_index] + R_OH * (cosal2 * u + sinal2 * v)
        coord_h2 = coordinates[oxy1_index] + R_OH * (cosal2 * u - sinal2 * v)

        # Add oxygen
        atoms.append((resnum,coordinates[oxy1_index],'O','OW'))
        # Add hydrogens
        atoms.append((resnum,coord_h1,'H','HW1'))
        atoms.append((resnum,coord_h2,'H','HW2'))
        # Add virtual point
        atoms.append((resnum,coordinates[oxy1_index],'M','MW'))
    
    return atoms


def output(ice,h,hinv,fileformat='xyz',filename='ice'):
    """
    """
    atoms = put_hydrogens(ice,h,hinv)
    # Write unit cell pdb

    f = open("output/ice.xyz", "w")

    # f.write("CRYST1   17.310   17.310   17.310  90.00  90.00  90.00 P 1           1\n")
    f.write("{}\n\n".format(len(atoms)))
    for atom in atoms:
        # #Collect PDB variables
        # id = str(i)
        # atom = thf_line['ELEMENT'].iloc[i]
        # resname = thf_line['MOL'].iloc[i]
        # resnum = thf_line['RESIDUE_NUM'].iloc[i]
        # x = thf_line['X'].iloc[i]
        # y = thf_line['Y'].iloc[i]
        # z = thf_line['Z'].iloc[i]

        # # Create PDB entry
        # line = "{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
        #     "ATOM", i+1, atom, resname, resnum, x, y, z, 0., 0.)
        line = "{},{},{},{}\n".format(atom[2],atom[1][0],atom[1][1],atom[1][2])
        f.write(line)

    f.close()
    print("Output written")

# ice, box_dim = parse_csv('input/s2-hydrate.csv')
# h, hinv = gen_matrices(box_dim)
# ice = neighbour_list(ice,h,hinv)
# ice = init_hydrogens(ice)
# ice = shake_bonds(ice)
# ice = adjust_bonds(ice)
# dipole = get_dipole(ice,h,hinv)

if __name__ == "__main__":
    dipole_target = 0.1

    ice, box_dim = parse_csv('input/s2-hydrate.csv')
    h, hinv = gen_matrices(box_dim)
    ice = neighbour_list(ice,h,hinv)
    ice = init_hydrogens(ice)
    ice = adjust_bonds(ice)
    dipole = get_dipole(ice,h,hinv)

    nround = 0
    while True:
        ice_old = ice
        dipole_old = dipole

        ice = shake_bonds(ice)
        ice = adjust_bonds(ice)

        dipole = get_dipole(ice,h,hinv)
        if (dipole < dipole_old):
            print("Round {}, dipole is now {:.5f}".format(nround,dipole))
        else:
            ice = ice_old
            dipole = dipole_old
        
        if (dipole < dipole_target):
            break
            
        nround += 1

    output(ice,h,hinv)