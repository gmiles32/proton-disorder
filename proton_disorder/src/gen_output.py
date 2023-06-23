from src.constants import *
from src.mc import get_bonds
import numpy as np
import os

def virtual_coord(x1, x2, x3):
    """
    Returns the virtual coordinate (or 4th point) for a TIP4P water using the following alogrithm:
    x4 = x1 + a*(x2-x1) + b*(x3-x1)
    where a,b = 0.128012065

    Refer to OPLSAA TIP4P itp file in GROMACS for more information
    """
    a = b = 0.128012065
    return x1 + a * (x2 -x1) + b * (x3 - x1)

def put_hydrogens(ice,h,hinv,tip4p):
    """
    Generates the coordinates for hydrogens in ice structure. Additionally will
    generate virtual point for TIP4P water. 

    Input
    - ice: ice array
    - h: 3x3 matrix with system dimensions
    - hinv: inverse matrix of h
    - tip4p: flag to generate virtual point for tip4p water. Defaults to True
             (You should be using TIP4P water dunce)

    Output: atoms array of 4-tuples ordered as (residue #, [x,y,z], atom name, pdb name)
    """

    coordinates = ice[OXY_COORD_INDEX]
    atoms = []
    resnum = 1

    for oxy1_index in range(len(coordinates)):
        bonds = get_bonds(ice, oxy1_index)
        h_coords = []
        for oxy2_index in bonds:
            # PBC correction
            r = coordinates[oxy2_index] - coordinates[oxy1_index]
            s = np.matmul(hinv,r)
            s = s - np.rint(s)
            r = np.matmul(h,s)
            # Generate coordinate for hydrogen atom
            norm = np.linalg.norm(r)
            unit_r = r / norm
            h_coord = coordinates[oxy1_index] + R_OH * unit_r
            h_coords.append(h_coord)

        # Atom tuples are ordered as (residue #, [x,y,z], atom name, pdb name)
        # Add oxygen
        atoms.append((resnum,coordinates[oxy1_index],'O','OW'))
        # Add hydrogens
        atoms.append((resnum,h_coords[0],'H','HW1'))
        atoms.append((resnum,h_coords[1],'H','HW2'))
        
        if tip4p:
            # Add virtual point
            atoms.append((resnum,virtual_coord(coordinates[oxy1_index],h_coords[0],h_coords[1]),'M','MW'))

        resnum += 1

    return atoms

def output(ice,h,hinv,filename='output/ice.xyz',tip4p=True):
    """
    Generates properly formatted ice structure from ice array.

    Input
    - ice: Ice array
    - tip4p: Flag to generate virtual point for tip4p water. Defaults to True
             (You should be using TIP4P you dunce)
    - filename: Filename for output. Must have file extension pdb, xyz, or csv

    Output: Formatted output file in specified directory
    """

    atoms = put_hydrogens(ice,h,hinv,tip4p)

    filename_parts = os.path.splitext(filename)
    fileformat = filename_parts[len(filename_parts) - 1]

    f = open(filename, "w")

    if 'pdb' in fileformat:
        f.write("CRYST1   {:2.3f}   {:2.3f}   {:2.3f}  90.00  90.00  90.00 P 1           1\n".format(
            h[0][0],h[1][1],h[2][2]))
        for i in range(len(atoms)):
            atom = atoms[i]
            id = str(i)
            element = atom[PDB_NAME_INDEX]
            resname = 'HOH'
            resnum = atom[RESNUM_INDEX]
            x = atom[ATOM_COORD_INDEX][0]
            y = atom[ATOM_COORD_INDEX][1]
            z = atom[ATOM_COORD_INDEX][2]

            # Create PDB entry
            line = "{:6s}{:5d} {:^4s} {:3s}  {:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
                "ATOM", i+1, element, resname, resnum, x, y, z, 0., 0.)
            f.write(line)

    elif 'xyz' in fileformat:
        f.write("{}\n\n".format(len(atoms)))
        for atom in atoms:
            line = "{}\t{}\t{}\t{}\n".format(atom[ATOM_NAME_INDEX],
                                             atom[ATOM_COORD_INDEX][0],
                                             atom[ATOM_COORD_INDEX][1],
                                             atom[ATOM_COORD_INDEX][2])
            f.write(line)

    elif 'csv' in fileformat:
        f.write("{}\n\n".format(len(atoms)))
        for atom in atoms:
            line = "{},{},{},{}\n".format(atom[ATOM_NAME_INDEX],
                                          atom[ATOM_COORD_INDEX][0],
                                          atom[ATOM_COORD_INDEX][1],
                                          atom[ATOM_COORD_INDEX][2])
        f.write(line)
    else:
        raise Exception("Output file format '{}' not supported".format(fileformat))

    f.close()
    print("Output file " + filename + " written")
