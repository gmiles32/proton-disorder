import numpy as np
import pandas as pd
from random import randint
from scipy.spatial.transform import Rotation as R
from multiprocessing import Pool

def virtual_coord(x1, x2, x3):
    """
    Returns the virtual coordinate (or 4th point) for a TIP4P water using the following alogrithm:
    x4 = x1 + a*(x2-x1) + b*(x3-x1)
    where a,b = 0.128012065

    Refer to OPLSAA TIP4P itp file in GROMACS for more information

    Author: Gabe Miles
    """
    a = b = 0.128012065
    return x1 + a * (x2 -x1) + b * (x3 - x1)

def rotate_hydrogens(hydrate, rot_hydrate, index):
    """
    This function randomizes the orientation of THF in space. It does this with a transform matrix
    provided by the scipy Rotation package. By using this package, THF can be rotated around it's
    center in 3D space, and not just around an axis. This funtion also provides random rotations,
    with a random number of degrees in each direction.

    For more information on the theory, visit:
    https://stackoverflow.com/questions/14607640/rotating-a-vector-in-3d-space

    Input: THF dataframe that includes columns for 'X', 'Y', and 'Z'

    Output: Identical THF dataframe with updated coordinates given a random rotation in 3D space

    Author: Gabe Miles
    """
    rot_hydrate = pd.DataFrame(columns=['ELEMENT', 'X', 'Y', 'Z', 'RESIDUE_NUM', 'MOL'])
    water = hydrate.loc[hydrate['RESIDUE_NUM'] == index].copy()
    # Center oxygen on (0,0,0)
    center = water[['X','Y','Z']].loc[water['ELEMENT'] == 'OW'].values[0]
    water[['X','Y','Z']] = water[['X','Y','Z']].values - center
    
    # Notate HW and MW points together a random distance
    # This maintains the proper bond angles
    x,y,z = randint(0,360), randint(0,360), randint(0,360)
    r = R.from_euler('zyx', [x,y,z], degrees=True)
    new_xyz = (r.as_matrix() @ water[['X','Y','Z']].values.T).T
    rot_water = pd.DataFrame({'ELEMENT': water['ELEMENT'],
                    'X': new_xyz[:,0],
                    'Y': new_xyz[:,1],
                    'Z': new_xyz[:,2],
                    'RESIDUE_NUM': water['RESIDUE_NUM'],
                    'MOL': water['MOL'],
                    'CHARGE': water['CHARGE']})
                    
    # Move water center back to where it should be, add water to new hydrate
    rot_water[['X','Y','Z']] = rot_water[['X','Y','Z']].values + center
    rot_hydrate = pd.concat([rot_hydrate, rot_water], axis=0, ignore_index=True)

    # return rot_hydrate

def dipole(hydrate):
    x_dipole = (hydrate['X'] * hydrate['CHARGE']).sum()
    y_dipole = (hydrate['Y'] * hydrate['CHARGE']).sum()
    z_dipole = (hydrate['Z'] * hydrate['CHARGE']).sum()
    return x_dipole + y_dipole + z_dipole

def gen_unit():
    hydrate=pd.read_csv('input/s2-hydrate.csv',index_col=0)
    hydrate_ox=pd.DataFrame()
    hydrate_ox[['X','Y','Z']]=hydrate[['OX','OY','OZ']]
    hydrate_ox.insert(0,'ELEMENT','OW')
    hydrate_hy=pd.DataFrame()
    hydrate_hy[['X','Y','Z']]=hydrate[['H1X','H1Y','H1Z']]
    hydrate_hy.insert(0, 'ELEMENT', 'HW1')
    hydrate_hz=pd.DataFrame()
    hydrate_hz[['X','Y','Z']]=hydrate[['H2X','H2Y','H2Z']]
    hydrate_hz.insert(0, 'ELEMENT', 'HW2')

    mw_x_coord = virtual_coord(hydrate_ox['X'], hydrate_hy['X'], hydrate_hz['X'])
    mw_y_coord = virtual_coord(hydrate_ox['Y'], hydrate_hy['Y'], hydrate_hz['Y'])
    mw_z_coord = virtual_coord(hydrate_ox['Z'], hydrate_hy['Z'], hydrate_hz['Z'])

    hydrate_mw = pd.DataFrame({'X': mw_x_coord,
                            'Y': mw_y_coord,
                            'Z': mw_z_coord})
    hydrate_mw.insert(0,'ELEMENT','MW')
    unit = pd.DataFrame(columns=['ELEMENT', 'X', 'Y', 'Z', 'RESIDUE_NUM'])

    residue_num = 1
    for row in range(len(hydrate_ox.index)):
        water = pd.concat([hydrate_ox.iloc[row], hydrate_hy.iloc[row], hydrate_hz.iloc[row], hydrate_mw.iloc[row]],
                        axis=1, ignore_index=True).T   
        water['RESIDUE_NUM'] = residue_num
        unit = pd.concat([unit, water], axis=0, ignore_index=True)
        residue_num += 1

    # Add mol and charge parameters
    unit['MOL'] = 'HOH'
    unit['CHARGE'] = 0.
    unit['CHARGE'].loc[unit['ELEMENT'] == 'HW1'] = 0.52 
    unit['CHARGE'].loc[unit['ELEMENT'] == 'HW2'] = 0.52  
    unit['CHARGE'].loc[unit['ELEMENT'] == 'MW'] = -1.04 

    return unit

def gen_hydrate():
    cell_size = 17.30
    num_cells = 3
    coord_offset = [[-cell_size, 0, cell_size],[0,0,cell_size],[cell_size, 0, cell_size],
                    [-cell_size, 0, 0],[cell_size, 0, 0],
                    [-cell_size, 0, -cell_size],[0, 0, -cell_size],[cell_size, 0, -cell_size]]
    
    unit = gen_unit()

    # Create 3x3 plane crystal
    hydrate_3_plane = unit.copy()
    for offset in coord_offset:
        # For hydrate
        residue_num = hydrate_3_plane['RESIDUE_NUM'].iloc[len(hydrate_3_plane.index) - 1]
        hydrate_unit = unit.copy()
        hydrate_unit[['X','Y','Z']] = hydrate_unit[['X','Y','Z']] + offset
        hydrate_unit['RESIDUE_NUM'] = hydrate_unit['RESIDUE_NUM'] + residue_num
        hydrate_3_plane = pd.concat([hydrate_3_plane, hydrate_unit], axis=0, ignore_index=True)

    # Stack 3x3x3 crystal
    hydrate_3_3 = hydrate_3_plane.copy()
    for i in range(1, num_cells):
        residue_num = hydrate_3_3['RESIDUE_NUM'].iloc[len(hydrate_3_3.index) - 1]
        hydrate_unit = hydrate_3_plane.copy()
        hydrate_unit['Y'] = hydrate_unit['Y'] + (cell_size * i)
        hydrate_unit['RESIDUE_NUM'] = hydrate_unit['RESIDUE_NUM'] + residue_num
        hydrate_3_3 = pd.concat([hydrate_3_3, hydrate_unit], axis=0, ignore_index=True)
    
    return hydrate_3_3

if __name__ == '__main__':
    hydrate = gen_hydrate()
    rot_hydrate = pd.DataFrame(['ELEMENT','X','Y','Z','RESIDUE_NUM','MOL',
                    'CHARGE'])
    while True:
        if np.abs(dipole(rot_hydrate)) >= 0.05:
            #rotate hydrogens
            with Pool(processes=3) as pool:
                p = pool.map(rotate_hydrogens, args=(hydrate, rot_hydrate, range(1,len(hydrate.index))))
            continue
        else:
            break
    
