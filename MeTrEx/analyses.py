"""This is the analyses module. 

The analyses module contains different methods for reading files and
computing analyses.

Methods:
------
distanceAtoms():
    Calculate the distance of selected atoms in the given universe.
speedMolecule():
    Calculate the speed of a selected molecule/atom in the given
    universe.
hydrogenBonds():
    Not implemented.
readXVGfile():
    Read an xvg file and return relevant information for a 
    visualization.
"""

__version__ = '1.0'
__author__ = 'Christiane Rohse'

import math
import numpy as np
#import MDAnalysis as mda
#from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

def distanceAtoms(universe, atom1, atom2): 
    """Calculate the euclidean distance between two atoms within a 
    given frame range and return the distances.

    Parameter:
    ------
    universe (Universe):
        MDAnalysis universe.
    atom1 (str):
        Selection string for the univers' first atom.
    atom2 (str):
        Selection string for the univers' second atom.

    Return:
    -----
    distance (array):
        Float values of distances for all frames of the universe, 
        in order of appearance.
    """
    atom_group1 = universe.select_atoms(atom1)
    atom_group2 = universe.select_atoms(atom2)
    pos1 = []
    pos2 = []
    # Position extraction
    for frame in universe.trajectory:
        pos1.append(atom_group1.positions[0])
        pos2.append(atom_group2.positions[0])
        # if you want to use molecules as input, you need to change the above 2 lines to:
        # pos1.append(atom_group1.center_of_mass())
        # pos2.append(atom_group2.center_of_mass())
    # Euclidean distance calculation
    distance = []
    i=0
    while i < len(pos1):
        dis = math.sqrt((pos1[i][0]-pos2[i][0])**2 + (pos1[i][1]-pos2[i][1])**2 + (pos1[i][2]-pos2[i][2])**2)
        distance.append(dis)
        i += 1
    return np.array(distance)

def speedMolecule(universe, molecule):
    """Calculate the speed (nm/ns) of a molecule by its center of mass.
    
    Parameter:
    ------
    universe (Universe):
        Reference to the Universe.
    molecule (str):
        Selection string of the molecule to be selected via MDAnalysis.
    
    Return:
    ------
    array ([speed, timesteps]):
        Calculated speed values (has number of frames-1 entries) in 
        nm/ns and corresponding time steps from the universe).
    """
    speed = []
    time_steps = [] # get time interval of frames from universe
    positions = []
    molecule_reference = universe.select_atoms(molecule)
    for ts in universe.trajectory:
        positions.append(molecule_reference.center_of_mass())
        time_steps.append(universe.trajectory.time)
    j = len(time_steps)-1
    last_time_step = time_steps[j] + (time_steps[j] - time_steps[j-1])
    time_steps.append(last_time_step)
    i=0
    imax = len(positions)-1
    while i < imax:
        dis = math.sqrt((positions[i][0]-positions[i+1][0])**2 + (positions[i][1]-positions[i+1][1])**2 + (positions[i][2]-positions[i+1][2])**2)
        time_interval = time_steps[i+1] - time_steps[i]
        speed.append(float(dis/time_interval)*1000) # change from ns to us (ns, us, ms) or: change from A to nm (* 1000)
        i += 1
    return [np.array(speed), np.array(time_steps[:-2])]

def hydrogenBonds(universe, molecule):
    """Calculate the number of hydrogenbonds for a given molecule at
    every time of the simulation.
    """
    hbonds = []
    # TODO: calculate number of hbonds
    return np.array(hbonds)

def readXVGfile(xvg):
    """Read an xvg file and return relevant information for the 
    visualisation.

    Parameter:
    -----

    Return:
    Tuple(title, x_label, y_label, columns, ind, secondary_labels):
        title (str): Title of the data
        x_label (str): Label of the x axis.
        y_label (str): Label of the y axis.
        columns (list of lists): List containing a list for each columns 
            values as float numbers.
        ind (list): List of consecutive integers, starting with 0,
            ending with maximum index of data.
        secondary_label (list of str): List containing the additional
            names of the columns. 
    """
    title = None
    x_label = 'x axis'
    y_label = 'y axis'
    columns = None
    ind = []
    secondary_labels = []
    i = 0
    f = open(xvg, 'r')
    for line in f.readlines():
        if line.startswith('#'):
            pass
        elif line.startswith('@'):
            line_strings = line.split()
            if line_strings[1] == 'title':
                title = line.split('"')[1]
            elif line_strings[1] == 'xaxis':
                x_label = line.split('"')[1]
            elif line_strings[1] == 'yaxis':
                y_label = line.split('"')[1]
            elif line_strings[1][0] == 's':
                    label = line.split('"')[1]
                    secondary_labels.append(label)
        else:
            ind.append(i)
            i += 1
            line_split = line.split()
            line_length = len(line_split)
            if columns is None: 
                columns = []
                for i in range(line_length):
                    columns.append([])
            for j in range(line_length):
                columns[j].append(float(line_split[j]))
    f.close()
    return (title, x_label, y_label, columns, ind, secondary_labels)