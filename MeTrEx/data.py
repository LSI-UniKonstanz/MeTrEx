"""This is the data module.

The data module contains different data related classes which can be
used to store molecular simulation data.

Classes:
------
Data : 
    Data object holding all relevant information to be shown in the
    main view. 
DataSub : 
    Collection of SingleDataSets. To be shown in one view.
SingleDataSet : 
    Object to store a single data set (2D-4D).
"""

__version__ = '1.0'
__author__ = 'Christiane Rohse'

from copy import deepcopy
from threading import Lock
import os
import numpy as np
import mdtraj as mdt
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder
import surface
import colormap as colm
from colormap import returnColors

from dialogs import SelectRepresentationAtom

class Data():
    """Data object holds information about the simulation data to be explored.

    For readability and runtime error reason: all possibly accessible variables are initialized with object creation 
    (this means: a variable exists, but its value might be either nonsense or 'None', where 'None' is preferred)

    Methods:
    -----
    distinctMol():
        Collect distinct molecules from the PDB file.
    extremesPos():
        Calculate extreme values for x, y, z positions in the PDB frame.
    reduceFrames():
        Data preprocessing. 
    membranePositions():
        Collect the positions of the phosphor atoms of the 
        lipid bilayers.
    changeColormap():
        Assigns new colors to the drug molecules.
    generateView():
        Setup all variables needed for the main view.
    changeSettings():
        Recalculate lipid abstraction.
    resetSettings():
        Reset positions to originally loaded positions.
    selectFrames():
        Store a given range of frame in a file.
    calcBox():
        Calculate the model box.
    calcLipidSurfaces():
        Calculate the surface abstraction.
    """

    def __init__(self):
        self.top = None
        self.data = None
        self.reduced_data = None # reduced simulation data, optional
        self.universe = mda.Universe.empty(0) # could be None instead
        self.k = 1 # every k-th frame will be selected
        self.n = 0 # number of frames to skip at the beginning for the visualization
        self.original_frame_number = 0 # NO index! needed when you reduce the number of frames at the preprocessing step
        self.original_max_frame = 0 # index, loaded and shown number of frames, no selection
        self.box = [0,0,0]
        self.polynomial = 5
        self.membrane_extension = 15
        self.color_list = []
        self.colormap = 'plasma'

        self.min_frame = 0 # index
        self.max_frame = 0 # index
        self.timesteps_mainview = None
        self.timesteps_offset = 0
        self.unit_distance = 'nm' # to be set at the beginning?
        self.unit_time = 'ns' # to be set at the beginning?
        self.n_atoms = 0
        self.atoms = [] # maybe use comma as separator?: 'moleculetype,residuenumber,atom'
        self.atoms_positions = []
        self.atoms_dict = {}

        # self.distinctMol() -> is called upon opening files
        self.distinct_molecule_types = set() # = name of molecules
        self.distinct_molecules = set() # = name + resid number
        self.anker_selection = False # selection of representative anker atom for molecules
        self.anker = None # anker atom for representation

        self.line_representation_molecule_types = set() # this is set as well upon opening a file, redundant, check if neccessary in the end, = name of molecules
        self.plane_representation_molecule_types = set() # = name of lipid molecules

        self.drug_molecules = {} # molecules identified by name and resid
        self.drug_molecules_reference = {} # atom groups 
        ## change settings
        self.drug_molecules_positions = {} # dictionary containing the positions of the drug molecules throughout the data: x, y, z
        self.original_drug_moleucles_positions = None
        ##
        self.drug_molecules_colors = {}
        self.drug_molecules_show = {} # mol: bool
        # Maping:
        self.drug_molecules_show_mapping = False # if true: show calculated data (line collection)
        self.drug_molecules_mseg = None # segments to be shown
        self.drug_molecules_mseg_original = None
        self.drug_molecules_linecol = None # line collection
        self.drug_molecules_mcmap = 'viridis_r'
        self.drug_molecules_mvalues = None # property to be shown
        self.drug_molecules_moffset = {} # offset for the different spheres: {'PMB393' : 220, 'name': int, ...}
        self.drug_molecules_morder = [] # order of moecules in the line collection
        self.drug_molecules_mlabel = '' 
        self.drug_molecules_text = {}
        self.mextreme = [] # 

        self.lipid_molecules = set() # molecules identified by name and resid 
        self.lipid_molecules_reference = {} # atom group for leaflet0 and leaflet1
        self.lipid_molecules_show = [True, True] # index 0 = leaflet0, index 1 = leaflet1
        ## change settings
        self.lipid_molecules_positions = {} # positions for leaflet0 and leaflet1 for the whole simulation time
        self.original_lipid_molecules_positions = None
        self.lipid_molecules_surface = {}
        self.original_lipid_molecules_surface = None
        ##
        self.lipid_molecules_show_abstraction = True # use this to decide if abstraction of original positions are shown
        self.leaflet_colors = {'leaflet0':[0.69114, 0.87318, 0.4097], 'leaflet1':[0.54228, 0.74134, 0.44047]} 

        self.label_x = 'x position'
        self.label_y = 'y position'
        self.label_z = 'z position'

    def distinctMol(self):
        """Store all unique molecules, all unique molecule types and
        all atoms with their coordinates in dictionaries.

        Parameters: ?
        ----
        pdb file path

        The gro file is used as follows:
        First two and last line are skipped. All other lines include 
        information about atoms, molecules, etc. 
        molecule name: line[5:10], residue name in gro file
        molecule number: line[0:5], residue number in gro file

        The pdb file is used as follows: 
        Lines starting with 'ATOM' are selected and split to set the 
        molecule types and distinct molecules.
        molecule name: Column 4 = position [17:20] is the residue name
        in pdb file
        molecule number: Column 5 = position [22:26] is the residue 
        sequence number in pdb file
        """
        if self.top[-3:]=='gro':
            f = open(self.top, 'r')
            for line in f.readlines()[2:-1]:
                molecule_type = line[5:10].strip()
                self.distinct_molecule_types.add(molecule_type)
                molecule = molecule_type + line[0:5].strip()
                self.distinct_molecules.add(molecule)
                a = line[10:15].strip()
                atom = molecule + a #line[10:15].strip()
                self.atoms.append(atom)
                x = float(line[20:28].strip())
                y = float(line[28:36].strip())
                z = float(line[36:44].strip())
                self.atoms_positions.append([x, y, z])
                try:
                    self.atoms_dict[molecule].append(a)
                except:
                    self.atoms_dict[molecule] = []
                    self.atoms_dict[molecule].append(a)
            f.close()
        if self.top[-3:]=='pdb':
            f = open(self.top, 'r')
            for line in f.readlines():
                if line.startswith('ATOM'):
                    molecule_types = line[17:20].strip()
                    self.distinct_molecule_types.add(molecule_types)
                    molecule = molecule_types + line[22:26].strip()
                    self.distinct_molecules.add(molecule)
                    a = line[12:16].strip()
                    atom = molecule + a #line[12:16].strip()
                    self.atoms.append(atom)
                    x = float(line[29:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    self.atoms_positions.append([x, y, z])
                    try:
                        self.atoms_dict[molecule].append(a)
                    except:
                        self.atoms_dict[molecule] = []
                        self.atoms_dict[molecule].append(a)
            f.close()

        for key, value in self.atoms_dict.items():
            value.sort()
            value.insert(0, 'all')
    
    def extremePos(self):
        """"Calculate the extreme positions of the pdb frame for a plot axis setting."""
        x = [] 
        y = [] 
        z = []
        for position in self.atoms_positions:
            x.append(position[0])
            y.append(position[1])
            z.append(position[2])
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        self.mextreme = [x.min(), x.max(), y.min(), y.max(), z.min(), z.max()]

    def reduceFrames(self):
        """Reduce the simulation data file by selecting every k-th 
        frame and skipping the first n frames. The reduced data file
        is saved as 'reducedData.xtc' file for further usage with the
        MDAnalysis library. 

        Parameter: ?
        ----
        param: Data object with topology file and data file
        param: number k to select k-th frame
        param: number n to skip n frames first

        Return:
        ----
        void
        Reduced data file is stored in objects self.reduced_data variable
        """
        if self.k > 1 or self.n > 0:
            both = False
            onlyk = False
            onlyn = False
            if self.k > 1 and self.n > 0:
                both = True
            elif self.k == 1 and self.n > 0:
                onlyn = True
            else: 
                onlyk = True

            reduced_traj = []
            j=1
            t = 1
            if onlyk or both:
                for chunk in mdt.iterload(self.data, chunk=100, top=self.top, stride=self.k):
                    reduced_traj.append(chunk)
                if both: 
                    # write intermediate result to file, because single
                    # frames cannot be accessed in reduced_traj
                    intermediate = reduced_traj[0]
                    while t < len(reduced_traj):
                        intermediate += reduced_traj[t]
                        t += 1
                    intermediate.save_xtc('intermediateData.xtc')
                    del intermediate
                    self.data = 'intermediateData.xtc'
                    onlyn = True
                    reduced_traj = []
                else: 
                    new_traj = reduced_traj[0]
            if onlyn:
                for chunk in mdt.iterload(self.data, chunk=100, top=self.top, skip=self.n):
                    reduced_traj.append(chunk)
                new_traj = reduced_traj[0]

            while j < len(reduced_traj):
                new_traj += reduced_traj[j]
                j += 1
            # Save file because cannot be read by MDAnalysis otherwise
            new_traj.save_xtc('reducedData.xtc') 
            self.reduced_data = 'reducedData.xtc'
            # remove all memory consuming objects no longer needed
            del reduced_traj
            del new_traj
            try:
                os.remove('intermediateData.xtc')
            except:
                # there is no file 'intermediateData.xtc'
                pass

    def membranePositions(self):
        """Calculate membrane positions for all frames and store them
        in the lipid_molecules_positions dicitonary.
        """
        # Only needed for recalculations for the selection of different frames
        self.lipid_molecules_positions['leaflet0'] = []
        self.lipid_molecules_positions['leaflet1'] = []
        phos0 = []
        phos1 = []
        for ts in self.universe.trajectory: 
            phos0.append(self.lipid_molecules_reference['leaflet0'].positions)
            phos1.append(self.lipid_molecules_reference['leaflet1'].positions)
            leaflet0 = np.array(phos0[0])
            leaflet1 = np.array(phos1[0])
            self.lipid_molecules_positions['leaflet0'].append(leaflet0)
            self.lipid_molecules_positions['leaflet1'].append(leaflet1)
            phos0 = []
            phos1 = []
        self.lipid_molecules_positions['leaflet0'] = np.array(self.lipid_molecules_positions['leaflet0'])
        self.lipid_molecules_positions['leaflet1'] = np.array(self.lipid_molecules_positions['leaflet1'])

    def changeColormap(self, colormap):
        """Change the drug_molecules_colors dictionary and 
        the colormap string according to given colormap.

        Parameter:
        ----
        colormap (string): 
            Which colormap to choose from. This can be any Matplotlip
            supported colormap or an equivalent colormap, e.g. from
            Crameri.

        Return:
        ----
        void
        """
        self.colormap = colormap
        n = len(self.drug_molecules_reference)
        if n > 255:
            self.color_list = returnColors(255, colormap)
        else: 
            self.color_list = returnColors(n, colormap)
        i=0
        for key, value in self.drug_molecules_reference.items():
                self.drug_molecules_colors[key] = self.color_list[i]
                i += 1
                if i > len(self.color_list)-1:
                    i = 0
        # if greyscale colormap selected: change also membrane colors to grey, else keep oroginal
        if colormap in colm.color_map_grey:
            self.leaflet_colors = {'leaflet0':[0.011219, 0.011219, 0.011219], 'leaflet1':[0.011219, 0.011219, 0.011219]}
        else:
            self.leaflet_colors = {'leaflet0':[0.69114, 0.87318, 0.4097], 'leaflet1':[0.54228, 0.74134, 0.44047]}
        # set new colors to be selected from when changing a molecules color
        self.color_list = returnColors(30, colormap)

    def generateView(self):
        """Calculate all data needed to show the main view. This 
        includes the MDAnalysis Universe, the molecule and membranes
        references, positions and colors as well as the abstraction of
        the membrane. Also, the maximum frame number is set.

        Parameter:
        ----
        -
        Return:
        ----
        void
        """
        # SELECT MOLECULES
        print('step 1: create universe')
        if self.reduced_data is not None:
            self.universe = mda.Universe(self.top, self.reduced_data)
        else:
            self.universe = mda.Universe(self.top, self.data)
        print('step 2: select molecules, add color')
        self.n_atoms = int(str(self.universe).split()[2])

        for molecule_type in self.line_representation_molecule_types:
            for molecule in self.distinct_molecules:
                if molecule.startswith(molecule_type):
                    self.drug_molecules[molecule] = [molecule_type, molecule[len(molecule_type):]]
###
        if self.anker_selection == True:
            atoms_dict = dict()
            for molecule, value in self.drug_molecules.items():
                atoms_dict[molecule] = self.universe.select_atoms(f'resname {value[0]} and resid {value[1]} and not name H*')
            atom_names = set([atom.name for atom in atoms_dict[molecule] for molecule in atoms_dict])
            
            dlg = SelectRepresentationAtom(atom_names)
            dlg.exec()
            
            if dlg.selected_atom is None:
                self.anker = 'C'
            else:
                self.anker = dlg.selected_atom
        else:
            self.anker = 'C'
        print(f'\tselected anker atom: {self.anker}')
        
        for molecule, value in self.drug_molecules.items():
            try:
                self.drug_molecules_reference[molecule] = self.universe.select_atoms(f'resname {value[0]} and resid {value[1]} and name {self.anker}')
            except:
                b = False
                try:
                    self.drug_molecules_reference[molecule] = self.universe.select_atoms('resname ' + value[0] + ' and resid ' + value[1] + ' and name C')
                    print(f'\tautomated anker atom: C for {molecule}')
                except:
                    b = True
                if b:
                    try:
                        self.drug_molecules_reference[molecule] = self.universe.select_atoms('resname ' + value[0] + ' and resid ' + value[1] + ' and name CA')
                        print(f'\tautomated anker atom: CA for {molecule}')
                    except:
                        # neither C nor CA were atom names within the line representation moleucles, no lines shown
                        print('\t no anker atom for {molecule}')
                        pass
    ###
#        for molecule, value in self.drug_molecules.items():
#            # names for the atoms shown in the line representation are 'C' or 'CA'
#            b = False
#            try:
#                self.drug_molecules_reference[molecule] = self.universe.select_atoms('resname ' + value[0] + ' and resid ' + value[1] + ' and name C')
#            except: 
#                b = True
#            if b:
#                try:
#                    self.drug_molecules_reference[molecule] = self.universe.select_atoms('resname ' + value[0] + ' and resid ' + value[1] + ' and name CA')
#                except:
#                    # neither C nor CA were atom names within the line representation moleucles, no lines shown
#                    pass

        self.drug_molecules_reference = dict(sorted(self.drug_molecules_reference.items()))
        
        n = len(self.drug_molecules_reference)
        if n > 255:
            self.color_list = returnColors(255, self.colormap)
        else: 
            self.color_list = returnColors(n, self.colormap)
        i=0
        for key, value in sorted(self.drug_molecules_reference.items()):
                self.drug_molecules_positions[key] = []
                self.drug_molecules_show[key] = True # fill dictionary with bool which indicates if this molecule shall be shown in the main view
                self.drug_molecules_colors[key] = self.color_list[i]
                i += 1
                if i > len(self.color_list)-1:
                    i = 0
        self.color_list = returnColors(30, self.colormap) 

        # setup for lipid molecules
        self.plane_representation_molecule_types = self.distinct_molecule_types - self.line_representation_molecule_types
        atomGroupLipids = mda.AtomGroup([],self.universe)
        for moleculeType in self.plane_representation_molecule_types:
            atomGroupLipids += self.universe.select_atoms('resname '+moleculeType+' and name P*') # phosphor atom for surface representation
        # leaflet finder finds only two leaflets (bilayer)

        print('step 3: search for leaflets')
        if bool(atomGroupLipids):
            L = LeafletFinder(self.universe, atomGroupLipids)
            leaflet0 = L.groups(0)
            leaflet1 = L.groups(1)
            self.lipid_molecules_reference['leaflet0'] = leaflet0
            self.lipid_molecules_reference['leaflet1'] = leaflet1
            self.lipid_molecules_positions['leaflet0'] = []
            self.lipid_molecules_positions['leaflet1'] = []
            phos0 = []
            phos1 = []
        else:
            self.lipid_molecules_show = [False, False]
        print('step 4: collect positions')
        timesteps = []
        # select data from all frames
        #m = 0
        for ts in self.universe.trajectory:
            #print('frame:', m)
            #m += 1
            for key, value in self.drug_molecules_reference.items(): 
                #print('key: ', key) # name
                #print('value: ', value) # x, y, z arrays
                #print(value.positions)
                try: 
                    self.drug_molecules_positions[key].append(value.positions[0])
                except: # molecule has no atom name 'C' 
                    pass
            timesteps.append(int(self.universe.trajectory.time / 1000)) # absolute timesteps are in [ps], dividing by 1000 changes to [ns]
            if bool(atomGroupLipids):
                phos0.append(self.lipid_molecules_reference['leaflet0'].positions)
                phos1.append(self.lipid_molecules_reference['leaflet1'].positions)
                leaflet0 = np.array(phos0[0])
                leaflet1 = np.array(phos1[0])
                self.lipid_molecules_positions['leaflet0'].append(leaflet0)
                self.lipid_molecules_positions['leaflet1'].append(leaflet1)
                phos0 = []
                phos1 = []
        self.max_frame = len(self.universe.trajectory)-1 # max index for frames
        self.original_max_frame = deepcopy(self.max_frame)
        self.timesteps_mainview = np.array(timesteps) 
        # store positions in numpy arrays
        if bool(atomGroupLipids):
            self.lipid_molecules_positions['leaflet0'] = np.array(self.lipid_molecules_positions['leaflet0'])
            self.lipid_molecules_positions['leaflet1'] = np.array(self.lipid_molecules_positions['leaflet1'])
        for molecule, positions in self.drug_molecules_positions.items():
            positions = np.array(self.drug_molecules_positions[molecule])
            self.drug_molecules_positions[molecule] = positions
        print('step 5: calculate lipid surface')
        # calculate polynomial represenation of membrane surface
        if bool(atomGroupLipids):
            self.box = self.calcBox()
            self.calcLipidSurfaces()
            self.lipid_molecules_surface['leaflet0'] = np.array(self.lipid_molecules_surface['leaflet0'])
            self.lipid_molecules_surface['leaflet1'] = np.array(self.lipid_molecules_surface['leaflet1'])
        print('step 6: finished calculations')

        if self.k > 1 or self.n > 0:
            # original_frame_number is NO index
            # self.max_frame: index; self.n, k: no index; + 1 convert self.max_frame to no index (= shown)
            # this calculation gives only rough results if n and k selected
            if self.k == 1 and self.n > 0:
                # only skipping n frames: we can calculate exact number of frames
                self.original_frame_number = self.max_frame + 1 + self.n
            else:
                self.original_frame_number = (self.max_frame + 1 + self.n) * self.k
        else:
            self.original_frame_number = self.max_frame + 1

    def changeSettings(self, frame_from, frame_to):
        """Copy original data for a later reverse of changes,
        calculate lipid surface abstractions and set new frame indices.
        """
        # copy original, prevent race conditions
        lock = Lock()
        with lock:
            self.original_drug_moleucles_positions = deepcopy(self.drug_molecules_positions) 
            self.original_lipid_molecules_positions = deepcopy(self.lipid_molecules_positions)
            self.original_lipid_molecules_surface = deepcopy(self.lipid_molecules_surface)
            self.original_timesteps_mainview = deepcopy(self.timesteps_mainview)
        # set current values
        for key, values in self.drug_molecules_positions.items():
            self.drug_molecules_positions[key] = values[frame_from:frame_to+1]
        self.timesteps_mainview = self.timesteps_mainview[frame_from:frame_to+1]
        for key, values in self.lipid_molecules_positions.items():
            self.lipid_molecules_positions[key] = values[frame_from:frame_to+1]
        self.lipid_molecules_surface = {}
        self.calcLipidSurfaces()
        self.min_frame = 0
        self.max_frame = frame_to - frame_from
    
    def resetSettings(self):
        """Reverse changeSettings() by copying back original data 
        and resetting the frame indices.
        """
        lock = Lock() # prevent race conditions
        with lock:
            self.drug_molecules_positions = deepcopy(self.original_drug_moleucles_positions) 
            self.lipid_molecules_positions = deepcopy(self.original_lipid_molecules_positions)
            self.lipid_molecules_surface = deepcopy(self.original_lipid_molecules_surface)
            self.timesteps_mainview = deepcopy(self.original_timesteps_mainview)
        self.min_frame = 0
        self.max_frame = deepcopy(self.original_max_frame)

    def selectFrames(self, frame_from, frame_to):
        """Read data framewise and store selected frames in a file.
        """
        self.universe.trajectory.close() # close open file
        selected_traj = []
        if frame_from == frame_to:
            selected_traj.append(mdt.load_frame(self.data, frame_from, top=self.top))
            new_traj = selected_traj[0]
        else:
            for i in range(frame_from, frame_to+1): # frame from and frame to are indices, but range() excludes end of loop
                selected_traj.append(mdt.load_frame(self.data, i, top=self.top))
            new_traj = selected_traj[0]
            i=1
            while i < len(selected_traj):
                new_traj += selected_traj[i]
                i += 1
        new_traj.save_xtc('selectedFrames.xtc') # cannot be read by MDAnalysis otherwise
        return ('selectedFrames.xtc')

    def calcBox(self):
        """Return the box coordinates from a topology file."""
        if self.top[-3:]=='gro':
            f = open(self.top, 'r')
            coordinates = f.readlines()[-1].split()
            x = float(coordinates[0])
            y = float(coordinates[1])
            z = float(coordinates[2])
            box = np.array([x,y,z])
            f.close()
        elif self.top[-3:]=='pdb':
            f = open(self.top, 'r')
            for line in f.readlines():
                if line.startswith('CRYST1'):
                    box = np.array([float(line[6:15]), float(line[15:24]), float(line[24:33])])
            f.close()
        else:
            return (np.array([0, 0, 0]))
        return box

    def calcLipidSurfaces(self):
        """Calculate polynomial representation for the membrane surface at all given timesteps."""
        self.lipid_molecules_surface['leaflet0'] = []
        self.lipid_molecules_surface['leaflet1'] = []
        for leaflet, values in self.lipid_molecules_positions.items():
            for frame in values:
                # curve_fit function changes frame data, 
                # therefore deepcopy is needed
                f = deepcopy(frame) 
                self.lipid_molecules_surface[leaflet].append(surface.curve_fit(f, self.box, self.polynomial, self.membrane_extension))


class DataSub():
    """Object holds information to be viewed in an additional view.
    
    Methods:
    -----
    setData(kwargs)
    addLine(SingleDataSet)
    addPointset(SingleDataSet)
    calculateExtremes()
    """
    def __init__(self): # , unit, unit_mpl, selection_string
        # Line plot
        self.id = ''
        self.analysis = ''
        self.unit = ''
        self.unit_matplotlib = ''
        self.unit_mpl_x = ''
        self.unit_mpl_y = ''
        self.unit_mpl_z = ''
        self.unit_mpl_t = ''
        
        # overall min, max, to determine the plotting axes limits
        self.min = {} # (x, y, z, t) # y:[x1, x2, ..]
        self.max = {} # (x, y, z, t)
        self.slider_range = (0,0)

        self.visualization = '' 
        self.lines = []
        self.pointsets_list = [] # not used currently
        self.pointsets = {}

    def setData(self, **kwargs):
        """Set the data of the DataSub object.
        
        Parameters
        --------
        labels(dict): dictionary for labels, e.g.: {pos : text},
            depends on plot to be used
        u(string): String for displaying the unit as text on a QLabel
            or similar object.
        u_mpl(string): String used to display the unit in a matplotlib
            figure object. 
        analysis(string): name of the analysis calculated
        """
        for key, value in kwargs.items():
            if key == 'labels':
                self.labels = value
            elif key == 'u':
                self.unit = value
            elif key == 'u_mpl':
                self.unit_matplotlib = value
            elif key == 'analysis':
                self.analysis = value
            elif key == 'id':
                self.id = value
            #elif key == 'vis':
            #    self.visualization = value

    def addLine(self, line):
        self.lines.append(line)

    def addPointset(self, pointset):
        self.pointsets[pointset.selection_string] = pointset

    def calculateExtremes(self):
        """min x, min y, min z over all lines / points. used for plotting axes limits"""
        x_mins = []
        x_maxs = []
        y_mins = []
        y_maxs = []
        z_mins = []
        z_maxs = []
        t_mins = []
        t_maxs = []
        
        if len(self.lines) == 0:
            reference = self.pointsets 
        else:
            reference = self.lines 

        try:
            for dp in reference:
                x_mins.append(dp.x_min)
                x_maxs.append(dp.x_max)
                y_mins.append(dp.y_min)
                y_maxs.append(dp.y_max)
                try:
                    z_mins.append(dp.z_min)
                    z_maxs.append(dp.z_max)
                except:
                    z_mins = [0]
                    z_maxs = [0]
                try: 
                    t_mins.append(dp.t_min)
                    t_maxs.append(dp.t_max)
                except: 
                    t_mins = [0]
                    t_maxs = [0]
        except:
            for name, dp in reference.items():
                x_mins.append(dp.x_min)
                x_maxs.append(dp.x_max)
                y_mins.append(dp.y_min)
                y_maxs.append(dp.y_max)
                try:
                    z_mins.append(dp.z_min)
                    z_maxs.append(dp.z_max)
                except:
                    z_mins = [0]
                    z_maxs = [0]
                try: 
                    t_mins.append(dp.t_min)
                    t_maxs.append(dp.t_max)
                except: 
                    t_mins = [0]
                    t_maxs = [0]

        x_min = np.array(x_mins).min()
        x_max = np.array(x_maxs).max()
        y_min = np.array(y_mins).min()
        y_max = np.array(y_maxs).max()
        z_min = np.array(z_mins).min()
        z_max = np.array(z_maxs).max()
        t_min = np.array(t_mins).min()
        t_max = np.array(t_maxs).max()

        self.min = (x_min, y_min, z_min, t_min)
        self.max = (x_max, y_max, z_max, t_max)


class SingleDataSet():
    """Data object, that holds information about a single set of data
    points, e.g. from a speed calculation.

    Parameter:
    ------
    x : list
        x values of the data (1st dimension)
    y : list
        y values of the data (2nd dimension)
    selection : string
        molecules/atoms identification
    color : [float] or string 
        color, specified as rgba (float) values
        or css color (string)

    optional:
    label(dict): tuple of position, text
    z(list): z values of the data (3rd dimension)
    t(list): t values of the data (4th dimension)
   
    Return:
    ----
    void

    Methods:
    ------
    addLabel(key, value): set labels(dict)with key: tuple of positions,
        value: text to show - e.g.: {(position x, position y) : text}
    """
    def __init__(self, x, y, selection, color, **kwargs):
        """Initialize SingleDataSet object with x, y, seleection and
       color attributes.
       """
        self.x = np.array(x)
        self.x_min = self.x.min()
        self.x_mini = list(np.where(self.x == self.x_min))[0]
        self.x_max = self.x.max()
        self.x_maxi = list(np.where(self.x == self.x_max))[0]
        self.x_avg = np.average(self.x) # sum(self.x)/len(self.x)
        self.y = np.array(y)
        self.y_min = self.y.min()
        self.y_mini = list(np.where(self.y == self.y_min))[0]
        self.y_max = self.y.max()
        self.y_maxi = list(np.where(self.y == self.y_max))[0]
        self.y_avg = np.average(self.y) # sum(self.y)/len(self.y)
        

        self.selection_string = selection
        self.color = color
        self.labels = {}

        for key, value in kwargs.items():
            if key == 'z':
                self.z = np.array(value)
                self.z_min = self.z.min()
                self.z_mini = list(np.where(self.z == self.z_min))[0]
                self.z_max = self.z.max()
                self.z_maxi = list(np.where(self.z == self.z_max))[0]
                self.z_avg = np.average(self.z) # sum(self.z)/len(self.z)
            elif key == 't':
                self.t = np.array(value)
                self.t_min = self.t.min()
                self.t_mini = list(np.where(self.t == self.t_min))[0]
                self.t_max = self.t.max()
                self.t_maxi = list(np.where(self.t == self.t_max))[0]
                self.t_avg = np.average(self.t) # sum(self.t)/len(self.t)

    def addLabel(self, key, value):
        """Add a label to the list of labels.
        
        Parameters:
        ------
        key: Tuple of positions
        value: String to display

        Return:
        ------
        void
        """
        self.labels[key] = value
