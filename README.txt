MeTrEx - Membrane Trajectory Exploration

MeTrEx is a program for the visual exploration of molecular simulation
data from membranes interacting with small molecules. 
It's main feature is to show an overview of the molecules' course through-
out the simulation with an abstract visualization of the membrane. This
overview of the data is shown on the 'main view', which is shown as soon
as data is loaded. 
Different analyses can be mapped onto the main view. These analyses can
also be shown in separate plots below the main view, in 'bottom views'.
Additionally, you can load other data files in 'sub windows', shown in
'sub plots'. 
Sliders and information panels give information about the currently
shown frame. Exporting data is provided for image, csv and xpdb files. 

HOW TO USE:

INSTALLATION

1. MeTrEx.exe:
Windows 10+ is required.
Open MeTrEx via double clicking the icon, no installation needed.

2. Code: 
Install the following dependencies:
moduls
python	3.7.12	conda-forge
geomdl	5.3.1	orbingol
Geomdl.cli
geomdl.shapes		
biopython	1.79	conda-forge
pyqtwebengine	5.15.7	conda-forge
networkx	2.5.1	conda-forge
pyqtchart	5.15.7	conda-forge
PyQt6	6.4.2	pypi
scipy	1.7.3	conda-forge
scikit-learn	0.22.1	conda-forge
pygemo	2.13.0	pypi
nurbspy	1.1.2	pypi
mdanalysis	2.1.0	conda-forge pypi?
mdtraj	1.9.7	conda-forge
matplotlib	3.5.3	conda-forge
cm-crameri	1.4	conda-forge

Run the main.py file: 
python main.py

USAGE

Download the folder and open the MeTrEx.exe file (maybe you need to wait
a moment for the program to start). 
In the provided data you find a topology file 'topology_6PMB.pdb' which
can be loaded with any of the trajectory files of the type '.xtc', 
e.g.: '100_frames.xtc'. 
Depending on your computer's available memory
and processor capabilities, the computation time for the first visuali-
sation differs. It might take a few seconds for small files and up to 
hours for the larger files. Therefore it might be best to start exploring
smaller files first. 
If you reduce the file size as a preprocessing step, the program will save
a reduced file as 'reducedData.xtc'.
For the exploration you can add a mapping, customize your visualizations
by choosing different colors or a different colormap. You can also change
the surface abstraction, by changing the calculation properties or don't
show an abstraction at all (= use positions of the phosphorus atoms). 

1. Open Files:
Menu > File > Open 
Crtl + O
Select structure and data file, select preprocessing. 
n = number of frames to skip at the beginning of the data
k = select every k-th frame to be shown
Molecule selection: Select the type of molecules that shall represented by its
trajectory line.
Remark: 
If a data reduction is selected, a file 'reducedData.xtc' is stored.

2. Show Analysis Mapping:
View > Map ... 

3. Show Analysis separately:
View > Show below > ...

4. Open external files:
Analysis > Show XY-XVG file

5. Save Images:
accessible via:
Above the main view: save-icon.
Right side of the BottomView: save button.
Right side of the additional windows: save button.
(When the check box in the analysis-overview-box is activatet:
also a csv file of the in the overview shown data is saved.)
The legend of the main view can be saved separately via:
Menu > File > Save legend

6. Save PDB files:
Menu > File > Save PDB
Menu > File > Save selection

7. Interactions:
Directly with the main view and the visualisations in the additional window.
Additional options are found next to the views. 

TEST DATA:

The provided data files contain simulation data from a biological
membrane interacting with six molecules (polymyxins).
Interesting properties of the interaction might be: the membrane thickness
and curvature, the speed of the molecules and the intramolecular distances
between 'C'/'C1' and 'C7'.
- topology file: 'topology_6PMB.pdb'
- data file: '500_frames.xtc'7
- xvg file: 'op_up_LP0_acyl_1.xvg'

LICENCE
MeTrEx can be licenced with 'GPL 3.0 or later'.

AUTHOR
© Christiane Rohse