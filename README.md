# MeTrEx - Membrane Trajectory Exploration

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

## Installation

### Requirements
1. You need python version 3.8 or higher 
2. You need [conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-python) package manager installed on your system.

### Downloading and Installing MeTrEx
3. Clone the  [GitHub repository](https://github.com/sa-ja/MeTrEx): `git clone https://github.com/sa-ja/MeTrEx`
4. Change directory: `cd MeTrEx`
5. Build the environment: `conda env create`
6. Activate the environment: `conda activate MeTrEx`

### Run MeTrEx
7. Start MeTrEx from console: `python MeTrEx/main.py`

## HOW TO USE:

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

## TEST DATA:

The provided data files contain simulation data from a biological
membrane interacting with six molecules (polymyxins).
Interesting properties of the interaction might be: the membrane thickness
and curvature, the speed of the molecules and the intramolecular distances
between 'C'/'C1' and 'C7'.
- topology file: 'topology_6PMB.pdb'
- data file: '500_frames.xtc'
- xvg file: 'op_up_LP0_acyl_1.xvg'

## LICENCE
MeTrEx can be licenced with 'GPL 3.0 or later'.

## AUTHOR
© Christiane Rohse & Beat Ehrmann @AG Schreiber @University of Konstanz
