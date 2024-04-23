"""This is the dialogs module.

The dialogs module contains classes for different kinds of dialogs.
Dialogs are initialized with parameters and have no return value(s).

Classes:
------
OpenDialog:
    Collect information to open files.
ChangeViewDataDialog:
    Collect a range of numbers.
ChangeColormapDialog:
    Collect a new color map (a new string from a list).
PreprocessingSelectionDialog:
    Collect values for the preprocessing of the data shown in
    the main view.
OpenXVGDialog:
    Collect file path and visualisation information for xvg files.
AboutDialog:
    Shows information about MeTrEx.
AnalysisSelectionDialog:
    Collect MDAnalysis selection strings for presented molecules and
    atoms.
Color:
    Form a square Widget showing a color.
ChangeColorDialog:
    Collects a molecule identifier and a new color.
ChangeNameDialog:
    Collects the data set to be changed and its new name.
"""

__version__ = '1.0'
__author__ = 'Christiane Rohse'

import os
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QPalette, QColor, QFont, QIntValidator, QValidator, QIcon
from PyQt6.QtWidgets import QLabel, QCheckBox, QDialog, QPushButton,\
    QDialogButtonBox, QVBoxLayout, QHBoxLayout, QLineEdit, QFileDialog, QMessageBox, QGridLayout, QWidget,\
    QGroupBox, QScrollArea, QSpinBox, QComboBox, QSizePolicy, QFrame,  QFormLayout, QListWidget
from visualisations import SubPlotCanvas

class OpenDialog(QDialog):
    """Create an open dialog.

    Methods:
    ------
    openFile():
        Check file input for validity, close dialog if valid.
    selectTopPath():
        Open file window for selection of topology file.
    selectTrajPath():
        Open file window for selection of trajectory / data file.
    rejectSelection():
        Reset variables, close dialog.
    """

    def __init__(self):
        """Show open dialog with selection possibilities for a topology
        and a data file.
        
        """
        super().__init__()
        self.setWindowTitle('Open files')

        qbtn = QDialogButtonBox.StandardButton.Open | QDialogButtonBox.StandardButton.Cancel
        self.buttonBox = QDialogButtonBox(qbtn)
        self.buttonBox.accepted.connect(self.openFile)
        self.buttonBox.rejected.connect(self.rejectSelection)

        layout = QVBoxLayout()
        label = QLabel('You need to specify a topology file and a file containing the simulation data:')

        layoutTop = QHBoxLayout()
        self.textTop = QLineEdit()
        self.textTop.setEnabled(False)
        self.textTop.setText('select a topology file')
        self.topPath = None
        selectButton1 = QPushButton('select')
        selectButton1.clicked.connect(self.selectTopPath)
        layoutTop.addWidget(self.textTop)
        layoutTop.addWidget(selectButton1)

        layoutTraj = QHBoxLayout()
        self.textTraj = QLineEdit()
        self.textTraj.setEnabled(False)
        self.textTraj.setText('select a data file')
        self.trajPath = None
        selectButton2 = QPushButton('select')
        selectButton2.clicked.connect(self.selectTrajPath)
        layoutTraj.addWidget(self.textTraj)
        layoutTraj.addWidget(selectButton2)

        layout.addWidget(label)
        layout.addLayout(layoutTop)
        layout.addLayout(layoutTraj)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

    def openFile(self):
        """Test if all needed variables are set and accept."""
        if self.topPath == None or self.trajPath == None:
            messageBox = QMessageBox(self)
            messageBox.setWindowTitle('Caution!')
            messageBox.setText('You need to specify a topology file and a data file.')
            messageBox.setStandardButtons(QMessageBox.StandardButton.Ok)
            self.buttonBox = messageBox.exec()   
        # check compatibility of files specified
        elif self.topPath[-3:]=='pdb' and self.trajPath[-3:]=='trr':  #self.topPath[-3:]=='gro' and self.trajPath[-3:]=='xtc' or self.topPath[-3:]=='pdb' and self.trajPath[-3:]=='trr'
            messageBox = QMessageBox(self)
            messageBox.setWindowTitle('Caution!')
            messageBox.setText('Your selected file types do not match.')
            messageBox.setStandardButtons(QMessageBox.StandardButton.Ok)
            self.buttonBox = messageBox.exec()  
        else:
            self.accept()

    def selectTopPath(self):
        """Open file window for selection of topology file."""
        file_filter = 'Structure File (*.pdb *.gro)'
        response = QFileDialog.getOpenFileName(
                    parent=self,
                    caption='Select a file',
                    directory=os.getcwd(),
                    filter=file_filter,
                    initialFilter='Structure File (*.pdb *.gro)'
                    )
        self.textTop.setText(str(response[0]))
        self.topPath = response[0]

    def selectTrajPath(self):
        """Open file window for selection of data file."""
        file_filter = 'Simulation Data File (*.xtc *.trr)'
        response = QFileDialog.getOpenFileName(
                    parent=self,
                    caption='Select a file',
                    directory=os.getcwd(),
                    filter=file_filter,
                    initialFilter='Simulation Data File (*.xtc *.trr)'
                    )
        self.textTraj.setText(str(response[0]))
        self.trajPath = response[0]

    def rejectSelection(self):
        """Set file paths to None and close the dialog."""
        self.trajPath = None
        self.topPath = None
        self.reject()


class ChangeViewDataDialog(QDialog):
    """Dialog to select needed parameters for updating the main view.
    
    Parameter:
    ------
    parent : QObject
        Parent of this dialog.
    from_value : int
        Selection range start index (will not be adjusted, use index
        to be presented to the user).
    to_value : int
        Selection range stop index (will not be adjusted, use index
        to be presented to the user).
    optional:
    save : bool
        Collect file name (e.g. to save the selection).
    current : int
        Currently presented frame in the main view. 

    Methods:
    ------
    checkTypeRange():
        Checks the validity of the input.
    acceptChanges():
        Accepts the input, closes dialog.
    """

    def __init__(self, parent, from_value, to_value, **kwargs):
        """Show selection dialog."""
        super().__init__(parent)
        self.parent = parent
        self.setWindowTitle('Choose settings')
        self.selection_accepted = False
        self.save = False
        self.from_value = None
        self.to_value = None
        self.f = from_value
        self.t = to_value
        self.current = '1'

        for key, value in kwargs.items():
            if key == 'save':
                self.save = value
            if key == 'current':
                self.current = value

        select_frames = QGroupBox('Select frames indices')
        select_frames.setAccessibleName('select_frames')
        combined_layout = QVBoxLayout()
        description_layout = QHBoxLayout()
        select_frames_layout = QHBoxLayout()

        label_from = QLabel('from')
        label_to = QLabel('to')
        self.frame_from = QLineEdit()
        self.frame_from.setEnabled(True)
        self.frame_from.editingFinished.connect(self.checkTypeRange)
        self.frame_to = QLineEdit()
        self.frame_to.setEnabled(True)
        self.frame_to.editingFinished.connect(self.checkTypeRange)
        self.frame_from.setText(self.current)
        self.frame_to.setText(self.current) 

        description_layout.addWidget(label_from)
        description_layout.addWidget(label_to)
        select_frames_layout.addWidget(self.frame_from)
        select_frames_layout.addWidget(self.frame_to)
        combined_layout.addLayout(description_layout)
        combined_layout.addLayout(select_frames_layout)
        select_frames.setLayout(combined_layout)

        if self.save: 
            qbtn = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        else:
            qbtn = QDialogButtonBox.StandardButton.Open | QDialogButtonBox.StandardButton.Cancel
        self.buttonBox = QDialogButtonBox(qbtn)
        self.buttonBox.accepted.connect(self.acceptChanges)
        self.buttonBox.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addWidget(select_frames)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

    def checkTypeRange(self):
        """Checks the type and range of an input."""
        validation_rule = QIntValidator(self.f, self.t)
        # 42: no idea what that is for
        if validation_rule.validate(self.frame_from.text(), 42)[0] == QValidator.State.Acceptable and validation_rule.validate(self.frame_to.text(), 42)[0] == QValidator.State.Acceptable:
            self.selection_accepted = True
        elif not validation_rule.validate(self.frame_from.text(), 42)[0] == QValidator.State.Acceptable:
            self.frame_from.setText('1')
            self.selection_accepted = False
        elif not validation_rule.validate(self.frame_to.text(), 42)[0] == QValidator.State.Acceptable:
            self.frame_to.setText('1')
            self.selection_accepted = False
        else: # probably not reachable
            self.frame_from.setText('1')
            self.frame_to.setText('1')
            self.selection_accepted = False
 
    def acceptChanges(self):
        """Check input values for validity and set return values."""
        if (int(self.frame_to.text()) - int(self.frame_from.text())) < 0:
            self.selection_accepted = False
            messageBox = QMessageBox(self)
            messageBox.setWindowTitle('Caution!')
            messageBox.setText('\'To\' index is too small.')
            messageBox.setStandardButtons(QMessageBox.StandardButton.Ok)
            self.buttonBox = messageBox.exec()
        elif not self.selection_accepted:
            self.reject()
        else:
            self.from_value = int(self.frame_from.text())
            self.to_value = int(self.frame_to.text())
            self.accept()    


class PreprocessingSelectionDialog(QDialog):
    """Dialog to select k, n, line representation molecules.
    
    Parameter:
    ------
    molecules : set
        Set, containing molecule names (str) to be chosen from.
    nkonly : bool
        True, if only n and k shall be chosen, else also molecules.

    Methods:
    ------
    setSelection():
        Accept the selection, close the dialog.    

    Notes:
    ------
    k : 
        Every k-th frame will be selected, default: 1. Max. 1000000.
    n : 
        n frames will be skipped, default: 0. Max. 1000000.
    """

    def __init__(self, molecules=set(), nkonly=False):
        """Show selection dialog for frame reduction and visualization settings."""
        super().__init__()
        self.setWindowTitle('Choose settings')
        self.n = 0
        self.k = 1
        self.selected_molecules = None
        self.nkonly = nkonly
        self.anker_selection = False
        self.anker_selection_text = 'select manually'


        preprocess_groupbox = QGroupBox('Data Reduction')
        wrap_layout = QVBoxLayout()
        preprocess_layout = QHBoxLayout()
        bold = QFont()
        bold.setBold(True)
        preprocessing_caution = QLabel('Be aware:')
        preprocessing_caution.setFont(bold)
        preprocessing_explanation = QLabel(
            'This following setting FIRST takes every k-th frame\n'\
            +'of your data file and AFTERWARDS skips n frames!'
            )
        italic = QFont()
        italic.setItalic(True)
        spin_k_layout = QVBoxLayout()
        spin_k_label = QLabel('Specify every k-th frame to select')
        spin_k_label.setFont(italic)
        spin_k = QSpinBox()
        spin_k.setAccessibleName('spin_k')
        spin_k.setRange(1, 1000000)
        spin_k.setValue(1)
        spin_k.setSingleStep(100)
        spin_k_layout.addWidget(spin_k_label)
        spin_k_layout.addWidget(spin_k)
        spin_n_layout = QVBoxLayout()
        spin_n_label = QLabel('Select the number of frames\n to skip at the beginning')
        spin_n_label.setFont(italic)
        spin_n = QSpinBox()
        spin_n.setAccessibleName('spin_n')
        spin_n.setRange(0, 1000000)
        spin_n.setValue(0)
        spin_n.setSingleStep(10)
        spin_n_layout.addWidget(spin_n_label)
        spin_n_layout.addWidget(spin_n)
        preprocess_layout.addLayout(spin_k_layout)
        preprocess_layout.addLayout(spin_n_layout)
        wrap_layout.addWidget(preprocessing_caution)
        wrap_layout.addWidget(preprocessing_explanation)
        wrap_layout.addLayout(preprocess_layout)
        preprocess_groupbox.setLayout(wrap_layout)
        
        if not self.nkonly:
            vis_selection_groupbox = QGroupBox('Visualization')
            vis_selection_layout = QVBoxLayout()
            label = QLabel('You need to specify the name of the molecules for the line representation:')
            vis_selection_layout.addWidget(label)
            for mol in sorted(molecules): #FIX by Beat, added sorting
                self.checkBtn = QCheckBox(mol)
                vis_selection_layout.addWidget(self.checkBtn)
            vis_selection_groupbox.setLayout(vis_selection_layout)
            
            anker_selection_groupbox = QGroupBox('Anker Atom')
            anker_selection_layout = QVBoxLayout()
            anker_label = QLabel('Do you want to manually select an anker atom for the line representation?')
            anker_selection_layout.addWidget(anker_label)
            self.checkAnker = QCheckBox(self.anker_selection_text)
            anker_selection_layout.addWidget(self.checkAnker)
            anker_selection_groupbox.setLayout(anker_selection_layout)
            
        qbtn = QDialogButtonBox.StandardButton.Ok
        self.button_box = QDialogButtonBox(qbtn)
        self.button_box.accepted.connect(self.setSelection)

        layout = QVBoxLayout()
        layout.addWidget(preprocess_groupbox)
        if not self.nkonly:
            layout.addWidget(vis_selection_groupbox)
            layout.addWidget(anker_selection_groupbox)
        layout.addWidget(self.button_box)
        self.setLayout(layout)

    def setSelection(self):
        """Process selections at confirmation."""
        self.selected_molecules = set()
        if not self.nkonly:
            for child in self.findChildren(QCheckBox):
                if child.isChecked() == True: #Check for selection
                    if child.text() == self.anker_selection_text: #Check if selection is for manual anker selection
                        self.anker_selection = True
                    else:
                        self.selected_molecules.add(child.text())
        for child in self.findChildren(QSpinBox):
            if child.accessibleName() == 'spin_n':
                self.n = child.value()
            elif child.accessibleName() == 'spin_k':
                self.k = child.value()
        self.accept()


class ChangeColormapDialog(QDialog):
    """Show a dialog for the user to select a colormap for all colors.
    
    Parameter:
    ------
    ccmap: string
        Currently used colour map.
    colormaps[]: list 
        Colour maps to choose from.    

    Methods:
    ------
    acceptChanges():
        Sets internal variable, closes the dialog.
    """
    def __init__(self, ccmap,colormaps=[]):
        super().__init__()
        title = 'Colormap selection'
        self.setWindowTitle(title)
        self.colormap = ''
        layout = QVBoxLayout()

        selection_groupbox = QGroupBox('Select a colormap')
        selection_layout = QVBoxLayout()
        current_colormap = QLabel('Current colormap:')
        current_name = QLabel(ccmap)
        current_name_font = QFont()
        current_name_font.setBold(True)
        current_name.setFont(current_name_font)
        new_colormap = QLabel('New colormap:')
        self.colormap_selection = QComboBox()
        self.colormap_selection.addItems(colormaps)

        selection_layout.addWidget(current_colormap)
        selection_layout.addWidget(current_name)
        selection_layout.addWidget(new_colormap)
        selection_layout.addWidget(self.colormap_selection)
        selection_groupbox.setLayout(selection_layout)

        qbtn = QDialogButtonBox.StandardButton.Open | QDialogButtonBox.StandardButton.Cancel
        self.buttonBox = QDialogButtonBox(qbtn)
        self.buttonBox.accepted.connect(self.acceptChanges)
        self.buttonBox.rejected.connect(self.reject)

        layout.addWidget(selection_groupbox)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)
        
    def acceptChanges(self):
        """Collect new color map."""
        self.colormap = self.colormap_selection.currentText()
        self.accept()
        

class OpenXVGDialog(QDialog):
    """Show a dialog to select a data file, set a string for the title 
    and select a referred to molecule.

    Parameter:
    ------
    molecules : dict
        Molecule identifiers (dict keys) will be shown to the user to 
        select from. 

    Methods:
    ------
    selectPath():
        Open a window to select a file path. File filter is used to
        show xvg files only.
    acceptChanges():
        Test validity of user input, save variables and close window
        if input correct, otherwise show a message and return to the
        dialog.
    rejectChanges():
        Reset variables, close dialog.
    """

    def __init__(self, molecules):
        """Show open dialog with selection possibilities for an 
        xvg file and optional reference to drug molecules.
        """
        super().__init__()
        title = 'XVG selection'
        self.setWindowTitle(title)
        self.file_path = None
        self.drug_molecules = []
        self.dim3 = False
        layout = QVBoxLayout()

        qbtn = QDialogButtonBox.StandardButton.Open | QDialogButtonBox.StandardButton.Cancel
        self.buttonBox = QDialogButtonBox(qbtn)
        self.buttonBox.accepted.connect(self.acceptChanges)
        self.buttonBox.rejected.connect(self.rejectChanges)
        # File path
        path_layout = QHBoxLayout()
        self.path_string = QLineEdit()
        self.path_string.setEnabled(False)
        self.path_string.setText('select an xvg file')
        selection_button = QPushButton('select')
        selection_button.clicked.connect(self.selectPath)
        path_layout.addWidget(self.path_string)
        path_layout.addWidget(selection_button)
        # File Title
        self.file_title = QLineEdit()
        self.file_title.setText('Fill in the title for the new window')
        # Use index = frame no as 3rd dimension
        dim = QCheckBox('use frame number as 3rd dimension')
        dim.setToolTip('Useful for Eigenvektor projections.\nAplicable if 2 or more colums in xvg file.')
        # Reference drug molecules for later choosing color
        # in two columns.
        drug_molecules_layout = QGridLayout()
        i=0
        j=0
        for key, value in molecules.items():
            molecule = QCheckBox(key)
            drug_molecules_layout.addWidget(molecule, i, j, Qt.AlignmentFlag.AlignLeft)
            j += 1
            if j == 2:
                i += 1
                j = 0
        layout.addLayout(path_layout)
        layout.addWidget(self.file_title)
        layout.addWidget(dim)
        layout.addLayout(drug_molecules_layout)
        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

    def selectPath(self):
        """Open file window for selection of xvg file."""
        file_filter = 'XVG File (*.xvg)'
        response = QFileDialog.getOpenFileName(
                    parent=self,
                    caption='Select a file',
                    directory=os.getcwd(),
                    filter=file_filter,
                    initialFilter='XVG File (*.xvg)'
                    )
        self.path_string.setText(str(response[0]))
        self.file_path = response[0]

    def acceptChanges(self):
        """Check if all needed information has been set (file path), if no title set use file name as title, close dialog."""
        if self.file_path == None or self.file_path == '':
            messageBox = QMessageBox(self)
            messageBox.setWindowTitle('Caution!')
            messageBox.setText('You need to specify an xvg file to show its data.')
            messageBox.setStandardButtons(QMessageBox.StandardButton.Ok)
            self.buttonBox = messageBox.exec()  
        else:
            for child in self.findChildren(QCheckBox):
                if child.isChecked() and child.text() != 'use frame number as 3rd dimension':
                    self.drug_molecules.append(child.text())
                elif child.isChecked() and child.text() == 'use frame number as 3rd dimension':
                    self.dim3 = True
            if self.file_title.text() == 'Fill in the title for the new window':
                string = self.file_path.split('/')
                self.file_title.setText(string[len(string)-1])
            self.accept()

    def rejectChanges(self):
        """Set self.file_path to None to indicate that no file shall be loaded and shown, close dialog."""
        self.file_path = None
        self.reject()

class AboutDialog(QDialog):
    """Dialog to display information about the program."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle('About MeTrEx')
        self.setFixedSize(250, 220)

        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout(self)
        self.setLayout(layout)

        # Name, version, and license labels
        name_label = QLabel('MeTrEx')
        version_label = QLabel('Version 1.0')
#        license_label = QLabel('GPL-3')
        license_label = QLabel()
        license_label.setPixmap(QIcon.fromTheme("text-x-generic").pixmap(16, 16))
        license_label.setText('<font color="#0074D9"><b>GPL-3</b></font>')
#        author_label = QLabel('\u00A9 C. Rohse & B. Ehrmann')

        # Set label alignments
        for label in (name_label, version_label, license_label):
            label.setAlignment(Qt.AlignmentFlag.AlignHCenter)

        # Add labels to a wrapper widget
        wrapper_widget = QWidget()
        wrapper_layout = QVBoxLayout(wrapper_widget)
        wrapper_layout.addSpacing(50)
        wrapper_layout.addWidget(name_label)
        wrapper_layout.addWidget(version_label)
        wrapper_layout.addWidget(license_label)
#        wrapper_layout.addWidget(author_label)
        wrapper_layout.addSpacing(40)

        # Ok button
        ok_button = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok)
        ok_button.accepted.connect(self.accept)

        # Add widgets to the main layout
        layout.addWidget(wrapper_widget)
        layout.addWidget(ok_button)


class DocumentationDialog(QDialog):
    """Dialog to display program documentation."""
    
    # Class variable for text content
    TEXT_CONTENT = (
            """
MeTrEx - Membrane Trajectory Exploration

MeTrEx is a program for the visual exploration of molecular simulation data from membranes interacting with small molecules.
Its main feature is to show an overview of the molecules' course throughout the simulation with an abstract visualization of the membrane.
This overview of the data is shown on the 'main view', which is shown as soon as data is loaded.
Different analyses can be mapped onto the main view. These analyses can also be shown in separate plots below the main view, in 'bottom views'.
Additionally, you can load other data files in 'sub windows', shown in 'sub plots'.
Sliders and information panels give information about the currently shown frame. Exporting data is provided for image, csv and xpdb files.

USAGE
Depending on your computer's available memory and processor capabilities, the computation time for the first visualisation differs. It might take a few seconds for small files and up to hours for the larger files. Therefore it might be best to start exploring smaller files first.
If you reduce the file size as a preprocessing step, the program will save a reduced file as 'reducedData.xtc'.
For the exploration you can add a mapping, customize your visualizations by choosing different colors or a different colormap. You can also changethe surface abstraction, by changing the calculation properties or don't show an abstraction at all (= use positions of the phosphorus atoms).

1. Open Files:
Menu > File > Open
Crtl + O
Select structure and data file, select preprocessing.
n = number of frames to skip at the beginning of the data
k = select every k-th frame to be shown
Molecule selection: Select the type of molecules that shall represented by its trajectory line.
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
(When the check box in the analysis-overview-box is activatet also a csv file of the in the overview shown data is saved.)
The legend of the main view can be saved separately via:
Menu > File > Save legend

6. Save PDB files:
Menu > File > Save PDB
Menu > File > Save selection

7. Interactions:
Directly with the main view and the visualisations in the additional window.\n
Additional options are found next to the views.

LICENCE
MeTrEx can be licensed with 'GPL 3.0 or later'.
AUTHOR
© Christiane Rohse & Beat Ehrmann
            """
        )

    def __init__(self):
        super().__init__()
        self.setWindowTitle('Documentation')
        self.resize(500, 400)  # Adjusted size for better fit

        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout(self)
        self.setLayout(layout)

        # Scroll area for text content
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)  # Allow the content to resize with the scroll area
        layout.addWidget(scroll_area)

        # Widget for the text content
        content_widget = QWidget()
        scroll_area.setWidget(content_widget)
        content_layout = QVBoxLayout(content_widget)

        # Info label
        info_label = QLabel(self.TEXT_CONTENT)
        info_label.setWordWrap(True)  # Allow text wrapping
        content_layout.addWidget(info_label)

        # Ok button
        button_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok)
        button_box.accepted.connect(self.accept)
        layout.addWidget(button_box, alignment=Qt.AlignmentFlag.AlignRight)

        
class AnalysisSelectionDialog(QDialog):
    """Selection Dialog for selecting molecule types, single molecules,
    atoms needed for analysis.

    Parameter:
    --------
    analysis: 'drug_molecules', any str. If 'drug_molecules' is 
        used, each molecule shown in the line representation is
        listed and no other molecules can be selected.
    **kwargs: 
        intra (bool): True: show additional field for second atom 
            of the same molecule.
        !atom (bool): True: set atom selection entries in the
            drop down menu in the __init__ function (static), 
            False: calculete selection entries dynamically with the
            selected residue numer.
        atoms (list): List of all atoms of the universe, given as 
            string: 'moleculetype'+'residuenumber'+'atom'
        atoms_dict (dict): Dictionary containing a list of atoms
            for each distinct Molecule. The list starts with 'all'.
        distinct_mol (set): Set containing the distinct molecules
            identified by resid name and resid number.
        distinct_mol_t (set):Set containing the distinct molecule
            types identified by resid name.
        drug_mol (dict): Dictionary containing the drug molecules
            (line represenation molecules).
        universe (Universe): MDAnalysis Universe object of the data
            from which to select.
        visual (bool): True: show the visual selection panel.
        rad (bool): True: show radius selection.
    
    Methods:
    ------
    acceptSelection():
        Check input vor validity, if correct close dialog, else inform
        the user and return to dialog.
    adjustAtomSelection():
        Adjusts the atom selection drop down menu on-the-fly depending
        on the residue number selected.
    addToSelection():
        Adds the selection string of a molecule to the summary.
    updateSummary():
        Updates the summary whenever a string is selected or deselected.
    removeFromSelection():
        Remove a string from the selection set.
    showMessage():
        Show a message to the user.
    rejectSelection():
        Reset variables, close the dialog.
    getCurrentSelectionString():
        Returns the selection string of a selection line, after
        selection via the '+' and '-' buttons.

    """
    def __init__(self, analysis, **kwargs):
        """Show selection dialog for molecules and atoms."""
        super().__init__()
        self.setWindowTitle('Choose analysis settings')
        self.setMinimumWidth(400)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding)
        self.intra = False # 2 atoms per molecule
        self.selected_molecules = set()
        self.resid = True # show resid numbers to select
        self.atom = False # 
        self.atoms = []
        self.atoms_dict = {}
        self.visual_selection = False # shall the visual selection be shown?
        self.radius = False # radius selection
        self.intra_all = False

        for arg, value in kwargs.items():
            if arg == 'intra':
                self.intra = value
            elif arg == 'atom': # ?
                self.atom = value
            elif arg == 'atoms':
                self.atoms = value
            elif arg == 'atoms_dict':
                self.atoms_dict = value
            elif arg == 'distinct_mol':
                self.distinct_molecules = value
            elif arg == 'distinct_mol_t':
                self.distinct_molecule_types = value
            elif arg == 'drug_mol':
                self.drug_molecules = value
            elif arg == 'universe':
                self.universe = value
            elif arg == 'visual':
                self.visual_selection = value
            elif arg == 'windowtitle':
                self.setWindowTitle(value)
            elif arg == 'rad':
                self.radius = value
            elif arg == 'intra_all':
                self.intra_all = value

        # MANUAL SELECTION
        manual_selection = QGroupBox('Manual selection')
        manual_selection.setAccessibleName('manual_selection')
        manual_selection_intermediate_layout = QVBoxLayout()
        manual_selection_scrollarea = QScrollArea()
        manual_selection_scrollarea.setFrameShape(QFrame.Shape.NoFrame)
        manual_selection_scrollarea_widget = QWidget()
        manual_selection_layout = QVBoxLayout(manual_selection_scrollarea_widget)
        # populate selection for each distinct type of molecule
        if analysis == 'drug_molecules':
            reference = self.drug_molecules
            self.resid = False
            self.atom = True
        else:
            reference = self.distinct_molecule_types
            self.resid = True
            self.atom = False
        if self.radius: # make sure you cannot select senseless stuff
            self.intra = False
        # populate selection combos an buttons
        for molecule in sorted(reference):
            molecules_layout = QHBoxLayout()
            molecule_type_name = QLabel(molecule)
            if self.resid:
                molecule_resid = QComboBox()
                molecule_resid_numbers = []
                for mol in self.distinct_molecules:
                    if mol.startswith(molecule):
                        molecule_resid_numbers.append(int(mol[len(molecule):]))
                molecule_resid_numbers.sort()
                molecule_resid_numbers = list(map(str, molecule_resid_numbers))
                molecule_resid_numbers.insert(0, 'all')
                molecule_resid.addItems(molecule_resid_numbers)
                molecule_resid.currentIndexChanged.connect(self.adjustAtomsSelection)
                molecule_resid.setAccessibleName(molecule+'resid')
                molecule_resid.setMinimumWidth(60)
                molecule_resid.adjustSize()
            # select atom
            molecule_atoms = QComboBox() # add atom names to select by selection of residue number
            molecule_atoms.addItems({'all'})
            molecule_atoms.setAccessibleName(molecule+'atom')
            molecule_atoms.setMinimumWidth(60)
            molecule_atoms.adjustSize()
            if self.atom:
                molecule_atoms.clear()
                molecule_atoms.addItems(self.atoms_dict[molecule])
            if self.intra:
                molecule_atoms2 = QComboBox() # add atom names to select by selection of residue number
                if analysis == 'drug_molecules':
                    molecule_atoms2.addItems(self.atoms_dict[molecule])
                else:
                    molecule_atoms2.addItems({'all'})
                molecule_atoms2.setAccessibleName(molecule+'atom2')
                molecule_atoms2.setMinimumWidth(60)
                molecule_atoms2.adjustSize() 
            # add and remove selection
            if self.radius:
                unit_text = '['.encode('utf-8') + u'\u00C5'.encode('utf-8') + ']'.encode('utf-8')
                text = 'r ' + unit_text.decode('utf-8') + ':'
                radius_label = QLabel(text)
                radius_label.setToolTip('Molecules within given \n radius will be selected.')
                radius_label.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
                radius = QSpinBox()
                radius.setAccessibleName(molecule + 'radius')
                radius.setMinimumWidth(60)
            button_add = QPushButton('+')
            button_add.setFixedWidth(25)
            button_add.setAccessibleName(molecule)
            button_add.setCheckable(True)
            button_add.clicked.connect(self.addToSelection)
            button_remove = QPushButton('-')
            button_remove.setFixedWidth(25)
            button_remove.setAccessibleName(molecule)
            button_remove.clicked.connect(self.removeFromSelection)
            # add all boxes and buttons to the layout
            molecules_layout.addWidget(molecule_type_name)
            if self.resid:
                molecules_layout.addWidget(molecule_resid)
            molecules_layout.addWidget(molecule_atoms)
            if self.intra:
                molecules_layout.addWidget(molecule_atoms2)
            if self.radius:
                molecules_layout.addWidget(radius_label)
                molecules_layout.addWidget(radius)
            molecules_layout.addWidget(button_add)
            molecules_layout.addWidget(button_remove)
            manual_selection_layout.addLayout(molecules_layout)
        manual_selection_scrollarea.setWidgetResizable(True)
        manual_selection_scrollarea.setWidget(manual_selection_scrollarea_widget)
        manual_selection_intermediate_layout.addWidget(manual_selection_scrollarea)
        manual_selection.setLayout(manual_selection_intermediate_layout)
        # manual_selection.setLayout(manual_selection_layout)
        j = len(reference)*40
        if j < 400:
            manual_selection.setMinimumSize(360, j)
        else: 
            manual_selection.setMinimumSize(360, 400)
        # atom list dictionary, containing molecule + list of corresponding atoms

        # VISUAL SELECTION
        vis_selection = QGroupBox('Visual Selection')
        if not self.visual_selection:
            vis_selection.setHidden(True)
        vis_selection.setAccessibleName('vis_selection')
        vis_selection_layout = QVBoxLayout()
        vis_check_button_layout = QHBoxLayout()
        check_btn_data = QCheckBox('Data')
        check_btn_data.setAccessibleName('vis_data')
        check_btn_upper = QCheckBox('Upper Leaflet')
        check_btn_upper.setAccessibleName('vis_upper')
        check_btn_lower = QCheckBox('Lower Leaflet')
        check_btn_lower.setAccessibleName('vis_lower')
        vis_check_button_layout.addWidget(check_btn_data)
        vis_check_button_layout.addWidget(check_btn_upper)
        vis_check_button_layout.addWidget(check_btn_lower)
        visc = SubPlotCanvas(self, width=2, height=2, dpi=200) # add needed vis here
        #visc.scatterPlotSelection(self.parent.data) 
        vis_selection_layout.addLayout(vis_check_button_layout)
        vis_selection_layout.addWidget(visc)
        vis_selection.setLayout(vis_selection_layout)

        # SELECTION SUMMARY
        selection_summary = QGroupBox('Summary')
        selection_summary_layout = QVBoxLayout()
        summary_label = QLabel('You selected the following molecules/atoms:')
        self.summary_list = QLabel() # change text according to selection
        selection_summary_layout.addWidget(summary_label, Qt.AlignmentFlag.AlignTop)
        selection_summary_layout.addWidget(self.summary_list)
        selection_summary_layout.addStretch(3)
        selection_summary.setLayout(selection_summary_layout)

        self.selection_layout = QGridLayout()
        self.selection_layout.addWidget(manual_selection, 0, 0)
        self.selection_layout.addWidget(selection_summary, 1, 0)
        self.selection_layout.addWidget(vis_selection, 0, 1, 2, 1)

        qbtn = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        self.button_box = QDialogButtonBox(qbtn)
        self.button_box.accepted.connect(self.acceptSelection)
        self.button_box.rejected.connect(self.rejectSelection)

        layout = QVBoxLayout()
        layout.addLayout(self.selection_layout)
        layout.addWidget(self.button_box)
        self.setLayout(layout)

    def acceptSelection(self):
        """Check input for validity and close the dialog if correct."""
        n = 0
        for child in self.findChildren(QPushButton):
            if child.text() == '+':
                n += 1
        if self.intra_all:
            if len(self.selected_molecules) < n:
                self.showMessageBox('You need to specify atoms for all molecules.')
                return
        self.accept()

    def adjustAtomsSelection(self):
        """Adjust atom selection list with respect of the selected
        residue number.

        """
        sender = self.sender()
        name = sender.accessibleName() # molecule + 'resid'
        molecule_no = sender.currentText() # resid number
        for box in self.findChildren(QGroupBox):
            if box.accessibleName() == 'manual_selection':
                kids = box.findChildren(QComboBox)
                molecule_name = name[:-5]
                for kid in kids:
                    if kid.accessibleName() == molecule_name + 'atom' or kid.accessibleName() == molecule_name + 'atom2':
                        kid.clear()
                        if molecule_no == 'all':
                            kid.addItem('all')
                        else:
                            kid.addItems(self.atoms_dict[molecule_name + molecule_no])

    def addToSelection(self):
        """Add moleule/atom to the selection."""
        selection_string = self.getCurrentSelectionString()
        if self.intra:
            if len(selection_string.split()) < 17 or selection_string.split()[16] == 'all' or selection_string.split()[16] == selection_string.split()[7]:
                self.sender().setChecked(False)
                self.showMessageBox('You need to specify a distinct molecule \n and two distinct atoms.')
                return
            if self.intra_all:
                self.sender().setEnabled(False)
        self.selected_molecules.add(selection_string)
        self.updateSummary()
        self.sender().setChecked(False)
    
    def updateSummary(self):
        """Update the summary after adding / removing molecules or 
        atoms to / from the selection. and show the new in the 
        corresponding label.

        """
        summary = ''
        for item in self.selected_molecules:
            string_list = item.split()
            #print('sel_list: ', string_list, '\n', 'length: ', len(string_list))
            length = len(string_list)
            if length == 2: # only resname
                summary += item[8:] + '\n'
            elif length == 5: # resname and resid
                summary += string_list[1] + ', ' + string_list[4] + '\n'
            elif length == 4: # resname + radius
                summary += string_list[3][:-1] + ', r: ' + string_list[1] + '\n'
            elif length == 7 : # resname and resid + radius
                summary += string_list[3] + ', ' + string_list[6][:-1] + ', r: ' + string_list[1] + '\n'
            elif length == 10: # resname, resid and atom + radius
                summary += string_list[3] + ', ' + string_list[6] + ', ' + string_list[9][:-1] + ', r: ' + string_list[1] + '\n'
            elif length == 8: # resname, resid and atom
                summary += string_list[1] + ', ' + string_list[4] + ', ' + string_list[7] + '\n'
            else: # length == 17: resname, resid and atom ; resname, resid and atom
                summary += string_list[1] + ', ' + string_list[4] + ', ' + string_list[7] + ' + ' + string_list[16] + '\n'
        self.summary_list.setText(summary)

    def removeFromSelection(self):
        """Remove molecule/atom form the selection."""
        selection_string = self.getCurrentSelectionString()
        update_add_button = False
        if selection_string in self.selected_molecules:
            update_add_button = True
        try:
            self.selected_molecules.remove(selection_string)
            self.updateSummary()
        except:
            pass
        self.sender().setChecked(False)
        if self.intra_all and update_add_button:
            name = self.sender().accessibleName()
            for child in self.findChildren(QPushButton):
                if child.text() == '+' and child.accessibleName() == name:
                    child.setEnabled(True)
    
    def showMessageBox(self, text):
        """Show a 'Caution' message box to the user.
        
        Parameter:
        ------
        text : string
            Instruction shown to the user.
        """
        message_box = QMessageBox.warning(
            self, 
            'Caution!',
            text,
            QMessageBox.StandardButton.Ok,
            QMessageBox.StandardButton.Ok
            )

    def rejectSelection(self):
        """Reset variables, close the dialog."""
        self.selected_molecules = set()
        self.reject()

    def getCurrentSelectionString(self):
        """Returns the currently selected string of a line after the
        '+' or '-' button was pressed.

        """
        sender = self.sender()
        name = sender.accessibleName()
        if not self.resid:
            mol = self.drug_molecules[name][0]
            resid_text = self.drug_molecules[name][1]
        else:
            mol = name
        atom_text = ''
        atom2_text = ''
        for box in self.findChildren(QGroupBox):
            if box.accessibleName() == 'manual_selection':
                kids = box.findChildren(QComboBox)
                resid = name+'resid'
                atom = name+'atom'
                atom2 = name+'atom2'
                for kid in kids:
                    if kid.accessibleName() == resid:
                        resid_text = kid.currentText()
                    if kid.accessibleName() == atom:
                        atom_text = kid.currentText()
                    if kid.accessibleName() == atom2:
                        atom2_text = kid.currentText()
                if self.radius:
                    radii = box.findChildren(QSpinBox)
                    for rad in radii:
                        if rad.accessibleName() == name + 'radius':
                            radius = rad.value()
                            #print('raidus is:', radius)
        if atom_text == 'all' and resid_text == 'all': 
            selection_string = 'resname ' + mol
        elif resid_text == 'all': # atom text != all
            selection_string = 'resname ' + mol + ' and name ' + atom_text
        elif atom_text == 'all':
            selection_string = 'resname ' + mol + ' and resid ' + resid_text
        elif atom2_text != '':
            selection_string = 'resname ' + mol + ' and resid ' + resid_text + ' and name ' + atom_text \
                + ' ; '+ 'resname ' + mol + ' and resid ' + resid_text + ' and name ' + atom2_text
        else:
            selection_string = 'resname ' + mol + ' and resid ' + resid_text + ' and name ' + atom_text
        if self.radius and radius > 0:
            selection_string = '(around ' + str(radius) + ' ' + selection_string + ')'
        return selection_string


class Color(QWidget):
    """Create a plain widget that displays a color specified as list 
    of rgb values.

    Parameter:
    ------
    color (list):
        Color of the Widget, given as list of float values.
    size (int):
        Width of the widget

    Methods:
    ------
    setAccessibleName():
        Set the name of the Color object.

    """

    def __init__(self, color, size):
        """Initialize QWidget with color and fixed width.
        
        """
        super(Color, self).__init__()
        self.setAutoFillBackground(True)
        self.color = color
        palette = self.palette()
        palette.setColor(QPalette.ColorRole.Window, QColor(color[0], color[1], color[2]))
        self.setPalette(palette)
        self.setFixedWidth(size)

    def setAccessibleName(self, name):
        """Set the name of the Color object."""
        self.accessibleName = name


class ChangeColorDialog(QDialog):
    """Open a dialog to ask the user for the molecule and the color to
    change. Returns molecule name (identifier) and color.

    Parameters: 
    ------
    molecules (dictionary): 
        Dictionary whichs keys are used to identify the molecules/
        data to list.
    colors (list): 
        List, containing the colors as rgba values to choose from.
    selection_text (string):
        Name of the Area/Box listing the molecules/data to select
        from. 

    Methods:
    ------
    acceptColor():
        Check input for validity, return a tuple containing the 
        molecule identifier and the selected color.
    rejectColor():
        Reset variables (=None), return variables, close the
        dialog.
    changeCheckStateMolecules():
        Update the possibility to select molecules (checkable / 
        not checkable).
    changeCheckStateColors():
        Update the possibility to select colors (checkable / 
        not checkable).

    """

    def __init__(self, molecules, colors, selection_text=''):
        """Create a dialog to change the line color of the molecules
        represented by a line in the main view.
 
        """
        super().__init__()
        self.setWindowTitle('Change Color')
        self.molecule = None
        self.color = None
        try:
            self.caution_text = selection_text.split()[1].lower()
        except:
            self.caution_text = 'a checkbox'

        # DIALOG BUTTONS
        qbtn = QDialogButtonBox.StandardButton.Apply | QDialogButtonBox.StandardButton.Cancel
        self.button_box = QDialogButtonBox(qbtn)
        self.button_box.rejected.connect(self.reject)
        self.button_box.clicked.connect(self.acceptColor)
        self.rejected.connect(self.rejectColor)
        
        # MOLECULES / NAMES CHECKBOXES
        layout = QVBoxLayout()
        molecules_groupbox = QGroupBox(selection_text)
        self.molecules_layout = QGridLayout()
        i = 0
        j = 0
        for molecule in molecules:
            check_btn = QCheckBox(molecule)
            check_btn.clicked.connect(self.changeCheckStateMolecules)
            self.molecules_layout.addWidget(check_btn, i, j) # arrange in two columns
            j += 1
            if j > 1:
                i += 1
                j = 0
        molecules_groupbox.setLayout(self.molecules_layout)
        
        # LABEL
        label = QLabel('Select New Color')

        # COLOR SCROLL AREA
        self.color_layout = QVBoxLayout()
        widget = QWidget()
        scroll_area = QScrollArea()
        i = 0
        for i in colors: # self.parent.data.color_list
            _color_r = int(i[0]*255)
            _color_g = int(i[1]*255)
            _color_b = int(i[2]*255)
            self.color_display = QHBoxLayout()
            self.check_btn = QCheckBox('')
            self.check_btn.setAccessibleName(str(i)) # string of the color-array
            self.check_btn.setCheckable(True)
            self.check_btn.clicked.connect(self.changeCheckStateColors)
            col = Color((_color_r, _color_g, _color_b), 20)
            col.setAccessibleName(str(i))
            self.color_display.addWidget(col)
            self.color_display.addWidget(self.check_btn)
            self.color_layout.addLayout(self.color_display)
        widget.setLayout(self.color_layout)
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(widget)

        # LAYOUT
        layout.addWidget(molecules_groupbox)
        layout.addWidget(label)
        layout.addWidget(scroll_area)
        layout.addWidget(self.button_box)
        self.setLayout(layout)

    def acceptColor(self, s):
        """Return molecule and corresponding new color.
        
        Parameter:
        ------
        s : Reference
            Caller of this method.

        Return:
        ----
        (molecule, color) : tuple
            Tuple containing the selected molecule identifier and 
            selected color.
            
        """
        if s.text() == 'Apply': # both dialog buttons call this function
            if self.molecule is None or self.color is None:
                message_box = QMessageBox(self)
                message_box.setWindowTitle('Caution!')
                message_box.setText('You need to specify a ' + self.caution_text + ' and a color.')
                message_box.setStandardButtons(QMessageBox.StandardButton.Ok)
                self.buttonBox = message_box.exec()
            else:
                self.accept()
                return (self.molecule, self.color)

    def rejectColor(self):
        """Return None for molecule and new color.
        
        Return:
        ----
        (molecule, color) : tuple
            Tuple containing the selected molecule identifier (=None)
            and selected color (=None).

        """
        self.molecule = None
        self.color = None
        return (self.molecule, self.color)

    def changeCheckStateMolecules(self, s):
        """Enable / Disable the checkability of molecule check boxes.
        
        Parameter:
        ------
        s : bool
            True if button is checked, False otherwise.

        """
        if s:
            for index in range(self.molecules_layout.count()):
                widget = self.molecules_layout.itemAt(index).widget()
                if not widget.isChecked():
                    widget.setEnabled(False)
                elif widget.isChecked():
                    self.molecule = widget.text()
        else:
            for index in range(self.molecules_layout.count()):
                widget = self.molecules_layout.itemAt(index).widget()
                widget.setEnabled(True)
            self.molecule = None

    def changeCheckStateColors(self, s):
        """Enable / Disable the checkability of color check boxes.
        
        Parameter:
        ------
        s : bool
            True if color is checked, False otherwise.

        """
        if s:
            for index in range(self.color_layout.count()):
                layout = self.color_layout.itemAt(index)
                for index in range(layout.count()):
                    widget = layout.itemAt(index).widget()
                    if widget.__class__.__name__ == 'QCheckBox':
                        if not widget.isChecked():
                            widget.setEnabled(False)
                        elif widget.isChecked():
                            self.color = widget.accessibleName()
                            for child in self.findChildren(Color):
                                if child.accessibleName == self.color:
                                    self.color_rgba = child.color
        else:
            for index in range(self.color_layout.count()):
                layout = self.color_layout.itemAt(index)
                for index in range(layout.count()):
                    widget = layout.itemAt(index).widget()
                    if widget.__class__.__name__ == 'QCheckBox':
                        widget.setEnabled(True)
            self.color = None


class ChangeNameDialog(QDialog):
    """Show a window, presenting the original names and input fields
    for the user to set a new name.
    
    Parameter: 
    ------
    current (list):
        List of names (str) to be renamed.

    Methods:
    ------
    acceptChanges():
        Store result and close the dialog.
    closeWindow():
        Set the variable containing the new selection to 'None'.
    rejectChanges():
        Reset the variables, closes the dialog.
    """
    def __init__(self, current):
        """Current []: List of strings to be renamed."""
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)
        self.new_names = {}
        scroll_area = QScrollArea()
        self.names_layout = QFormLayout()
        scroll_area.setLayout(self.names_layout)
        for string in current:
            self.new_names[string] = ''
            old = QLabel(string)
            new = QLineEdit(string)
            self.names_layout.addRow(old, new)
        qbtn = QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        self.button_box = QDialogButtonBox(qbtn)
        self.button_box.accepted.connect(self.acceptChanges)
        self.button_box.rejected.connect(self.rejectChanges)
        self.rejected.connect(self.closeWindow)

        layout.addWidget(scroll_area)
        layout.addWidget(self.button_box)

    def acceptChanges(self):
        """Collect the user input and stores it in a dictionary. Closes
        the dialog afterwards.
        
        """
        for i in range(self.names_layout.rowCount()):
            result = self.names_layout.takeRow(0) # takeRow() removes the row
            old = result.labelItem.widget().text()
            new = result.fieldItem.widget().text()
            self.new_names[old] = new
        self.accept()

    def closeWindow(self):
        """Sets the variable new_names to 'None'."""
        self.new_names = None

    def rejectChanges(self):
        """Sets the variable new_names to 'None' and closes the 
        dialog.
        
        """
        self.new_names = None
        self.reject()
        
class SelectRepresentationAtom(QDialog):
    """Open a dialog to select and set the atom for line representation for each molecule.
    
    Parameter:
    ------
    atom_names (list/set): list of atom names for line representation
    selected_atom (str or None): The selected atom name. None if no atom is selected.

    Methods:
    ------
    initUI():
        Initialise UI.
    onItemClicked():
        Slot method called when an item in the list is clicked.
    """
    def __init__(self, atom_names):
        super().__init__()
        self.atom_names = sorted(atom_names)
        self.selected_atom = None
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()
        
        instruction_label = QLabel("Please select an atom (of the selected molecule) for line representation:")
        layout.addWidget(instruction_label)
        
        self.list_widget = QListWidget()
        self.list_widget.addItems(self.atom_names)
        self.list_widget.itemClicked.connect(self.onItemClicked)
        layout.addWidget(self.list_widget)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self.setLayout(layout)
        self.setWindowTitle("Select anker atom")

    def onItemClicked(self, item):
        self.selected_atom = item.text()
