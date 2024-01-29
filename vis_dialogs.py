"""This is the vis_dialogs module.

The vis_dialogs module contains dialog classes to be used with the
visualisations of the visualisations module. 

Classes:
------
AnnotationDialog:
    Dialog to collect relevant information to add an annotation to
    a visualisation. 
"""

__version__ = '1.0'
__author__ = 'Christiane Rohse'

from PyQt6.QtWidgets import QDialog, QDialogButtonBox, QVBoxLayout, QHBoxLayout


class AnnotationDialog(QDialog): 
    """Show a dialog to add, change and remove annotations from a 
    dictionary of annotations (e.g. used in the below main view plotting).

    Methods:
    ------
    acceptChanges():
        Accept changes and store information in variables, close dialog.
    rejectChanges():
        Store 'None' in all variables, close dialog.
    """

    def __init__(self, title, annotations):
        """Show open dialog with selection possibilities for a topology and a data file."""
        super().__init__()
        title = 'Annotations: '+ title
        self.setWindowTitle(title)
        self.annotations_changed = False

        qbtn = QDialogButtonBox.StandardButton.Open | QDialogButtonBox.StandardButton.Cancel
        self.buttonBox = QDialogButtonBox(qbtn)
        self.buttonBox.accepted.connect(self.acceptChanges)
        self.buttonBox.rejected.connect(self.rejectChanges)

        layout = QVBoxLayout()
        for key, value in annotations.items():
            annotation_layout = QHBoxLayout()
            # TODO: add all labels, buttons, etc.
            layout.addLayout(annotation_layout)

        layout.addWidget(self.buttonBox)
        self.setLayout(layout)

    def acceptChanges(self):
        """"""
        # TODO: implement.

    def rejectChanges(self):
        """"""
        # TODO: implement.
