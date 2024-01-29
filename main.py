"""This is the main module. 

The main module contains basic program settings and starts the 
application.
"""

__version__ = '1.0'
__author__ = 'Christiane Rohse'

import sys
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QApplication
from gui import MainWindow

app = QApplication(sys.argv)
# set Icon for all windows
app.setWindowIcon(QIcon('./windowIcon.png'))

w = MainWindow()
w.show()
try:
    sys.exit(app.exec())
except SystemExit:
    print('end of execution')