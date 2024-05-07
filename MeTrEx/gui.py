"""This is the gui module.

The gui module contains all parameters and classes relevant for the
functionalities and visualisation of the gui.


Classes:
------
MainWindow:
    The main window of MeTrEx.
SubWindow:
    Separate window to show any kind of data,
    e.g. external loaded data.
BottomView:
    Object for showing analysed parameters of the data.
BottomViewGroup:
    Container for BottomViews. Extendable for dran-and-drop option.
DocumentationWindow:
    Shows information about MeTrEx.
MainSlider:
    Slider as GUI element.
MenuAction:
    Interface for GUI elements, visualisation, analyses and data
    objects.

Parameter:
----
slider_style:
    Slider properties in CSS format.
"""

__version__ = '1.0'
__author__ = 'Christiane Rohse'

from cmath import sqrt
from copy import deepcopy
import os
import platform #FIX by Beat to check for MacOS
import psutil
import random
import csv
from PyQt6.QtCore import Qt, pyqtSignal, pyqtSlot, QCoreApplication, QSize
from PyQt6.QtGui import QAction, QKeySequence, QFont, QGuiApplication
from PyQt6.QtWidgets import QMainWindow, QLabel, QCheckBox, QPushButton,\
    QVBoxLayout, QHBoxLayout, QFileDialog, QMessageBox, QGridLayout,\
    QWidget, QSlider, QGroupBox, QStyle, QScrollArea, QSpinBox,\
    QSizePolicy, QFrame
from data import Data, DataSub, SingleDataSet
import analyses
import colormap
from colormap import returnColors
import numpy as np
import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import MDAnalysis as mda
from visualisations import MainPlotCanvas, SubPlotCanvas
from dialogs import OpenDialog, ChangeViewDataDialog,\
    PreprocessingSelectionDialog, ChangeColormapDialog, OpenXVGDialog,\
    AboutDialog, ChangeNameDialog, AnalysisSelectionDialog, ChangeColorDialog, \
    DocumentationDialog


slider_style = """
QSlider {{
    margin:{_margin}px;
}}
QSlider::groove:horizontal {{
    border-radius: {_bg_radius}px;
    /*border: 0px solid #424242; */
    height: {_bg_size}px;
    margin: 0px;
    background-color: {_bg_color};
}}
QSlider::grove:horizontal:hover{{
    background-color: {_bg_color_hover};
}}
QSlider::handle:horizontal {{
    border: none;
    height: {_handle_size}px;
    width: {_handle_size}px;
    margin: {_handle_margin}px;
    border-radius: {_handle_radius}px;
    background-color: {_handle_color};
}}
QSlider::handle:horizontal:hover {{
    background-color: {_handle_color_hover};
}}
QSlider::handle:horizontal:pressed {{
    background-color: {_handle_color_pressed};
}}
QSlider::add-page:horizontal {{
    background: thistle;
}}
QSlider::sub-page:horizontal {{
    background: purple;
}}
"""


class MainSlider(QSlider):
    """Create a horizontal slider with custom settings."""

    def __init__(
            self,
            margin=0,
            bg_size=10,
            bg_radius=4,
            bg_color='red',
            bg_color_hover='blue',
            handle_margin=0,
            handle_size=11,
            handle_radius=5,
            handle_color='fuchsia',
            handle_color_hover='black',
            handle_color_pressed='red'
            ):
        """Create a QSlider with a customized style sheet."""
        super(MainSlider, self).__init__()
        self.setOrientation(Qt.Orientation.Horizontal)
        adjust_style = slider_style.format(
            _margin=margin,
            _bg_size=bg_size,
            _bg_radius=bg_radius,
            _bg_color=bg_color,
            _bg_color_hover=bg_color_hover,
            _handle_margin=handle_margin,
            _handle_size=handle_size,
            _handle_radius=handle_radius,
            _handle_color=handle_color,
            _handle_color_hover=handle_color_hover,
            _handle_color_pressed=handle_color_pressed
            )
        self.setStyleSheet(adjust_style)


class DocumentationWindow(QWidget):
    """Show information about this program.
    
    Parameter:
    ------
    parent (QMainWindow):
        Parent of this widget.
    """
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.setWindowTitle('MeTrEx')
        self.setMinimumSize(340, 300) # w, h

        info = QLabel('Version 1.0')
        layout = QVBoxLayout()
        layout.addWidget(info)
        self.setLayout(layout)


class SubWindow(QWidget):
    """SubWindow class shows given data and emits signals when values
    change.

    Parameters:
    -------
    sliderChanged(pyqtSignal):
        Accepted signal for updating the slider position.

    Methods:
    -------
    showAnimation():
        Show an animation of the data in the sub window.
    linkSliders():
        Update variable to show if sliders are linked.
    updateLinkedSliders():
        Update the slider position and the view.
    jumpToFrame():
        Update slider position, view and information panelc according to
        user specification.
    sliderChanged():
        Update view, information panel and jump_to spin box while the
        slider is moved.
    sliderChangedFinished():
        Update view, information panel and jump_to spin box after the
        slider has been moved.
    setInformation():
        Populate the information panel.
    setSetting():
        Populate the settings panel.
    changeColor():
        Change the color of a shown data set.
    changename():
        Change the name of a shown data set.
    showHideSphere():
        Show or hide the sphere in the view.
    saveView():
        Save an .png image file of the view, including the legend.

    """
    sliderChanged = pyqtSignal(int)

    def __init__(self, parent, **kwargs):
        super().__init__()
        # ADJUST APPEARANCE
        self.parent = parent
        self.setWindowTitle('MeTrEx')
        self.setMinimumSize(650, 520) # w, h
        self.data = None
        self.large_data = False
        self.slider_position = 1
        self.sliders_linked = False
        self.show_spheres = True
        self.show_animation = False
        self.colormap = 'plasma'

        for key, value in kwargs.items():
            if key == 'title':
                self.setWindowTitle(value)
            if key == 'data':
                self.data = value
            if key == 'large':
                self.large_data = value

        layout = QGridLayout()

        # PLOTTING
        plotting_layout = QVBoxLayout()
        self.sc = SubPlotCanvas(self, width=7.5, height=10, dpi=100)
        plotting_layout.addWidget(self.sc)

        # SLIDER
        self.slider = MainSlider()
        self.slider.setAccessibleName('sub_slider')
        self.slider.setMinimum(self.data.slider_range[0])
        self.slider.setMaximum(self.data.slider_range[1])
        self.slider.setSingleStep(1)
        if self.sc.dim == '3d' or self.data.slider_range[1] > 50000:
            self.slider.sliderReleased.connect(self.sliderChangedFinished)
        else:
            self.slider.valueChanged[int].connect(self.sliderChanged)
        control_slider_layout =QVBoxLayout()
        from_to_layout = QHBoxLayout()
        from_to_font = QFont()
        from_to_font.setBold(True)
        from_to_font.setPixelSize(20)
        slider_from = QLabel(str(self.data.slider_range[0]))
        slider_from.setFont(from_to_font)
        self.slider_to = QLabel(str(self.data.slider_range[1]))
        self.slider_to.setFont(from_to_font)
        from_to_layout.addWidget(slider_from, 0, Qt.AlignmentFlag.AlignLeft)
        from_to_layout.addWidget(self.slider_to, 0, Qt.AlignmentFlag.AlignRight)
        
        # ANIMATION & LINK SLIDERS
        self.animation_layout = QGridLayout()
        self.animation_button = QPushButton('')
        pixmap_media_play = QStyle.StandardPixmap.SP_MediaPlay
        animation_icon = self.parent.style().standardIcon(pixmap_media_play)
        self.animation_button.setIcon(animation_icon)
        self.animation_button.setCheckable(True)
        self.animation_button.clicked.connect(self.showAnimation)
#        self.link_sliders = QCheckBox('link sliders')
#        self.link_sliders.setToolTip('Link this slider to the main view slider\nif frames match x-axis data.')
#        self.link_sliders.setChecked(False)
#        self.link_sliders.stateChanged.connect(self.linkSliders)
#        self.parent.sliderChanged.connect(self.updateLinkedSliders)
        jump_label = QLabel('jump to frame number:')
        self.jump_to = QSpinBox()
        self.jump_to.setRange(self.data.slider_range[0], self.data.slider_range[1])
        self.jump_to.setAccessibleName('jump_to')
        self.jump_to.editingFinished.connect(self.jumpToFrame)
        
        control_slider_layout.addLayout(from_to_layout)
        control_slider_layout.addWidget(self.animation_button, 0, Qt.AlignmentFlag.AlignLeft)
#        control_slider_layout.addWidget(self.link_sliders, 0, Qt.AlignmentFlag.AlignRight)
        control_slider_layout.addWidget(jump_label, 0, Qt.AlignmentFlag.AlignRight)
        control_slider_layout.addWidget(self.jump_to, 0, Qt.AlignmentFlag.AlignRight)
        plotting_layout.addWidget(self.slider, 1)
        plotting_layout.addLayout(control_slider_layout, 1)

        # INFORMATION GRID
        self.info_grid_layout = QVBoxLayout()
        # STATIC AND DYNAMIC INFO
        self.controls_layout = QVBoxLayout()
        self.setInformation()
        # SETTINGS
        self.settings_layout = QVBoxLayout()
        self.setSettings()
        self.info_grid_layout.addLayout(self.controls_layout)
        self.info_grid_layout.addLayout(self.settings_layout)
        self.info_grid_layout.addStretch()

        # LAYOUT SETUP
        layout.addLayout(plotting_layout, 0, 0)
        layout.addLayout(self.info_grid_layout, 0, 1)
        self.setLayout(layout)

    def showAnimation(self, s):
        """Play animation."""
        start = self.slider_position
        end = self.slider.maximum()
        self.show_animation = s
        pixmap_media_pause = QStyle.StandardPixmap.SP_MediaPause
        animation_icon_pause = self.parent.style().standardIcon(pixmap_media_pause)
        QCoreApplication.processEvents() #
        if s:
            self.animation_button.setIcon(animation_icon_pause)
        while self.show_animation:
            self.sliderChanged(start)
            QCoreApplication.processEvents()
            start += 1
            if start > end and s:
                start = 0
            self.slider.setValue(start)
        pixmap_media_play = QStyle.StandardPixmap.SP_MediaPlay
        animation_icon_play = self.style().standardIcon(pixmap_media_play)
        self.animation_button.setIcon(animation_icon_play)
        self.animation_button.setChecked(False)

#    def linkSliders(self, link):
#        """Update link_sliders variable according to selection."""
#        if link == 0: # sliders not linked
#            self.sliders_linked = False
#            return
#        if link == 2: # sliders to be linked
#            self.sliders_linked = True
#            return
#
#    @pyqtSlot(int)
#    def updateLinkedSliders(self, pos):
#        """Update the slider position of the view according to the
#        signal from the main window.
#
#        Parameter:
#        -----
#        position (int):
#            Position of the slider.
#
#        """
#        if self.sliders_linked:
#            if pos > self.data.slider_range[1]:
#                pos = self.data.slider_range[1]
#            self.sliderChanged(pos)
#            self.slider.setValue(pos)
#        else:
#            return

    def jumpToFrame(self):
        """Update the shown data view, info panel and slider according
        to the user selection.
        """
        new_pos = self.sender().value()
        self.slider.setValue(new_pos)
        self.sliderChanged(new_pos)

    def sliderChanged(self, s):
        """Updates the shown data view, info panel and jump_to QSpinBox
        according to the sliders current position.

        Parameter:
        -----
        position (int):
            Position of the slider.

        Note:
        ----
        Not implemented: Bidirectional linkage of sliders.
        """
        self.sc.updateSpheresScatterPlot(s-1, self.show_spheres)
        self.slider_position = s
        self.slider_position_label.setText(str(s))
        self.sc.draw()
        self.jump_to.setValue(s)

        # bidirectional linkage
        if self.sliders_linked:
            # Here you can implement bidirectional linkage of sliders.
            pass

    def sliderChangedFinished(self):
        """Update the shown data view, info panel and jump_to QSpinBox
        according to the sliders position when it was released.
       
        This function is particularly useful for large data, when a
        continuous update of the view becomes too slow.
        """
        s = self.sender().value()
        self.sliderChanged(s)

    def setInformation(self):
        """Setup information panel and show relevant information of the
        data.
        """
        self.overview_data_groupbox = QGroupBox('General Information')
        self.overview_data_groupbox.setAccessibleName('information')
        self.overview_data_layout = QVBoxLayout()
        self.overview_data_groupbox.setMinimumWidth(180)

        font = QFont()
        font.setBold(True)
        font.setPixelSize(30)
        self.slider_position_label = QLabel(str(self.slider_position))
        self.slider_position_label.setFont(font)
        self.slider_position_label.setStyleSheet('color: fuchsia')
        self.slider_position_label.setToolTip('Current position of the slider.\n(Line index of data lines)')
        
        lines_shown = QLabel(str(self.data.slider_range[1]) + ' lines')
        title_shown = QLabel(self.data.analysis)

        self.overview_data_layout.addWidget(self.slider_position_label, 0, Qt.AlignmentFlag.AlignRight)
        self.overview_data_layout.addWidget(lines_shown, 0, Qt.AlignmentFlag.AlignRight)
        self.overview_data_layout.addWidget(title_shown, 0, Qt.AlignmentFlag.AlignRight)

        self.overview_data_groupbox.setLayout(self.overview_data_layout)
        self.controls_layout.addWidget(self.overview_data_groupbox)

    def setSettings(self):
        """Setup settings grid for this view.
        This includes: change color, change name, save and hide sphere
        buttons.
        """
        self.settings_groupbox = QGroupBox('Settings & Options')
        self.settings_groupbox.setAccessibleName('settings')
        self.settings_groupbox_layout = QVBoxLayout()
        self.settings_groupbox.setMinimumWidth(180)
        # Settings
        change_color = QPushButton('Color')
        change_color.clicked.connect(self.changeColor)
        change_color.setToolTip('Change the colors of the shown data sets.')
        change_name = QPushButton('Name')
        change_name.clicked.connect(self.changeName)
        change_name.setToolTip('Change the names of the shown data sets.')
        save = QPushButton('Save')
        save.clicked.connect(self.saveView)
        self.show_hide_sphere = QPushButton('Hide Sphere')
        self.show_hide_sphere.setCheckable(True)
        self.show_hide_sphere.clicked.connect(self.showHideSphere)
        # Add to layout
        self.settings_groupbox_layout.addWidget(change_color)
        self.settings_groupbox_layout.addWidget(change_name)
        self.settings_groupbox_layout.addWidget(save)
        self.settings_groupbox_layout.addWidget(self.show_hide_sphere)
        self.settings_groupbox.setLayout(self.settings_groupbox_layout)
        self.settings_layout.addWidget(self.settings_groupbox)

    def changeColor(self):
        """Change the color of a single data set."""
        current_names = {}
        for name, ps in self.data.pointsets.items():
            current_names[ps.selection_string] = ''
        n = 30
        dlg = ChangeColorDialog(current_names, returnColors(n, self.colormap), 'Select Dataset')
        dlg.exec()
        if dlg.molecule is None or dlg.color is None:
            return
        else:
            color_string = list(dlg.color[1:-1].split(' '))
            color_string = list(filter(lambda a: a != '', color_string))
            new_color = np.array(list(map(float, color_string)))
            for plot, reference in self.sc.pointsets.items():
                if plot == dlg.molecule:
                    reference.set(color=new_color)
                    self.sc.data.pointsets[plot].color = new_color
            self.sc.ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
            self.sc.draw()

    def changeName(self):
        """Change the name(s) of the shown data set(s).

        This includes opening a selection dialog, changing the data
        object if needed and updating the view.
        """
        current_names = self.sc.pointsets.keys()
        dlg = ChangeNameDialog(current_names)
        dlg.exec()
        if dlg.new_names is None:
            return
        else:
            for on, nn in dlg.new_names.items():
                self.sc.data.pointsets[nn] = self.sc.data.pointsets.pop(on)
                self.sc.data.pointsets[nn].selection_string = nn
            self.sc.updatePointsets()

    def showHideSphere(self, s):
        """Change the sphere visibility in the data view according to
        the user selection.

        Parameter:
        ------
        s (bool):
            True: hide sphere(s), False: show sphere(s).
        """
        if s:
            self.show_hide_sphere.setText('Show Sphere')
            for sphere, reference in self.sc.spheres.items():
                reference.remove()
            self.sc.draw()
            self.show_spheres = False
        else:
            self.show_hide_sphere.setText('Hide Sphere')
            self.show_spheres = True
            if self.sc.dim == '3d':
                for name, ds in self.sc.data.pointsets.items():
                    self.sc.spheres[ds.selection_string] = self.sc.ax.scatter(
                        ds.x[self.slider_position-1],
                        ds.y[self.slider_position-1],
                        ds.z[self.slider_position-1],
                        color='fuchsia',
                        s=80
                        )
            elif self.sc.dim == '2d':
                for name, ds in self.sc.data.pointsets.items():
                    self.sc.spheres[ds.selection_string] = self.sc.ax.scatter(
                        ds.x[self.slider_position-1],
                        ds.y[self.slider_position-1],
                        color='fuchsia',
                        s=80
                        )
            else:
                return
            self.sc.draw()

    def saveView(self):
        """Save the view of this widget to a png file."""
        file_filter = 'Image File (*.png)'
        response = QFileDialog.getSaveFileName(
            parent=self,
            caption='image.png',
            directory=os.getcwd(),
            filter=file_filter,
            initialFilter='Image File (*.png)'
            )
        image_path = response[0]
        self.sc.fig.savefig(image_path)


class BottomView(QWidget):
    """Visualize analysis data.

    Methods
    ---
    updateInfoPanel(int):
        Update the information panel sfter the slide in the view
        has changed.
    setInfoPanel():
        
    saveView()
    updateVLine()
    changeColor()
    """
    def __init__(self, parent, name, color, data):
        super(BottomView, self).__init__()
        self.parent = parent
        self.setContentsMargins(0,0,0,0)
        self.layout = QHBoxLayout()
        self.layout.setObjectName(name)
        self.layout.setContentsMargins(0,0,0,0)
        self.name = name
        self.data = data
        r = color[0]*255
        g = color[1]*255
        b = color[2]*255

        self.id_label = QLabel('') # because you cannot write the name vertical easily, don't use any name here
        self.id_label.setStyleSheet('background-color: rgb('+str(r)+','+str(g)+','+str(b)+')')
        self.id_label.setAccessibleName(self.name)
        self.id_label.setMaximumWidth(20)
        self.id_label.setMinimumHeight(180)
        
        self.view = SubPlotCanvas(self, 7.5, 8, 80)
        self.view.setAccessibleName(self.name)
        self.view.vline_changed.connect(self.updateInfoPanel)

        info_panel_name = self.data.analysis + ' ' + self.data.unit
        self.information_panel_groupbox = QGroupBox(info_panel_name) # self.data.analysis
        self.information_panel_groupbox.setAccessibleName(self.name)
#        self.information_panel_groupbox.setFixedWidth(250)
        self.information_panel_groupbox.setMinimumWidth(250)
        self.information_panel_groupbox.setMaximumWidth(300)

        self.control_panel_groupbox = QGroupBox('Controls')
        self.control_panel_groupbox.setAccessibleName(self.name)
        self.control_panel_groupbox.setFixedWidth(70)
        controls_layout = QVBoxLayout()
        jump_to = QSpinBox()
        jump_to.setAccessibleName('jump_to')
        jump_to.setRange(1, self.data.slider_range[1])
        jump_to.setValue(1)
        jump_to.setSingleStep(10)
        jump_to.editingFinished.connect(self.updateVLine)
        show_hide_labels = QPushButton('s/h')
        show_hide_labels.setCheckable(True)
        show_hide_labels.setToolTip('Show or hide the minimum and maximum labels.')
        show_hide_labels.clicked.connect(self.showHideLabels)
        print_button = QPushButton('save')
        print_button.setAccessibleName(self.name)
        print_button.clicked.connect(self.saveView)
        remove_button = QPushButton('-')
        remove_button.setMaximumWidth(25)
        remove_button.setAccessibleName(self.name)
        remove_button.setCheckable(True)
        remove_button.clicked.connect(self.parent.deleteNewView) #self.parent.deleteView
        controls_layout.addWidget(jump_to)
        if len(self.data.lines) == 1:
            controls_layout.addWidget(show_hide_labels)
        controls_layout.addWidget(print_button)
        controls_layout.addWidget(remove_button)
        controls_layout.addStretch(1)
        self.control_panel_groupbox.setLayout(controls_layout)

        self.layout.addWidget(self.id_label)
        self.layout.addWidget(self.view)
        self.layout.addWidget(self.information_panel_groupbox)
        self.layout.addWidget(self.control_panel_groupbox)
        self.setLayout(self.layout)

    @pyqtSlot(int)
    def updateInfoPanel(self, pos):
        """Update the info panel when slider in view changes."""
        # update grid
        if pos > self.data.slider_range[1]:
            for element in self.info_panel_current_values:
                element[0].setText(str('{:.2f}'.format(element[1][self.data.slider_range[1]-1])))
                element[2].setText(str(self.data.slider_range[1]))
        else:
            for element in self.info_panel_current_values:
                element[0].setText(str('{:.2f}'.format(element[1][pos-1])))
                element[2].setText(str(pos))
        # update jump to spin box
        for child in self.findChildren(QSpinBox):
            if child.accessibleName() == 'jump_to':
                child.setValue(pos)

    def setInfoPanel(self):
        """Populate the information panel."""
        information_panel_layout = QVBoxLayout()
        # CHECKBOX TO SAVE CSV
        self.save_values = QCheckBox('save static values')
        self.save_values.setChecked(True)
        self.save_values.setLayoutDirection(Qt.LayoutDirection.RightToLeft)
        information_panel_layout.addWidget(self.save_values, 0, Qt.AlignmentFlag.AlignLeft.AlignTop)

        # Font to highlight grid values
        font_highlight_grid = QFont()
        font_highlight_grid.setBold(True)
        # PUT ALL VALUES IN A SCROLL AREA
        values_grid_wrapper = QWidget()
        values_grid_scroll = QScrollArea()
        values_grid_scroll.setFrameShape(QFrame.Shape.NoFrame)
        values_grid = QGridLayout()
        values_grid.setContentsMargins(0,0,0,0)
        name_label = QLabel('')
        font_curr = QFont()
        font_curr.setBold(True)
        current_label = QLabel('curr')
        current_label.setFont(font_curr)
        current_label.setToolTip('Value at slider position.')
        max_label = QLabel('max')
        min_label = QLabel('min')
        avg_label = QLabel('avg')
        values_grid.addWidget(name_label, 0, 0, Qt.AlignmentFlag.AlignLeft)
        values_grid.addWidget(current_label, 0, 1, Qt.AlignmentFlag.AlignRight)
        values_grid.addWidget(max_label, 0, 2, Qt.AlignmentFlag.AlignRight)
        values_grid.addWidget(min_label, 0, 3, Qt.AlignmentFlag.AlignRight)
        values_grid.addWidget(avg_label, 0, 4, Qt.AlignmentFlag.AlignRight)
        # Add a line for each line shown in the plot
        i = 1
        self.info_panel_current_values = []
        for line in self.data.lines:
            name_line = QLabel(line.selection_string)
            name_line.setFont(font_highlight_grid)
            curr_val = QLabel(str('{:.2f}'.format(line.y[0]))) # vline starts at index 0
            curr_val.setAccessibleName(line.selection_string)
            curr_val.setFont(font_highlight_grid)
            max_val = QLabel(str('{:.2f}'.format(line.y_max)))
            min_val = QLabel(str('{:.2f}'.format(line.y_min)))
            avg_val = QLabel(str('{:.2f}'.format(line.y_avg)))
            values_grid.addWidget(name_line, i, 0, Qt.AlignmentFlag.AlignLeft)
            values_grid.addWidget(curr_val, i, 1, Qt.AlignmentFlag.AlignRight)
            values_grid.addWidget(max_val, i, 2, Qt.AlignmentFlag.AlignRight)
            values_grid.addWidget(min_val, i, 3, Qt.AlignmentFlag.AlignRight)
            values_grid.addWidget(avg_val, i, 4, Qt.AlignmentFlag.AlignRight)
            i += 1
            frame_no = QLabel('Frame:')
            curr_frame = QLabel('1') # vline starts at pos 1
            curr_frame.setAccessibleName(line.selection_string)
            # Frame indices are 1, ..., n+1
            if len(line.y_maxi) > 1:
                max_frame = QLabel(str(line.y_maxi[0]+1)+',...')
            else:
                max_frame = QLabel(str(line.y_maxi[0]+1))
            if len(line.y_mini) > 1:
                min_frame = QLabel(str(line.y_mini[0]+1)+',...')
            else:
                min_frame = QLabel(str(line.y_mini[0]+1))
            values_grid.addWidget(frame_no, i, 0, Qt.AlignmentFlag.AlignLeft)
            values_grid.addWidget(curr_frame, i, 1, Qt.AlignmentFlag.AlignRight)
            values_grid.addWidget(max_frame, i, 2, Qt.AlignmentFlag.AlignRight)
            values_grid.addWidget(min_frame, i, 3, Qt.AlignmentFlag.AlignRight)
            self.info_panel_current_values.append([curr_val, line.y, curr_frame]) # this makes updating more easy
            i += 1
        # assemble layout
        values_grid_wrapper.setLayout(values_grid)
        values_grid_scroll.setWidget(values_grid_wrapper)
        information_panel_layout.addWidget(values_grid_scroll)
        self.information_panel_groupbox.setLayout(information_panel_layout)

    def showHideLabels(self, s):
        """Show or hide the sphere.
        
        Parameter:
        -----
        s (bool):
            When checked: True = hide sphere, otherwise show sphere.
        """
        if s:
            self.view.label_min.remove()
            self.view.label_max.remove()
            self.view.draw()
        else:
            self.view.addMinMaxLabels2DLine()
            self.view.draw()

    def saveView(self):
        """Save a view as a file type supported my Matplotlib."""
        # for some positions the index need to be adjusted (display: 1, ..., n+1)
        file_filter = 'Image File (*.png)'
        response = QFileDialog.getSaveFileName(
            parent=self,
            caption='image.png',
            directory=os.getcwd(),
            filter=file_filter,
            initialFilter='Image File (*.png)'
            )
        image_path = response[0]
        self.view.fig.savefig(image_path)
        if self.save_values.isChecked():
            # get information to be saved:
            data = []
            analysis_unit = [self.data.analysis + ' ' + self.data.unit]
            labels = ['name', 'value_current', 'frame_current', 'value_min', 'frame_min', 'value_max', 'frame_max', 'value_avg']
            data.append(analysis_unit)
            data.append(labels)
            for line in self.data.lines:
                kids = self.findChildren(QSpinBox)
                for kid in kids:
                    if kid.accessibleName() == 'jump_to':
                        curr = kid.value()
                name = line.selection_string
                curr_val = str('{:.2f}'.format(line.y[curr-1]))
                curr_frame = str(curr)
                max_val = str('{:.2f}'.format(line.y_max))
                max_frame = ''
                for i in range(len(line.y_maxi)):
                    max_frame += str(line.y_maxi[i] + 1) + ' '
                min_val = str('{:.2f}'.format(line.y_min))
                min_frame = ''
                for j in range(len(line.y_mini)):
                    min_frame += str(line.y_mini[j] + 1) + ' '
                avg_val = str('{:.2f}'.format(line.y_avg))
                data.append([name, curr_val, curr_frame, min_val, min_frame, max_val, max_frame, avg_val])
            # write to csv file
            csv_path = image_path[:-3] + 'csv'
            with open(csv_path, 'w', newline='') as csv_file:
                writer = csv.writer(csv_file, delimiter=',')
                for line in data:
                    writer.writerow(line)

    def updateVLine(self):
        """Update the vertical line in the view (slider, vline) and
        information panel.
        """
        sender = self.sender()
        value = sender.value()
        if value > self.data.slider_range[1]:
            # if slider out of range: use maximal slider value
            value = self.data.slider_range[1]
        self.view.vline.remove()
        self.view.vline = self.view.ax.axvline(value, c='fuchsia', linewidth=2.5)
        if self.view.single:
            label_text = ' frame: '+ str(value) + '\n' + ' '+str('{:.2f}'.format(self.view.data.lines[0].y[value-1])) # TODO: test index
        else:
            label_text = ' frame: '+ str(value)
        self.view.vline_info.remove()
        x_pos = float(value / self.view.ax.get_xbound()[1])
        self.view.vline_info = self.view.ax.text(x_pos, 1.05, label_text, transform=self.view.ax.transAxes)
        self.view.draw()
        # update table:
        self.updateInfoPanel(value)

    def changeColor(self, d): # dict = {line:color}
        """Change the color of a data set in the view."""
        for line, color in d.items():
            try:
                self.view.lines[line].set_color(color)
                self.view.ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
                self.view.draw()
                if self.id_label.accessibleName() == line:
                    r = color[0]*255
                    g = color[1]*255
                    b = color[2]*255
                    self.id_label.setStyleSheet('background-color: rgb('+str(r)+','+str(g)+','+str(b)+')')
            except: # no such line in plot
                pass
        

class BottomViewGroup(QScrollArea):
    """Holds all BottomView objects.
    
    Parameters
    ------
    parent: parent object, MainView

    Methods:
    ------
    newView(BV):
        Add a new BottomView (BV) to the area below the main view.
    deleteNewView():
        Delete a BottomView.
    """

    def __init__(self, parent):
        super(BottomViewGroup, self).__init__()
        self.parent = parent
        widget = QWidget()
        self.layout = QVBoxLayout()
        widget.setLayout(self.layout)
        self.setMinimumSize(720, 200) # >= 7.5 x 1 in
        self.setMaximumHeight(400)
        self.setWidgetResizable(True)
        self.setAcceptDrops(True)
        self.setWidget(widget)
        self.setContentsMargins(0,0,0,0)

        self.bottom_views = {}

    def newView(self, view):
        """Adds a BottomView object to the layout, registers it in the
        list of views and updates size of BottomViewGroup if neccessary.
        """
        self.layout.addWidget(view)
        self.bottom_views[view.name] = view
        i = len(self.bottom_views)
        if i < 3:
            self.setMinimumSize(720, 200*i)

    def deleteNewView(self):
        sender = self.sender()
        name = sender.accessibleName()
        reply = QMessageBox.question(
            self,
            'Remove analysis',
            'Are you sure you want to remove this analysis?',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes
            )
        if reply == QMessageBox.StandardButton.Yes:
            reply = True
        else:
            reply = False
        if reply:
            self.bottom_views[name].setParent(None)
            del self.bottom_views[name]
            i = len(self.bottom_views)
            if i < 3:
                if i == 0:
                    for child in self.parent.findChildren(BottomViewGroup):
                        child.setHidden(True)
                else:
                    self.setMinimumSize(720, 200*i)
        else:
            sender.setChecked(False)

 
class MainWindow(QMainWindow):
    """Main window class sets up view and controls of the UI and
    connects them to their specific functionality.

    Parameter:
    ------
    sliderChanged (pyqtSignal):
        Signal emitted when the sliders are changed.

    Methods:
    ------
    linkSliders():
        Link both sliders of the main view and, if applicable, the
        BottomView sliders. They will spdate together.
    showAnimation():
        Show an animation of the data in the main view.
    setOverviewLayout():
        Populate the information and settings panel at the right side
        of the main view.
    updateCurrentData():
        Update the values shown in the information panel.
    updateView():
        Show the visualisation of the data in the main view.
    controlButtons():
        Populate buttons for the interaction with the main view.
    changeColorMap():
        Change the colour map for all elements of the GUI and view.
    updatePolynom():
        Recalculate the surface abstraction.
    changeMembraneRepresentation():
        Change the membrane representation from a simple plain to
        a triangulated surface and back.
    setFrame():
        Change the slider position and update the view.
    sliderChangedMolecules():
        Update the view according to slider position.
    updateBottomViewSliders():
        Update the position of the sliders in the BottomViews.
    sliderChangedMembrane():
        Update the membrane state in the main view according to the
        slider position.
    showHideMolecule():
        Show or hide a trajectory in the main view.
    showHideMembraneLeaflet()
        Show or hide a membrane leaflet (upper, lower).
    changeColor():
        Change the color of a trajectory line in the main view.
        If applicable: also change the color in BottomViews.
    clearControls():
        Clear all panels next to the main view.
    redrawMainView():
        Show the main view as shown for the first time.
    changeColorBottomViews():
        Update the color code in the BottomViews.
    redraw():
        Update shperes and leaflets in the main view.
    redrawShowHideMapping():
        Update elements of the main view when a mapping is shown.
    redrawShowHide():
        Update the view when a trajectory shall be shown / hidden.
    redrawColor():
        Update color of trajectory line.
    addSubWindow():
        Show a new Window.
    closeEvent():
        Close all windows shown when the application is exited.
    """
    sliderChanged = pyqtSignal(int)

    def __init__(self):
        super(MainWindow, self).__init__()
        # ADJUST APPEARANCE
        self.setWindowTitle('MeTrEx')
        
        # Get the screen geometry
        screen_geometry = QGuiApplication.primaryScreen().geometry()
        
        # Calculate the minimum size based on a percentage of the screen size
        width_ratio = 0.6
        height_ratio = 0.75
        initial_width = screen_geometry.width() * width_ratio
        initial_height = screen_geometry.height() * height_ratio
        
        # Set the minimum size of the main window
        min_width = 725
        min_height = 810
        self.setMinimumSize(min_width, min_height)# w, h
        
        # Set the initial size of the main window (larger than the minimum size)
        if min_width < initial_width and min_height < initial_height:
            self.resize(initial_width, initial_height)# w, h
        size = self.size()
        self.width = size.width()  # Get width
        self.height = size.height()  # Get height

        # maximal width of the info panel of the main view
        self.size_info_max = 260
        self.size_analysis_unit_label = 90

        # WINDOW VARIABLES
        self.data = Data()
        self.slider_position_molecules = 0
        self.slider_position_membrane = 0
        self.sliders_linked = True
        self.show_animation = False
        self.show_selection = False
        self.selection_from = None
        self.selection_to = None
        self.color_changed_drug_molecules = False
        self.bottom_views = {}
        self.bottom_views_scales = {}
        # register all sub windows here to close them if the application is closed
        self.sub_windows = []

        # SET MAIN LAYOUT, ADD TO CENTRAL WIDGET
        self.layout = QGridLayout()
        widget = QWidget()
        widget.setLayout(self.layout)
        self.setCentralWidget(widget)

        # SLIDER LINES
        self.slider_molecules = MainSlider()
        self.slider_molecules.setAccessibleName('molecules_slider')
        self.slider_molecules.setMinimum(self.data.min_frame)
        self.slider_molecules.setMaximum(self.data.max_frame)
        self.slider_molecules.setHidden(True)
        self.slider_molecules.setSingleStep(1)
        self.slider_molecules.valueChanged[int].connect(self.sliderChangedMolecules)
        self.slider_molecules.sliderReleased.connect(self.updateBottomViewSliders)
        # SLIDER MEMBRANE
        self.slider_membrane = MainSlider()
        self.slider_membrane.setAccessibleName('membrane_slider')
        self.slider_membrane.setMinimum(self.data.min_frame)
        self.slider_membrane.setMaximum(self.data.max_frame)
        self.slider_membrane.setHidden(True)
        self.slider_membrane.setSingleStep(1)
        self.slider_membrane.valueChanged[int].connect(self.sliderChangedMembrane)
        self.slider_membrane.sliderReleased.connect(self.updateBottomViewSliders)
        # SLIDER FROM TO (no ticks possible with stylesheet)
        font_slider = QFont()
        font_slider.setBold(True)
        font_slider.setPixelSize(20)
        self.slider_from_to_layout = QHBoxLayout()
        self.slider_from = QLabel()
        self.slider_from.setHidden(True)
        self.slider_from.setFont(font_slider)
        self.slider_to = QLabel()
        self.slider_to.setHidden(True)
        self.slider_to.setFont(font_slider)
        self.slider_description = QLabel('top slider: line representation molecules   bottom slider: membrane')
        self.slider_description.font().setPointSize(9)
        self.slider_description.setHidden(True)
        self.slider_from_to_layout.addWidget(self.slider_from, alignment=Qt.AlignmentFlag.AlignLeft)
        self.slider_from_to_layout.addWidget(self.slider_description, alignment=Qt.AlignmentFlag.AlignHCenter)
        self.slider_from_to_layout.addWidget(self.slider_to, alignment=Qt.AlignmentFlag.AlignRight)
        # ANIMATION & LINK SLIDERS
        self.animation_layout = QGridLayout()
        self.animation_button = QPushButton('')
        pixmap_media_play = QStyle.StandardPixmap.SP_MediaPlay
        animation_icon = self.style().standardIcon(pixmap_media_play)
        self.animation_button.setIcon(animation_icon)
        self.animation_button.setCheckable(True)
        self.animation_button.clicked.connect(self.showAnimation)
        self.connect_sliders = QCheckBox('link sliders')
        self.connect_sliders.clicked.connect(self.linkSliders)
        self.connect_sliders.setChecked(True)
        self.animation_button.setHidden(True)
        self.connect_sliders.setHidden(True)
        self.animation_layout.addWidget(self.animation_button, 0, 0, Qt.AlignmentFlag.AlignLeft)
        self.animation_layout.addWidget(self.connect_sliders, 0, 0, Qt.AlignmentFlag.AlignRight)

        # JUMP TO FRAME
        self.set_frame_layout = QVBoxLayout()
        self.set_frame_label = QLabel('jump to frame number:')
        self.set_frame_label.setHidden(True)
        self.set_frame = QSpinBox()
        self.set_frame.setAccessibleName('set_frame')
        self.set_frame.setHidden(True)
        self.set_frame.setSingleStep(1)
        self.set_frame.editingFinished.connect(self.setFrame)
        self.set_frame_layout.addWidget(self.set_frame_label, 0, Qt.AlignmentFlag.AlignRight)
        self.set_frame_layout.addWidget(self.set_frame, 0, Qt.AlignmentFlag.AlignRight)

        # CONTROLS LAYOUT
        self.controls_layout = QVBoxLayout()
        self.layout.addLayout(self.controls_layout, 0, 1)

        # MAIN VIEW
        layout_main_view = QVBoxLayout()
        self.sc = MainPlotCanvas(width=7.5, height=8, dpi=100)
        self.toolbar = NavigationToolbar2QT(self.sc, self)
        unwanted_buttons = ['configure_subplots', 'edit_parameters']
        for action in self.toolbar._actions:
            if action in unwanted_buttons:
                self.toolbar.removeAction(self.toolbar._actions[action])
        layout_main_view.addWidget(self.toolbar)
        layout_main_view.addWidget(self.sc)
        layout_main_view.addWidget(self.slider_molecules)
        layout_main_view.addWidget(self.slider_membrane)
        layout_main_view.addLayout(self.slider_from_to_layout)
        layout_main_view.addLayout(self.animation_layout) ## TODO: solve popup windows
        layout_main_view.addLayout(self.set_frame_layout)
        # AREA FOR ADDITIONAL ANALYSIS
        self.additional_data_view = BottomViewGroup(self)
        self.additional_data_view.setHidden(True)
        layout_main_view.addWidget(self.additional_data_view)
        self.layout.addLayout(layout_main_view, 0, 0)

        # MENU BAR
        open_file = MenuAction(self, name='&Open', key=QKeySequence('Ctrl+O'), function='self.openFile')
        save_file = MenuAction(self, name='&Save to XPDB', key=QKeySequence('Ctrl+S'), function='self.saveFile')
        save_selection = MenuAction(self, name='&Save selection to XPDB', function='self.saveFile')
        save_selection.setAccessibleName('save_selection')
        save_legend_main_view = MenuAction(self, name='&Save Figure Legend', function='self.savelegendMainView')
        change_frames = MenuAction(self, name='&Select frames', function='self.changeFrames')
        map_distance = MenuAction(self, name='&Map intramolecular distance', function='self.mapIntraDistance')
        map_speed = MenuAction(self, name='&Map Speed', function='self.mapSpeed')
        map_position = MenuAction(self, name='&Map Position', function='self.mapPosition')
        reset_view = MenuAction(self, name='&Reset', function='self.resetView')
        below_speed = MenuAction(self, name='&Speed', function='self.belowSpeed')
        below_speed.setAccessibleName('single')
        below_speed_multiple = MenuAction(self, name='&Multiple Speed', function='self.belowSpeed')
        below_speed_multiple.setAccessibleName('multiple')
        below_distance = MenuAction(self, name='&Distance', function='self.belowDistance')
        below_distance.setAccessibleName('single')
        below_distance_multiple = MenuAction(self, name='&Multiple Distance', function='self.belowDistance')
        below_distance_multiple.setAccessibleName('multiple')
        show_external_xvg = MenuAction(self, name='&Show XY-XVG file', function='self.showExternalXVG')
        about_action = MenuAction(self, name='&About MeTrEx', function='self.showAbout')
        help_action = MenuAction(self, name='&Show Documentation', function='self.showDocumentation')
        

        menu = self.menuBar()
        m = menu.addMenu('&File')
        m.addAction(open_file)
        m.addSeparator()
        m.addAction(save_file)
        m.addAction(save_selection)
        m.addAction(save_legend_main_view)
        m = menu.addMenu('&View')
        m.addAction(change_frames)
        m.addSeparator()
        m.addAction(map_position)
        m.addAction(map_distance)
        m.addAction(map_speed)
        m.addSeparator()
        m.addAction(reset_view)
        m.addSeparator()
        m_sub = m.addMenu('&Show below')
        m_sub.addAction(below_speed)
        m_sub.addAction(below_speed_multiple)
        m_sub.addAction(below_distance)
        m_sub.addAction(below_distance_multiple)
        m = menu.addMenu('&Analysis')
        m.addAction(show_external_xvg)
        m = menu.addMenu('&Help')
        m.addAction(help_action)
        m.addAction(about_action)

    def get_current_window_size(self):
        """Get current window size of main window"""
        size = self.size()  # Get current window size
        self.width = size.width()  # Get width
        self.height = size.height()  # Get height
    
    def linkSliders(self, s):
        """Link and unlink the line representation and membrane sliders."""
        self.sliders_linked = s

    def showAnimation(self, s):
        """Play animation."""
        start = self.slider_position_molecules
        end = self.data.max_frame
        self.show_animation = s
        pixmap_media_pause = QStyle.StandardPixmap.SP_MediaPause
        animation_icon_pause = self.style().standardIcon(pixmap_media_pause)
        QCoreApplication.processEvents() #
        if s:
            self.sliders_linked = True
            self.connect_sliders.setChecked(True)
            self.animation_button.setIcon(animation_icon_pause)
        while self.show_animation:
            self.sliderChangedMolecules(start)
            QCoreApplication.processEvents()
            start += 1
            if start > end and s:
                start = 0
        pixmap_media_play = QStyle.StandardPixmap.SP_MediaPlay
        animation_icon_play = self.style().standardIcon(pixmap_media_play)
        self.animation_button.setIcon(animation_icon_play)
        self.animation_button.setChecked(False)

    def setOverviewLayout(self):
        self.overview_data_groupbox = QGroupBox('General Information')
        self.overview_data_groupbox.setAccessibleName('overview_data')
        self.overview_data_groupbox.setFixedHeight(340)
        self.overview_data_layout = QVBoxLayout()
        self.overview_data_groupbox.setMaximumWidth(self.size_info_max)

        # DYNAMIC INFO - TIME MOLECUES SLIDER
        time_layout = QGridLayout()
        self.time_label = QLabel('{:,}'.format(self.data.timesteps_mainview[self.slider_position_molecules]).replace(',',' '))
        font = QFont()
        font.setBold(True)
        font.setPixelSize(30)
        self.time_label.setFont(font)
        self.time_label.setStyleSheet('color: fuchsia')
        time_unit = QLabel(str(self.data.unit_time))
        time_layout.addWidget(self.time_label, 0, 0)
        time_layout.addWidget(time_unit, 0, 1)
        time_layout.setAlignment(time_unit, Qt.AlignmentFlag.AlignBottom)
        time_layout.setAlignment(Qt.AlignmentFlag.AlignRight)
        # DYNAMIC INFO - FRAME MOLECULES SLIDER
        self.frame_label = QLabel(str(self.slider_position_molecules+1)) # from 1 to n+1
        font = QFont()
        font.setBold(True)
        font.setPixelSize(20)
        self.frame_label.setFont(font)
        time_layout.addWidget(self.frame_label, 1, 0, Qt.AlignmentFlag.AlignRight)

        # DYNAMIC INFO - TIME MEMBRANE SLIDER
        time_layout_membrane = QGridLayout()
        self.time_label_membrane = QLabel('{:,}'.format(self.data.timesteps_mainview[self.slider_position_membrane]).replace(',',' '))
        font_membrane = QFont()
        font_membrane.setBold(True)
        font_membrane.setPixelSize(30)
        self.time_label_membrane.setFont(font_membrane)
        self.time_label_membrane.setStyleSheet('color: fuchsia')
        time_unit_membrane = QLabel(str(self.data.unit_time))
        time_layout_membrane.addWidget(self.time_label_membrane, 0, 0)
        time_layout_membrane.addWidget(time_unit_membrane, 0, 1)
        time_layout_membrane.setAlignment(time_unit_membrane, Qt.AlignmentFlag.AlignBottom)
        time_layout_membrane.setAlignment(Qt.AlignmentFlag.AlignRight)
        # DYNAMIC INFO _ FRAME MEMBRANE SLIDER
        self.frame_label_membrane = QLabel(str(self.slider_position_membrane+1)) # from 1 to n+1
        font = QFont()
        font.setBold(True)
        font.setPixelSize(20)
        self.frame_label_membrane.setFont(font)
        time_layout_membrane.addWidget(self.frame_label_membrane, 1, 0, Qt.AlignmentFlag.AlignRight)

        # STATIC INFO
        shown_frames_layout = QVBoxLayout()
        self.shown_frames = QLabel('{:,}'.format(self.data.max_frame+1).replace(',',' ')+' frames shown')
        total_frames = QLabel('{:,}'.format(self.data.original_frame_number).replace(',',' ')+' frames total')
        frames_skipped = QLabel(str(self.data.n)+' frame(s) skipped')
        atoms = QLabel('{:,}'.format(self.data.n_atoms).replace(',',' ')+' atoms')
        molecules = QLabel('{:,}'.format(len(self.data.distinct_molecules)).replace(',',' ')+' molecules')
        anker_atom = QLabel(f'{(self.data.anker)} as anker atom')
        if self.data.max_frame == 0:
            total_frames.setToolTip('This value might not be correct.\nMaybe youre selection of \'k\' or \'n\' in the preprocessing step\nwas not appropriate for your data.')
            frames_skipped.setToolTip('This value might not be correct.\nMaybe youre selection of \'k\' or \'n\' in the preprocessing step\nwas not appropriate for your data.')
        shown_frames_layout.addWidget(self.shown_frames)
        shown_frames_layout.addWidget(total_frames)
        shown_frames_layout.addWidget(frames_skipped)
        shown_frames_layout.addWidget(atoms)
        shown_frames_layout.addWidget(molecules)
        shown_frames_layout.addWidget(anker_atom)
        shown_frames_layout.setAlignment(self.shown_frames, Qt.AlignmentFlag.AlignRight)
        shown_frames_layout.setAlignment(total_frames, Qt.AlignmentFlag.AlignRight)
        shown_frames_layout.setAlignment(frames_skipped, Qt.AlignmentFlag.AlignRight)
        shown_frames_layout.setAlignment(atoms, Qt.AlignmentFlag.AlignRight)
        shown_frames_layout.setAlignment(molecules, Qt.AlignmentFlag.AlignRight)
        shown_frames_layout.setAlignment(anker_atom, Qt.AlignmentFlag.AlignRight)

        self.overview_data_layout.addLayout(time_layout)
        self.overview_data_layout.addLayout(time_layout_membrane)
        self.overview_data_layout.addLayout(shown_frames_layout)
        self.overview_data_groupbox.setLayout(self.overview_data_layout)
        self.controls_layout.addWidget(self.overview_data_groupbox)

    def updateCurrentData(self):
        """Update the overview panel showing the current frame number, etc.."""
        if self.show_selection:
            offset = self.selection_from - 1
        else:
            offset = 0
        self.time_label.setText('{:,}'.format(self.data.timesteps_mainview[self.slider_position_molecules]).replace(',',' '))
        self.frame_label.setText(str(self.slider_position_molecules+1+offset)) # from 1 to n+1
        self.time_label_membrane.setText('{:,}'.format(self.data.timesteps_mainview[self.slider_position_membrane]).replace(',',' '))
        self.frame_label_membrane.setText(str(self.slider_position_membrane+1+offset)) # from 1 to n+1

    def updateView(self):
        """Setup/update the main view according to self.data and self.sliderPosition"""
        if self.data.top is not None:
            if len(self.data.distinct_molecule_types) == 0:
                pass
            else:
                self.data.generateView()
                self.slider_molecules.setHidden(False)
                self.slider_molecules.setMaximum(self.data.max_frame) # Minimum value of slider is 0
                self.slider_membrane.setHidden(False)
                self.slider_membrane.setMaximum(self.data.max_frame) # Minimum value of slider is 0
                self.slider_from.setHidden(False)
                self.slider_from.setText(str(self.data.min_frame+1)) # from 1 to n+1
                self.slider_to.setHidden(False)
                self.slider_to.setText(str(self.data.max_frame+1)) # from 1 to n+1
                self.slider_description.setHidden(False)
                self.set_frame.setMinimum(self.data.min_frame+1) # index + 1
                self.set_frame.setMaximum(self.data.max_frame+1) # index + 1
                self.set_frame.setHidden(False)
                self.set_frame_label.setHidden(False)
                self.animation_button.setHidden(False)
                self.connect_sliders.setHidden(False)
                self.setOverviewLayout()
                self.controlButtons()
            self.redrawMainView()
        else:
            pass

    def controlButtons(self):
        """Add buttons on the MainWindow to steer the visibility of the drug molecules in the main visualization."""
        # Settings change buttons
        self.change_settings_groupbox = QGroupBox('Settings')
        self.change_settings_groupbox.setAccessibleName('settings_change_box')
        self.change_settings_groupbox.setFixedHeight(250)
        self.change_settings_groupbox.setMaximumWidth(self.size_info_max)
        self.change_settings_layout = QVBoxLayout()
        # CHANGE COLORMAP
        change_colormap = QPushButton('Colormap')
        change_colormap.setToolTip('Change the colormap for all line represenatation molecules.\nIf greyC selected, also membrane color changes.')
        change_colormap.clicked.connect(self.changeColorMap)
        # CHANGE COLOR
        change_color = QPushButton('Color')
        change_color.clicked.connect(self.changeColor)
        # CHANGE MEMBRANE REPRESENTATION
        change_membrane_representation = QPushButton('Membrane original')
        change_membrane_representation.setToolTip('original: phosphor atoms positions \n abstraction: regression of phosphor atom positions')
        change_membrane_representation.setCheckable(True)
        change_membrane_representation.clicked.connect(self.changeMembraneRepresenation)
        # CHANGE POLYNOM FOR ABSTRACT REPRESENTATION
        polynom_layout = QHBoxLayout()
        membrane_changes_layout = QVBoxLayout()
        change_polynom_layout = QHBoxLayout()
        polynom_label = QLabel('Polynomial:')
        polynom_label
        polynom_label.setToolTip('Change the polynom of \n the regression function \n which is used to compute \n the abtract membrane surface.')
        
        self.polynom_spin = QSpinBox()
        self.polynom_spin.setAccessibleName('polynom_change')
        self.polynom_spin.setMinimum(3)
        self.polynom_spin.setMaximum(15)
        self.polynom_spin.setValue(self.data.polynomial)
        self.polynom_spin.setMinimumSize(50, self.polynom_spin.minimumSizeHint().height())
        self.polynom_current = QLabel('')
        self.polynom_current.setText(str(self.data.polynomial))
        self.polynom_current.setMinimumSize(50, self.polynom_spin.minimumSizeHint().height())

        change_polynom_layout.addWidget(self.polynom_current)
        change_polynom_layout.addWidget(self.polynom_spin)
        membrane_extension_layout = QHBoxLayout()
        membrane_extension_label = QLabel('Expansion:')
        membrane_extension_label.setToolTip('Set the expansion of \n the membrane surface in percentage \n to prevent an undesired decline \n of the membrane at its edges.')

        self.expansion_spin = QSpinBox()
        self.expansion_spin.setAccessibleName('membrane_expansion')
        self.expansion_spin.setRange(5, 30)
        self.expansion_spin.setValue(self.data.membrane_extension)
        self.expansion_spin.setMinimumSize(50, self.expansion_spin.minimumSizeHint().height())
        self.membrane_extension_current = QLabel('' + ' %')
        self.membrane_extension_current.setText(str(self.data.membrane_extension) + ' %')
        self.membrane_extension_current.setMinimumSize(50, self.expansion_spin.minimumSizeHint().height())

        membrane_extension_layout.addWidget(self.membrane_extension_current)
        membrane_extension_layout.addWidget(self.expansion_spin)
        membrane_changes_layout.addWidget(polynom_label)
        membrane_changes_layout.addLayout(change_polynom_layout)
        membrane_changes_layout.addWidget(membrane_extension_label)
        membrane_changes_layout.addLayout(membrane_extension_layout)
        polynom_button = QPushButton('Apply')
        polynom_button.clicked.connect(self.updatePolynom)
        polynom_button.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)
        polynom_button.setMaximumHeight(90)
        polynom_layout.addLayout(membrane_changes_layout)
        polynom_layout.addWidget(polynom_button)
        #polynom_layout.setSizeConstraint()
        # ADD WIDGETS AND LAYOUTS
        self.change_settings_layout.addWidget(change_colormap)
        self.change_settings_layout.addWidget(change_color)
        self.change_settings_layout.addWidget(change_membrane_representation)
        self.change_settings_layout.addLayout(polynom_layout, 0)
        self.change_settings_groupbox.setLayout(self.change_settings_layout)
        self.controls_layout.addWidget(self.change_settings_groupbox)

        # Show / hide drug molecules and leaflets in main view buttons
        self.show_hide_groupbox = QGroupBox('Disable/Enable Visibility')
        self.show_hide_groupbox.setAccessibleName('show_hide_box')
        self.show_hide_groupbox.setMaximumWidth(self.size_info_max)              ####
        self.show_hide_groupbox.setContentsMargins(0, 7, 0, 0)

        self.wrapper_wrapper_layout = QVBoxLayout()
        self.show_hide_info_wrapper = QWidget()
        self.show_hide_info_scroll = QScrollArea()
        self.show_hide_info_scroll.setFrameShape(QFrame.Shape.NoFrame)

        self.show_hide_layout = QGridLayout()
        self.mapping_unit_label = QLabel('Analysis, Unit')
        self.mapping_unit_label.setAccessibleName('analysis_unit_label')
        self.mapping_unit_label.setMinimumWidth(self.size_analysis_unit_label)
        self.mapping_unit_label.setHidden(False)
        self.show_hide_layout.addWidget(self.mapping_unit_label, 0, 0, Qt.AlignmentFlag.AlignLeft)
        self.mapping_min_label = QLabel('min')
        self.mapping_min_label.setAlignment(Qt.AlignmentFlag.AlignRight)
        self.mapping_min_label.setMinimumWidth(35)#40
        self.mapping_min_label.setAccessibleName('analysis_min_label')
        self.mapping_min_label.setHidden(False)
        self.mapping_max_label = QLabel('max')
        self.mapping_max_label.setAlignment(Qt.AlignmentFlag.AlignRight)
        self.mapping_max_label.setMinimumWidth(35)#40
        self.mapping_max_label.setAccessibleName('analysis_max_label')
        self.mapping_max_label.setHidden(False)
        self.show_hide_layout.addWidget(self.mapping_min_label, 0, 1, Qt.AlignmentFlag.AlignRight)
        self.show_hide_layout.addWidget(self.mapping_max_label, 0, 2, Qt.AlignmentFlag.AlignRight)
        i = 1
        for key, value in sorted(self.data.drug_molecules_colors.items()):
            button = QPushButton(key)
            button.setCheckable(True)
            button.setAccessibleName('show_hide_pushbutton') # needed to change connection
            _color_r = int(self.data.drug_molecules_colors[key][0]*255)
            _color_g = int(self.data.drug_molecules_colors[key][1]*255)
            _color_b = int(self.data.drug_molecules_colors[key][2]*255)
            font_brightness = int(sqrt(0.299*_color_r**2 + 0.587*_color_g**2 + 0.144*_color_b**2).real)
            if font_brightness < 128:
                button.setStyleSheet('background-color: rgb('+str(_color_r)+','+str(_color_g)+','+str(_color_b)+'); color: white;')
            else:
                button.setStyleSheet('background-color: rgb('+str(_color_r)+','+str(_color_g)+','+str(_color_b)+'); color: black;')
            button.clicked.connect(self.showHideMolecule)
            self.show_hide_layout.addWidget(button, i, 0, Qt.AlignmentFlag.AlignLeft)
            layout_min_label = QVBoxLayout()
            min_label = QLabel('')
            min_label.setAccessibleName(key+'min')
            min_label.setHidden(False)
            min_label_frame = QLabel('')
            min_label_frame.setToolTip('Frame number of minima.')
            min_label_frame.setAccessibleName(key+'minframe')
            min_label_frame.setHidden(False)
            layout_min_label.addWidget(min_label, 0, Qt.AlignmentFlag.AlignRight)
            layout_min_label.addWidget(min_label_frame, 0, Qt.AlignmentFlag.AlignRight)
            self.show_hide_layout.addLayout(layout_min_label, i, 1, Qt.AlignmentFlag.AlignRight)
            layout_max_label = QVBoxLayout()
            max_label = QLabel('')
            max_label.setAccessibleName(key+'max')
            max_label.setHidden(False)
            max_label_frame = QLabel('')
            max_label_frame.setToolTip('Frame number of maxima.')
            max_label_frame.setAccessibleName(key+'maxframe')
            max_label_frame.setHidden(False)
            layout_max_label.addWidget(max_label, 0, Qt.AlignmentFlag.AlignRight)
            layout_max_label.addWidget(max_label_frame, 0, Qt.AlignmentFlag.AlignRight)
            self.show_hide_layout.addLayout(layout_max_label, i, 2, Qt.AlignmentFlag.AlignRight)
            i += 1
        self.wrapper_wrapper_layout.addWidget(self.show_hide_info_scroll)

        self.show_hide_info_wrapper.setLayout(self.show_hide_layout)
        self.show_hide_info_scroll.setWidget(self.show_hide_info_wrapper)
        if len(self.data.lipid_molecules_reference) > 0:
            leaflet0 = QPushButton('upper leaflet')
            leaflet0.setCheckable(True)
            leaflet0.clicked.connect(self.showHideMembraneLeaflet)
            leaflet1 = QPushButton('lower leaflet')
            leaflet1.setCheckable(True)
            leaflet1.clicked.connect(self.showHideMembraneLeaflet)
            self.wrapper_wrapper_layout.addWidget(leaflet0, 0, Qt.AlignmentFlag.AlignLeft)
            self.wrapper_wrapper_layout.addWidget(leaflet1, 0, Qt.AlignmentFlag.AlignLeft)
        
        self.show_hide_groupbox.setLayout(self.wrapper_wrapper_layout)
        self.controls_layout.addWidget(self.show_hide_groupbox)

    def changeColorMap(self):
        """"""
        colors = colormap.color_maps + colormap.color_map_grey
        color = ChangeColormapDialog(self.data.colormap, colors)
        color.exec()
        if len(color.colormap) == 0:
            return
        else:
            # set new color map
            self.data.changeColormap(color.colormap)
            # changes color for line representation molecules in main view
            self.redrawColor()
            # change leaflets
            self.redraw()
            # change colors of the bottom views and control buttons
            buttons = self.findChildren(QPushButton)
            for molecule, color in self.data.drug_molecules_colors.items():
                self.changeColorBottomViews(molecule, color)
                for button in buttons:
                    if button.text() == molecule:
                        _color_r = int(color[0]*255)
                        _color_g = int(color[1]*255)
                        _color_b = int(color[2]*255)
                        font_brightness = int(sqrt(0.299*_color_r**2 + 0.587*_color_g**2 + 0.144*_color_b**2).real)
                        if font_brightness < 128:
                            button.setStyleSheet('background-color: rgb('+str(_color_r)+','+str(_color_g)+','+str(_color_b)+'); color: white;')
                        else:
                            button.setStyleSheet('background-color: rgb('+str(_color_r)+','+str(_color_g)+','+str(_color_b)+'); color: black;')

    def updatePolynom(self):
        """Update the polynom for the abstract membrane representation."""
        reply = QMessageBox.question(
            self,
            'Change Membrane Abstraction',
            'Are you sure you want to recalculate the membrane abstraction?\n\n'+\
                'Depending on the chosen polynom, the percentage of membrane\n'+\
                'expansion, the numer of lipid molecules and the number of frames\n'+\
                'this might take some time to compute.',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes
            )
        if reply == QMessageBox.StandardButton.Yes:
            for child in self.findChildren(QSpinBox):
                if child.accessibleName() == 'polynom_change':
                    polynom = child.value()
                if child.accessibleName() == 'membrane_expansion':
                    extension = child.value()
            if self.data.polynomial == polynom and self.data.membrane_extension == extension:
                QMessageBox.information(
                    self,
                    'Calculation Finished',
                    'You did not specify new values.\nNo calculation needed.',
                    QMessageBox.StandardButton.Ok,
                    QMessageBox.StandardButton.Ok
                    )
                return
            else:
                self.data.polynomial = polynom
                self.data.membrane_extension = extension
                self.data.calcLipidSurfaces()
                self.data.lipid_molecules_surface['leaflet0'] = np.array(self.data.lipid_molecules_surface['leaflet0'])
                self.data.lipid_molecules_surface['leaflet1'] = np.array(self.data.lipid_molecules_surface['leaflet1'])
                self.redraw()
                self.polynom_current.setText(str(polynom))
                self.membrane_extension_current.setText(str(extension) + ' %')
                QMessageBox.information(
                    self,
                    'Calculation Finished',
                    'The membrane surface has been recalculated.',
                    QMessageBox.StandardButton.Ok,
                    QMessageBox.StandardButton.Ok
                    )
        else:
            return

    def changeMembraneRepresenation(self, checked):
        """Change the membrane represenation and update plot."""
        if checked:
            self.sender().setText('Membrane abstraction')
            self.data.lipid_molecules_show_abstraction = False
            self.redraw()
        else:
            self.sender().setText('Membrane original')
            self.data.lipid_molecules_show_abstraction = True
            self.redraw()

    def setFrame(self):
        """Update view according to selected frame number."""
        if self.show_selection:
            offset = self.selection_from - 1
        else:
            offset = 0
        for child in self.findChildren(QSpinBox):
            if child.accessibleName() == 'set_frame':
                position = child.value()-1 - offset
                self.slider_position_molecules = position
                self.slider_position_membrane = position
                self.slider_molecules.setValue(position)
                self.slider_membrane.setValue(position)
                if self.sliders_linked:
                    self.updateBottomViewSliders()
                self.updateCurrentData()
                self.redraw()

    def sliderChangedMolecules(self, s):
        """Set slider variable for the drug molecules to current position and call redrawMainView()."""
        self.slider_position_molecules = s
        if self.sliders_linked:
            self.slider_position_membrane = s
            self.slider_membrane.setValue(s)
            if self.show_selection:
                self.sliderChanged.emit(s+self.selection_from)
            else:
                self.sliderChanged.emit(s+1) # TODO: for slider linkage with SubWindow
        self.updateCurrentData()
        self.redraw()
        if self.show_selection:
            self.set_frame.setValue(s+self.selection_from)
        else:
            self.set_frame.setValue(s+1)

    def updateBottomViewSliders(self):
        """Sets the current slider position value to all QSpinBox
        objects and emits a editingFinished signal.
        """
        if self.show_selection:
            offset = self.selection_from - 1
        else:
            offset = 0
        if self.sliders_linked:
            if len(self.additional_data_view.bottom_views) > 0:
                for child in self.additional_data_view.findChildren(QSpinBox):
                    if child.accessibleName() != 'set_frame': #== 'jump_to': #does not work this way?!?
                        child.setValue(self.slider_position_molecules + 1 + offset) # for correct index
                        child.editingFinished.emit()

    def sliderChangedMembrane(self, s):
        """Set slider variable for the membrane to current position and call redraw()."""
        self.slider_position_membrane = s
        if self.sliders_linked:
            self.slider_position_molecules = s
            self.slider_molecules.setValue(s)
        self.updateCurrentData()
        self.redraw()
        if self.show_selection:
            self.set_frame.setValue(s+self.selection_from)
        else:
            self.set_frame.setValue(s+1)

    def showHideMolecule(self):
        """Set show variable of line representation molecules and call redraw()."""
        pushButtons = self.findChildren(QPushButton)
        for button in pushButtons:
            if button.text() in self.data.drug_molecules_show:
                if button.isChecked():
                    self.data.drug_molecules_show[button.text()] = False
                if not button.isChecked():
                    self.data.drug_molecules_show[button.text()] = True
        if self.data.drug_molecules_show_mapping:
            self.redrawShowHideMapping()
        else:
            self.redrawShowHide()
    
    def showHideMembraneLeaflet(self):
        """Set show variable of line representation molecules and call redraw()."""
        pushButtons = self.findChildren(QPushButton)
        for button in pushButtons:
            if button.text() == 'upper leaflet':
                if button.isChecked():
                    self.data.lipid_molecules_show[0] = False
                if not button.isChecked():
                    self.data.lipid_molecules_show[0] = True
            if button.text() == 'lower leaflet':
                if button.isChecked():
                    self.data.lipid_molecules_show[1] = False
                if not button.isChecked():
                    self.data.lipid_molecules_show[1] = True
        self.redrawShowHide()

    def clearControls(self):
        """Clear controls from control grid."""
        widgets_to_clear = self.findChildren(QGroupBox)
        for i in widgets_to_clear:
            if i.accessibleName() == 'settings_change_box' or i.accessibleName() == 'show_hide_box' or i.accessibleName() == 'overview_data':
                i.setParent(None)
    
    def redrawMainView(self):
        """This function clears and redraws the main view without any
        mapping. The main view will look the same like when drawn for
        the first time with the given data.
        """
        self.toolbar.update()
        self.sc.ax.cla()
        self.sc.updateView(
            False,
            positionMolecules=self.slider_position_molecules,
            positionMembrane=self.slider_position_membrane,
            data=self.data
            )
        self.sc.draw()
        self.toolbar.home()
    
    # access self.sc directly
    def changeColorBottomViews(self, name, color):
        """Change the color of the identification label and plotted
        line in the bottom views by calling the changeColor funciton
        of the BottomViews.

        Parameter:
        -------
        name ():
            Identification of the line to change. Usually this is the
            name and residue number of one of the line-representation-
            molecules.
        color ():
            The new color.
        """
        try:
            for bv, reference in self.additional_data_view.bottom_views.items():
                d = {}
                d[name] = color
                reference.changeColor(d)
        except:
            # There were no BottomViews.
            return
    
    # access self.sc directly
    def redraw(self):
        """Remove and redraw spheres and leaflets of the main view."""
        for sphere, reference in self.sc.spheres.items():
            positions = self.data.drug_molecules_positions[sphere]
            xpos = positions[:,0]
            ypos = positions[:,1]
            zpos = positions[:,2]
            reference.remove()
            self.sc.spheres[sphere] = self.sc.ax.scatter(
                xpos[self.slider_position_molecules],
                ypos[self.slider_position_molecules],
                zpos[self.slider_position_molecules],
                color='fuchsia',
                s=50
                )
        for leaflet, reference in self.sc.leaflets.items():
            try:
                reference.remove()
                col = self.data.leaflet_colors[leaflet]
                if self.data.lipid_molecules_show_abstraction:
                    self.sc.leaflets[leaflet] = self.sc.ax.plot_surface(
                        self.data.lipid_molecules_surface[leaflet][self.slider_position_membrane][0],
                        self.data.lipid_molecules_surface[leaflet][self.slider_position_membrane][1],
                        self.data.lipid_molecules_surface[leaflet][self.slider_position_membrane][2],
                        color=col,
                        antialiased=True,
                        alpha=0.2
                        )
                else:
                    self.sc.leaflets[leaflet] = self.sc.ax.plot_trisurf(
                        self.data.lipid_molecules_positions[leaflet][self.slider_position_membrane][:,0],
                        self.data.lipid_molecules_positions[leaflet][self.slider_position_membrane][:,1],
                        self.data.lipid_molecules_positions[leaflet][self.slider_position_membrane][:,2],
                        color=col,
                        linewidth=0.2,
                        antialiased=True,
                        alpha=0.2,
                        edgecolor='grey'
                        )
            except:
                pass
        self.sc.draw()
    
    # access self.sc directly
    def redrawShowHideMapping(self):
        """Remove and redraw elements of the main view while a mapping
        is applied.
        """
        # setup NaN values for the line to skip, if any
        hide = np.empty((1, 2, 3))
        hide[:] = np.NaN

        # keep original, but only if not existing jet
        if self.data.drug_molecules_mseg_original is None:
            self.data.drug_molecules_mseg_original = deepcopy(self.data.drug_molecules_mseg)
        else:
            self.data.drug_molecules_mseg = deepcopy(self.data.drug_molecules_mseg_original)

        # remove all elements not needed from the view
        # (spheres, labels, markers from extremes, colorbar,
        # line collection)
        for sphere, reference in self.sc.spheres.items():
            reference.remove()
        self.sc.spheres = {}
        for label, reference in self.sc.line_labels.items():
            reference.remove()
        self.sc.line_labels = {}
        for ex, refs in self.sc.extreme_labels.items():
            try:
                # if you set min and may labels, change here
                refs[0][0].remove()
                #refs[0][1].remove()
                refs[1][0].remove()
                #refs[1][1].remove()
            except:
                # no labels in this view
                pass
        self.sc.extreme_labels = {}
        # colorbar needs to be removed because it refers to data that
        # will change by this function
        self.sc.colorbar.remove()
        self.sc.linecoll.remove()

        # set needed elements of view
        for mol, value in self.data.drug_molecules_show.items():
            if value == False:
                # exchange linecoll to skip not shown lines
                position = self.data.drug_molecules_morder.index(mol)
                start_pos = position * (self.data.original_frame_number - 1)
                end_pos = start_pos + self.data.max_frame
                while start_pos < end_pos :
                    self.data.drug_molecules_mseg[start_pos] = hide
                    start_pos +=1
            else:
                # draw sphere, label, extremes
                positions = self.data.drug_molecules_positions[mol]
                xpos = positions[:,0]
                ypos = positions[:,1]
                zpos = positions[:,2]
                self.sc.line_labels[mol] = self.sc.ax.text(xpos[0], ypos[0], zpos[0], mol, None)
                self.sc.spheres[mol] = self.sc.ax.scatter(
                    xpos[self.slider_position_molecules],
                    ypos[self.slider_position_molecules],
                    zpos[self.slider_position_molecules],
                    color='fuchsia',
                    s=50
                    )
                self.sc.extreme_labels[mol] = []
                try:
                    n = len(self.data.drug_molecules_text[mol][0])
                    i = 0
                    while i < n:
                        ind = self.data.drug_molecules_text[mol][0][i]
                        # Change here, if you want to have min and max text labels
                        #txt = 'min: '+str(int(self.data.drug_molecules_text[mol][1]))
                        #label = self.sc.ax.text(xpos[ind], ypos[ind], zpos[ind], txt, None)
                        #label.set_fontsize('small')
                        marker = self.sc.ax.scatter(xpos[ind], ypos[ind], zpos[ind], marker='d',color='fuchsia', s=20)
                        #self.sc.extreme_labels[mol].append([label, marker])
                        self.sc.extreme_labels[mol].append([marker])
                        i += 1
                    j = 0
                    n = len(self.data.drug_molecules_text[mol][2])
                    while j < n:
                        ind = self.data.drug_molecules_text[mol][2][j]
                        # Change here, if you want to have min and max text labels
                        #txt = 'max: '+str(int(self.data.drug_molecules_text[mol][3]))
                        #label = self.sc.ax.text(xpos[ind], ypos[ind], zpos[ind], txt, None)
                        #label.set_fontsize('small')
                        marker = self.sc.ax.scatter(xpos[ind], ypos[ind], zpos[ind], marker='P',color='fuchsia', s=20)
                        #self.sc.extreme_labels[mol].append([label, marker])
                        self.sc.extreme_labels[mol].append([marker])
                        j += 1
                except:
                    # no extreme values available (e.g. map position)
                    pass
        self.data.drug_molecules_linecol = Line3DCollection(self.data.drug_molecules_mseg, cmap=plt.get_cmap('viridis')) # TODO: chekc if color is shown correctly
        # lines will be colored according to set_array
        self.data.drug_molecules_linecol.set_array(self.data.drug_molecules_mvalues)
        self.sc.linecoll = self.sc.ax.add_collection3d(self.data.drug_molecules_linecol)
        self.sc.colorbar = self.sc.fig.colorbar(self.sc.linecoll, ax=self.sc.ax, shrink=0.7, pad=0.13, label=self.data.drug_molecules_mlabel)

        # update plot
        self.sc.draw()
    
    # access self.sc directly
    def redrawShowHide(self):
        """Update view."""
        # show/hide button clicked
        for sphere, reference in self.sc.spheres.items():
            reference.remove()
        self.sc.spheres = {}
        for molecule, value in self.data.drug_molecules.items():
            try:
                reference = self.sc.lines[molecule][0].pop()
                reference.remove()
                reference_text = self.sc.lines[molecule][1]
                reference_text.remove()
                del self.sc.lines[molecule]
            except:
                pass
            if self.data.drug_molecules_show[molecule]:
                try:
                    positions = self.data.drug_molecules_positions[molecule]
                    xpos = positions[:,0]
                    ypos = positions[:,1]
                    zpos = positions[:,2]
                    line = self.sc.ax.plot(xpos, ypos, zpos, color=self.data.drug_molecules_colors[molecule])
                    text = self.sc.ax.text(xpos[0], ypos[0], zpos[0], molecule, None)
                    self.sc.lines[molecule] = [line, text]
                    self.sc.spheres[molecule] = self.sc.ax.scatter(
                        xpos[self.slider_position_molecules],
                        ypos[self.slider_position_molecules],
                        zpos[self.slider_position_molecules],
                        color='fuchsia',
                        s=50
                        )
                except: # there are molecules without position information
                    pass
            else:
                pass
        for leaflet, reference in self.sc.leaflets.items():
            try:
                reference.remove()
            except:
                pass
            if leaflet == 'leaflet0':
                ind = 0
            else:
                ind = 1
            if self.data.lipid_molecules_show[ind]:
                col = self.data.leaflet_colors[leaflet]
                if self.data.lipid_molecules_show_abstraction:
                    self.sc.leaflets[leaflet] = self.sc.ax.plot_surface(
                        self.data.lipid_molecules_surface[leaflet][self.slider_position_membrane][0],
                        self.data.lipid_molecules_surface[leaflet][self.slider_position_membrane][1],
                        self.data.lipid_molecules_surface[leaflet][self.slider_position_membrane][2],
                        color=col,
                        antialiased=True,
                        alpha=0.2
                        )
                else:
                    self.sc.leaflets[leaflet] = self.sc.ax.plot_trisurf(
                        self.data.lipid_molecules_positions[leaflet][self.slider_position_membrane][:,0],
                        self.data.lipid_molecules_positions[leaflet][self.slider_position_membrane][:,1],
                        self.data.lipid_molecules_positions[leaflet][self.slider_position_membrane][:,2],
                        color=col,
                        linewidth=0.2,
                        antialiased=True,
                        alpha=0.2,
                        edgecolor='grey'
                        )
        self.sc.draw()
    
    # access self.sc directly
    def redrawColor(self):
        """Apply the color change for a line represenation molecule."""
        for line, value in self.sc.lines.items():
            positions = self.data.drug_molecules_positions[line]
            xpos = positions[:,0]
            ypos = positions[:,1]
            zpos = positions[:,2]
            reference = value[0].pop()
            reference.remove()
            value[0] = self.sc.ax.plot(xpos, ypos, zpos, color=self.data.drug_molecules_colors[line])
        self.sc.draw()

    def changeColor(self):
        """Change color of all objects representing a molecule:
        line in main view and bottom views and the push button.
        """
        dlg = ChangeColorDialog(self.data.drug_molecules, self.data.color_list, 'Select Molecule')
        dlg.exec()
        try:
            color_string = list(dlg.color[1:-1].split(' '))
            color_string = list(filter(lambda a: a != '', color_string))
            color = np.array(list(map(float, color_string)))
            self.data.drug_molecules_colors[dlg.molecule] = color
            # update control panel button
            buttons = self.findChildren(QPushButton)
            for button in buttons:
                if button.text() == dlg.molecule:
                    _color_r = int(color[0]*255)
                    _color_g = int(color[1]*255)
                    _color_b = int(color[2]*255)
                    font_brightness = int(sqrt(0.299*_color_r**2 + 0.587*_color_g**2 + 0.144*_color_b**2).real)
                    if font_brightness < 128:
                        button.setStyleSheet('background-color: rgb('+str(_color_r)+','+str(_color_g)+','+str(_color_b)+'); color: white;')
                    else:
                        button.setStyleSheet('background-color: rgb('+str(_color_r)+','+str(_color_g)+','+str(_color_b)+'); color: black;')
            # update main view and bottom views
            self.redrawColor()
            self.changeColorBottomViews(dlg.molecule, color)
        except: # no color selected / canceled
            return

    def addSubWindow(self, **kwargs):
        """Initialize a SubWindow, register it in the list of
        subwindows and show it.

        Parameters:
        ------
        optional:
        title(string): Name of the SubWindow.
        """
        title = 'test'
        for key, value in kwargs.items():
            if key == 'title':
                title = value
            if key == 'data':
                data = value
        sub_window = SubWindow(self, title=title, data=data)
        self.sub_windows.append(sub_window)
        sub_window.show()

    def closeEvent(self, event):
        """Ensure closing of the main window is intended."""
        reply = QMessageBox.question(
            self,
            'Close',
            'Are you sure you want to close the window?',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes
            )
        if reply == QMessageBox.StandardButton.Yes:
            for i in range(len(self.sub_windows)):
                try:
                    self.sub_windows[i].close()
                except:
                    pass # maybe window was closed already
            event.accept()
        else:
            event.ignore()


class MenuAction(QAction):
    """Menu class shows menu actions.

    Methods:
    ------
    setAccessibleName():
        Set the name of the MenuAction.
    accessibleName():
        Get the name of the MenuAction.
    showExternalXVG():
        Show selected data in an additional window.
    saveFile():
        Save a PDB file from the main view.
    openFile():
        Show selected data in the main view.
    mapPosition():
        Show the time/frame index on the trajectory line.
    mapIntraDistance():
        Show a distance calcualtion on the trajectory lines.
    belowDistance():
        Show a distance calculation in a BottomView.
    showDocumentation():
        Show the documentation of MeTrEx.
    showAbout():
        Show Information about MeTrEx.
    mapSpeed():
        Show a speed calculation on the trajectory lines.
    updateInfoPanelMapping():
        Update the information panel of the main window when a
        mapping is shown.
    belowSpeed():
        Show a speed calculation in a BottomView.
    changeFrames():
        Changes the shown frame selection for the main view.
    resetView():
        Resets the main view to the original view shown.
    saveLegendMainView():
        Saves the legen of the main view.
    showOKCancelMessageBox:
        Confirmation message box is shown.
    showCautionMessageBox():
        Warning is shown.
    """

    def __init__(self, parent, **kwargs):
        """Show actions and connect them to their specific function."""
        super().__init__('&ToBeSpecified', parent)
        self.parent = parent
        self.name = ''

        for arg in kwargs:
            if arg == 'name':
                self.setText(kwargs[arg])
            elif arg == 'key':
                self.setShortcut(kwargs[arg])
            elif arg == 'function':
                self.triggered.connect(eval(kwargs[arg]))

    def setAccessibleName(self, name):
        self.name = name

    def accessibleName(self):
        return self.name

    def showExternalXVG(self):
        """Collect user information, calculate data to be shown,
        setup secondary window and show it.
        """
        dlg = OpenXVGDialog(self.parent.data.drug_molecules)
        dlg.exec()
        if dlg.file_path is None:
            return # rejected, nothing to do
        else:
            data = analyses.readXVGfile(dlg.file_path)
            dim3 = dlg.dim3
            title = data[0]
            x_label = data[1]
            y_label = data[2]
            columns = data[3]
            ind = data[4]
            secondary_y_labels = data[5]
            # setup data object
            d = DataSub()
            d.analysis = title
            d.unit_mpl_x = x_label
            d.unit_mpl_y = y_label
            # create as much colors as columns - 1 to be selected
            colors = returnColors(len(columns), 'plasma')
            if len(dlg.drug_molecules) > 0:
                colors = []
                for e in dlg.drug_molecules:
                    colors.append(self.parent.data.drug_molecules_colors[e])
                if len(colors) < len(columns):
                    n =  len(columns)-len(colors)+1
                    colors.extend(returnColors(n, 'plasma'))
                    colors = np.array(colors)
            # setup singe data sets, add to data object
            i = 1 # columns start at index 1
            j = 0 # colors start at index 0
            while i < len(columns):
                # Create data sets and add them to d
                if len(secondary_y_labels) > 0:
                    try:
                        name = secondary_y_labels[i-1]
                    except:
                        name = 'column '+str(i)
                else:
                    name = 'column '+str(i)
                if dim3:
                    ds = SingleDataSet(columns[0], columns[i], name, colors[j], z=ind)
                    ds.addLabel('z', 'Frame')
                    d.unit_mpl_z = 'Frame'
                else:
                    #c = colors[j]
                    ds = SingleDataSet(columns[0], columns[i], name, colors[j])
                d.addPointset(ds)
                i += 1
                j += 1
            d.calculateExtremes()
            d.slider_range = (1, len(columns[0]))
            # setup sub window
            sub_window = SubWindow(self.parent, title=dlg.file_title.text(), data=d)
            sub_window.slider.setMaximum(len(columns[0]))
            sub_window.slider.setMinimum(1)

            # visualize data
            if dim3:
                sub_window.sc.scatterPlot3D(d)
            else:
                sub_window.sc.scatterPlot2D(d)
            # show window
            sub_window.show()
            self.parent.sub_windows.append(sub_window)

    def saveFile(self):
        """Show needed dialogs for the selection of parameters to save a pdb file."""
        if self.parent.data.top is None:
            self.showCautionMessageBox()
        else:
            selection_strings = set()
            selection_strings.add('all')
            from_frame = self.parent.slider_position_molecules
            to_frame = from_frame + 1
            if self.accessibleName() == 'save_selection': #FIX by Beat, removed .sender()
                dlg = AnalysisSelectionDialog(
                    'all',
                    atoms=self.parent.data.atoms,
                    atoms_dict=self.parent.data.atoms_dict,
                    distinct_mol=self.parent.data.distinct_molecules,
                    distinct_mol_t=self.parent.data.distinct_molecule_types,
                    drug_mol=self.parent.data.drug_molecules,
                    universe=self.parent.data.universe,
                    windowtitle='Save Selection',
                    rad=True
                    )
                dlg.exec()
                selection_strings = dlg.selected_molecules
                if len(selection_strings) == 0:
                    return
                dlg = ChangeViewDataDialog(
                    self.parent,
                    self.parent.data.min_frame + 1,
                    self.parent.data.max_frame + 1,
                    current=str(self.parent.slider_position_molecules + 1),
                    save = True
                    )
                dlg.exec()
                if dlg.from_value is None:
                    return
                from_frame = dlg.from_value
                to_frame = dlg.to_value
            file_filter = 'XPDB file (*.pdb)'
            response = QFileDialog.getSaveFileName(
                parent=self.parent,
                caption='structure.pdb',
                directory=os.getcwd(),
                filter=file_filter,
                initialFilter='XPDB file (*.pdb)'
                )
            try:
                file_path = response[0]
                universe = self.parent.data.universe
                atoms_selection = universe.select_atoms('')
                # apply the following to a newly created pdb file from the trajectory at from_frame-1,
                # in order to select the at the beginning relevant molecules
                unique_selection_strings = set()
                for string in selection_strings:
                    if string[0] == '(':
                        for res in self.parent.data.distinct_molecule_types:
                            res_selection = universe.select_atoms('resname ' + res + ' and ' + string)
                            res_uniques = res_selection.residues.resids
                            for n in res_uniques:
                                unique_selection_strings.add('resname ' + str(res) + ' and resid ' + str(n))
                    else:
                        unique_selection_strings.add(string)
                # select atoms in the Data object's universe and write into file
                for us in unique_selection_strings:
                    atoms_selection += universe.select_atoms(us)
                with mda.Writer(file_path, multiframe=True) as W:
                    # python slicing: [m:n] = start with m, end with not including n
                    for ts in universe.trajectory[from_frame-1:to_frame]:
                        W.write(atoms_selection)
            except:
                return

    def openFile(self):
        """Show open dialogs needed to select file paths and molecules. Update main view."""
        dlg = OpenDialog()
        dlg.exec()
        if dlg.topPath == None or dlg.trajPath == None:
            return
        else:
            # reset slider to position 1
            self.parent.set_frame.setValue(1)
            # initialize data object with given files
            self.parent.data = Data()
            self.parent.data.top = dlg.topPath
            self.data_file_original = deepcopy(dlg.trajPath)
            self.parent.data.data = dlg.trajPath
            self.parent.data.distinctMol() # to be shown in the next dialog

            # collect size of data file, check if fits into memory, inform user about preprocessing
            data_file_size = os.stat(self.parent.data.data).st_size # [Bytes]
            cpu_info = psutil.cpu_freq().current #CPU freq in MHz
            memory_available = psutil.virtual_memory().available # [Bytes]
            if memory_available < data_file_size:
                text = 'Your data file is larger than your availbale memory capacitiy.\n\n'\
                    + 'You should reduce your data generously.\n'\
                    + 'Otherwise this program or even your system might collapse.\n\n'\
                    + 'Guidline: try a file size containing approximately 1000 frames.'
                self.showCautionMessageBox(text)
            elif (memory_available / 2) < data_file_size:
                text = 'Your data file is very large compared to your available memory capacity.\n'\
                    + 'Processing your file this program or even your system might collapse.\n'\
                    + 'Maybe you want to select or skip some frames?\n\n'\
                    + 'Guidline: try a file size containing approimately 1000 frames.'
                self.showCautionMessageBox(text, QMessageBox.Icon.Critical)
            else:
                # file has acceptable size
                pass

            # ask user to specify preprocessing and visualization settings
            do = True
            while do:
                dlg2 = PreprocessingSelectionDialog(molecules=self.parent.data.distinct_molecule_types)
                dlg2.exec()
                if dlg2.selected_molecules is None:
                    # user declined loading data by closing selection window.
                    self.parent.data = Data()
                    return
                # set data collected from dialog
                self.parent.data.line_representation_molecule_types = dlg2.selected_molecules
                self.parent.data.anker_selection = dlg2.anker_selection
                self.parent.data.n = dlg2.n
                self.parent.data.k = dlg2.k
                try:
                    # calculate preprocessing
                    if self.parent.data.k > 1 or self.parent.data.n > 0:
                        print('k:' ,self.parent.data.k, 'n:', self.parent.data.n)
                        self.parent.data.reduceFrames()
                        do = False
                    else:
                        # k = 1 and n = 0, no preprocessing selected
                        do = False
                except MemoryError:
                    print('Data file is still too large.')
                    # inform user about MemoryError - this might not work, as memory is full, though
                    text = 'Your file still seems to be too large to be processed.\nPlease load files again and select\na higher reduction, e.g. increase \'k\'.'
                    self.showCautionMessageBox(text)
                    return
                except OSError as e:
                    print('os error: ', e)
                    text = 'It seems that you selected a number of frames to skip\nwhich is equal or higher the number of frames in your file.\nPlease try again.'
                    self.showCautionMessageBox(text)
                    self.parent.data.data = deepcopy(self.data_file_original)

            if self.parent.data.k > 1 or self.parent.data.n > 0:
                text = 'You reduced the number of frames.\nBe aware that the frames shown are always numbered with\nconsecutive numbers starting from \'1\'!'
                self.showCautionMessageBox(text, QMessageBox.Icon.Information)
            
            # Get CPU frequency for speed estimation
            cpu_freq = psutil.cpu_freq().current + 0.0001
            # caluclate aproximal read time in seconds
            read_time = round((data_file_size * 700) / (cpu_freq * 10**6))
            print(f'Estimated read time: {read_time} seconds')
 
            if read_time < 60:
                text = 'Calculation time largely depends on the number of\nframes and lipid molecules.\n\nYour data will be shown in less than 1 min, approximately.'
            else:
                read_time = round(read_time/60)
                if read_time < 60:
                    if read_time < 5:
                        text = f'Calculation time largely depends on the number of\nframes and lipid molecules.\n\nYour data will be shown in about {read_time} min.'
                    else:
                        read_time_min = round(read_time - read_time * 0.1)
                        read_time_max = round(read_time + read_time * 0.1)
                        text = f'Calculation time largely depends on the number of\nframes and lipid molecules.\n\nYour data will be shown in about {read_time_min} - {read_time_max} min.'
                else:
                    read_time = round(read_time/60)
                    text = f'Calculation time largely depends on the number of\nframes and lipid molecules.\n\nYour data will be shown in about {read_time} hour/s.'
            
            # show message about calculaiton time depending on file size
            # ranges were set due to test data from one simulation
#            if data_file_size < 1e7:
#                text = 'Calculation time largely depends on the number of\nframes and lipid molecules.\n\nYour data will be shown in less than 3 min, approximately.'
#            elif data_file_size < 5e7:
#                text = 'Calculation time largely depends on the number of\nframes and lipid molecules.\n\nYour data will be shown in about 3 - 7 min.'
#            elif data_file_size < 1e8:
#                text = 'Calculation time largely depends on the number of\nframes and lipid molecules.\n\nYour data will be shown in about 4 - 10 min.'
#            else:
#                text = 'Calculation time largely depends on the number of\nframes and lipid molecules.\n\nUnfortenately, we cannot tell how long\nthe calculations will take for your data.'
            self.showCautionMessageBox(text, QMessageBox.Icon.Information)

            # reset the MainView if mapping is shown.
            # (Needed to remove colorbar)
            if self.parent.data.drug_molecules_show_mapping:
                self.resetView(message=False)
            # update main window
            self.parent.clearControls()
            self.parent.updateView()
#FIX by Beat, disables reloead of data after inital loading to avoid crash on MacARM architecture
            if platform.system() == 'Darwin' and platform.machine() == 'x86_64':
                self.setEnabled(False)

    def mapPosition(self):
        """Map the position on the line representation molecules on
        the MainView.
        """
        if self.parent.data.top is None:
            self.showCautionMessageBox()
        else:
            self.resetView(message=False)
            indices = np.arange(self.parent.data.min_frame, self.parent.data.max_frame, 1)

            # Create visualization
            segments = {}
            molecule_order = []
            indices_all = np.array([])
            seg_all = np.array([]).reshape(0, 2, 3)
            i = 0
            try:
                for molecule, positions in self.parent.data.drug_molecules_positions.items():
                    x = positions[:, 0]
                    y = positions[:, 1]
                    z = positions[:, 2]
                    points = np.array([x, y, z]).transpose().reshape(-1, 1, 3)
                    mol_segments = np.concatenate([points[:-1],points[1:]], axis=1)
                    segments[molecule] = mol_segments
                    self.parent.data.drug_molecules_moffset[molecule] = i
                    i += self.parent.data.max_frame
                    molecule_order.append(molecule)
                    indices_all = np.concatenate((indices_all, indices), axis=0)
                    seg_all = np.concatenate((seg_all, mol_segments), axis=0)
            except IndexError:
                self.showCautionMessageBox('It seems your data does not contain information about\nall molecules listed in the structure file.\nWe cannot show this mapping.')
                return
            except:
                return
            self.parent.data.drug_molecules_linecol = Line3DCollection(seg_all, cmap=plt.get_cmap('viridis_r'))
            self.parent.data.drug_molecules_linecol.set_array(indices_all)
            self.parent.data.drug_molecules_show_mapping = True
            self.parent.data.drug_molecules_mseg = seg_all
            self.parent.data.drug_molecules_mvalues = indices_all
            self.parent.data.drug_molecules_mlabel = 'frame number'
            self.parent.data.drug_molecules_morder = molecule_order

            # clear and redraw MainView
            self.parent.sc.ax.cla()
            self.parent.sc.fig.canvas.draw()
            # draw figure to ensure constrained layout can be computed
            self.parent.sc.fig.set_constrained_layout(True)
            self.parent.sc.updateView(
                True,
                positionMolecules=self.parent.slider_position_molecules,
                positionMembrane=self.parent.slider_position_membrane,
                data=self.parent.data
                )
            self.parent.sc.ax.set_anchor('C')
            self.parent.sc.draw()
            # update toolbar to set current "home"
            self.parent.toolbar.update()

            # update info panel
            for child in self.parent.findChildren(QLabel):
                if child.accessibleName() == 'analysis_unit_label':
                    child.setText('Position')
            self.showCautionMessageBox('Frame number mapping has been calculated.', QMessageBox.Icon.Information, 'Information')

    def mapIntraDistance(self):
        """Calculate the distances between selected atoms."""
        if self.parent.data.top is None:
            self.showCautionMessageBox()
        else:
            if self.parent.show_selection:
                self.showOkCancelMessageBox('Mapping', 'This functionality is not implemented currently.\nPlease reset the view to the original frame selection first.')
                return
            # Select which distances / molecule to plot
            dlg = AnalysisSelectionDialog(
                'drug_molecules',
                intra_all = True,
                intra=True,
                atoms=self.parent.data.atoms,
                atoms_dict=self.parent.data.atoms_dict,
                distinct_mol=self.parent.data.distinct_molecules,
                distinct_mol_t=self.parent.data.distinct_molecule_types,
                drug_mol=self.parent.data.drug_molecules,
                windowtitle='Distance Selection',
                universe=self.parent.data.universe
                )
            dlg.exec()
            selection_strings = dlg.selected_molecules
            # return if aborted by the user
            if len(dlg.selected_molecules) == 0:
                return
            self.resetView(message=False)
            # Calculate distances
            distances = {}
            extremes = {}
            for item in selection_strings:
                selection_list = item.split()
                molecule = selection_list[1] + selection_list[4] # + selection_list[7] + selection_list[16]
                atoms = item.split(';')
                atom1 = atoms[0][:-1]
                atom2 = atoms[1][1:]
                distances[molecule] = analyses.distanceAtoms(self.parent.data.universe, atom1, atom2)
                minval = distances[molecule].min()
                imin = list(np.where(distances[molecule] == distances[molecule].min()))[0]
                maxval = distances[molecule].max()
                imax = list(np.where(distances[molecule] == distances[molecule].max()))[0]
                extremes[molecule] = [imin, minval, imax, maxval]

            # Create visualization
            segments = {}
            molecule_order = []
            distance_all = np.array([])
            seg_all = np.array([]).reshape(0, 2, 3)
            i = 0
            try:
                for molecule, positions in self.parent.data.drug_molecules_positions.items():
                    x = positions[:, 0]
                    y = positions[:, 1]
                    z = positions[:, 2]
                    points = np.array([x, y, z]).transpose().reshape(-1, 1, 3)
                    mol_segments = np.concatenate([points[:-1],points[1:]], axis=1)
                    segments[molecule] = mol_segments
                    self.parent.data.drug_molecules_moffset[molecule] = i
                    i += self.parent.data.max_frame
                    molecule_order.append(molecule)
                    distance_all = np.concatenate((distance_all, distances[molecule]), axis=0)
                    seg_all = np.concatenate((seg_all, mol_segments), axis=0)
            except IndexError:
                self.showCautionMessageBox('It seems your data does not contain information about\nall molecules listed in the structure file.\nWe cannot show this mapping.')
                return
            except:
                return

            self.parent.data.drug_molecules_linecol = Line3DCollection(seg_all, cmap=plt.get_cmap('viridis_r'))
            self.parent.data.drug_molecules_linecol.set_array(distance_all)
            self.parent.data.drug_molecules_show_mapping = True
            self.parent.data.drug_molecules_mseg = seg_all
            self.parent.data.drug_molecules_mvalues = distance_all
            self.parent.data.drug_molecules_mlabel = 'distance [$\mathrm{\AA}$]'
            self.parent.data.drug_molecules_morder = molecule_order
            self.parent.data.extremePos()
            self.parent.data.drug_molecules_text = extremes

            # clear and redraw MainView
            self.parent.sc.ax.cla()
            # draw figure to ensure constrained layout can be computed
            self.parent.sc.fig.canvas.draw()
            self.parent.sc.fig.set_constrained_layout(True)
            self.parent.sc.updateView(
                True,
                positionMolecules=self.parent.slider_position_molecules,
                positionMembrane=self.parent.slider_position_membrane,
                data=self.parent.data
                )
            self.parent.sc.ax.set_anchor('C')
            self.parent.sc.draw()
            # update toolbar to set current "home"
            self.parent.toolbar.update()

            # Update info panel
            unit_text = '['.encode('utf-8') + u'\u00C5'.encode('utf-8') + ']'.encode('utf-8')
            unit_text = unit_text.decode('utf-8')
            unit_text = 'Distance, ' + unit_text
            self.updateInfoPanelMapping(extremes, unit_text)
            self.showCautionMessageBox('Distance mapping has been calculated.', QMessageBox.Icon.Information, 'Information')

    def belowDistance(self):
        """ask for atoms / molecules  to calculate distance, set up BottomView"""
        if self.parent.data.top is None:
            self.showCautionMessageBox()
        else:
            single_view = True
            if self.accessibleName() == 'multiple': #FIX by Beat, removed .sender()
                single_view = False
            dlg = AnalysisSelectionDialog(
                'all',
                intra=True,
                atoms=self.parent.data.atoms,
                atoms_dict=self.parent.data.atoms_dict,
                distinct_mol=self.parent.data.distinct_molecules,
                distinct_mol_t=self.parent.data.distinct_molecule_types,
                drug_mol=self.parent.data.drug_molecules,
                windowtitle='Distance Selection',
                universe=self.parent.data.universe
                )
            dlg.exec()
            selected_molecules = dlg.selected_molecules
            if len(selected_molecules) == 0: # don't do anything, just return
                return
            else:
                self.parent.additional_data_view.setHidden(False) # show bottom view group if not shown already
            # Setup for shared scales
            if not 'distance' in self.parent.bottom_views_scales.keys():
                self.parent.bottom_views_scales['distance'] = [0, 0] # y min, y max
            y_maxs = []
            # Collect data
            below_distances = {}
            colors = {}
            time_steps = {}
            for item in selected_molecules:
                selection_list = item.split()
                molecule = selection_list[1] + selection_list[4] + selection_list[7] + selection_list[16]
                atoms = item.split(';')
                atom1 = atoms[0][:-1]
                atom2 = atoms[1][1:]
                dist = analyses.distanceAtoms(self.parent.data.universe, atom1, atom2)
                # dist
                below_distances[molecule] = dist
                # color
                mol = selection_list[1]+selection_list[4]
                if mol in self.parent.data.drug_molecules_colors.keys():
                    color = self.parent.data.drug_molecules_colors[mol]
                else:
                    cols = returnColors(256, self.parent.data.colormap)
                    if self.parent.data.colormap in colormap.color_map_grey:
                        rand = random.randint(3, 255)
                        color = cols[rand]
                    else:
                        rand = random.randint(0, 255)
                        color = cols[rand]
                if not single_view:
                    n = 5
                    while any((x==color).all() for x in colors.values()):
                        if self.parent.data.colormap in colormap.color_map_grey:
                            color = returnColors(n, self.parent.data.colormap)[random.randint(1,n-1)]
                        else:
                            color = returnColors(n, self.parent.data.colormap)[random.randint(0,n-1)]
                        n += 7
                        if n >= 256:
                            n = 3
                colors[molecule] = color

            unit_text = '['.encode('utf-8') + u'\u00C5'.encode('utf-8') + ']'.encode('utf-8')
            unit_text = unit_text.decode('utf-8')
            unit_mpl='distance [$\mathrm{\AA}$]'

            lines = []
            for mol, distance in below_distances.items():
                line = SingleDataSet(list(range(len(distance)+1))[1:], distance, mol, colors[mol])
                for index in line.y_mini:
                    line.addLabel(index, 'min')
                for index in line.y_maxi:
                    line.addLabel(index, 'max')
                if single_view:
                    data = DataSub()
                    data.setData(u=unit_text, u_mpl=unit_mpl, analysis='Distance')
                    data.slider_range = (1, len(line.x))
                    data.addLine(line)
                    data.calculateExtremes()
                    y_maxs.append(data.max[1])
                    # use unique name
                    data.id = mol
                    for name in self.parent.additional_data_view.bottom_views.keys():
                        if data.id == name:
                            data.id = data.id + '_' + str(len(self.parent.additional_data_view.bottom_views))
                    bv = BottomView(self.parent.additional_data_view, data.id, line.color, data)
                    self.parent.additional_data_view.newView(bv)
                    bv.view.linePlot2D(True, data)
                    bv.setInfoPanel()
                else:
                    lines.append(line)
            if not single_view:
                data = DataSub()
                data.setData(u=unit_text, u_mpl=unit_mpl, analysis='Distance')
                for l in lines:
                    data.addLine(l)
                data.slider_range = (1, len(lines[0].x)) # all lines have same length here
                data.calculateExtremes()
                y_maxs.append(data.max[1])
                data.id = 'multipleDistance'
                # use unique name
                for name in self.parent.additional_data_view.bottom_views.keys():
                    if data.id == name:
                        data.id = data.id + '_' + str(len(self.parent.additional_data_view.bottom_views))
                bv = BottomView(self.parent.additional_data_view, data.id, [0, 0, 0], data)
                self.parent.additional_data_view.newView(bv)
                bv.view.linePlot2D(False, data)
                bv.setInfoPanel()

            current_y_max = self.parent.bottom_views_scales['distance'][1]
            y_maxs.append(current_y_max)
            self.parent.bottom_views_scales['distance'][1] = np.array(y_maxs).max()+0.1
            for key, value in self.parent.additional_data_view.bottom_views.items():
                if value.data.analysis == 'Distance':
                    value.view.ax.set_ylim(0, self.parent.bottom_views_scales['distance'][1])
            self.parent.get_current_window_size()
            if self.parent.width < 980:
                self.parent.resize(720+self.parent.size_info_max,self.parent.height)

    def showDocumentation(self):
        """Show a documentation for this program."""
        dlg = DocumentationDialog()
        dlg.exec()

    def showAbout(self):
        """Show the About Dialog."""
        dlg = AboutDialog()
        dlg.exec()

    def mapSpeed(self):
        """Calculate all variables needed to generate a line collection to show in the main view."""
        if self.parent.data.top is None:
            self.showCautionMessageBox()
        else:
            if self.parent.show_selection:
                self.showOkCancelMessageBox('Mapping', 'This functionality is not implemented currently.')
                return
            if self.parent.data.max_frame == 0:
                self.showCautionMessageBox('You cannot calculate a speed for only one given frame.', QMessageBox.Icon.Information)
                return
            self.resetView(message=False)
            speeds = {}
            extremes = {}
            for molecule, name in self.parent.data.drug_molecules.items():
                _selection_string = 'resname '+name[0]+' and resid '+str(name[1])
                speeds[molecule] = analyses.speedMolecule(self.parent.data.universe, _selection_string)[0]
                minval = speeds[molecule].min()
                imin = list(np.where(speeds[molecule] == speeds[molecule].min()))[0]
                maxval = speeds[molecule].max()
                imax = list(np.where(speeds[molecule] == speeds[molecule].max()))[0]
                extremes[molecule] = [imin, minval, imax, maxval]
             
            segments = {}
            molecule_order = []
            speed_all = np.array([])
            seg_all = np.array([]).reshape(0, 2, 3)
            i = 0
            try:
                for molecule, positions in self.parent.data.drug_molecules_positions.items():
                    x = positions[:, 0]
                    y = positions[:, 1]
                    z = positions[:, 2]
                    points = np.array([x, y, z]).transpose().reshape(-1, 1, 3)
                    mol_segments = np.concatenate([points[:-1],points[1:]], axis=1)
                    segments[molecule] = mol_segments
                    self.parent.data.drug_molecules_moffset[molecule] = i
                    i += self.parent.data.max_frame
                    molecule_order.append(molecule)
                    speed_all = np.concatenate((speed_all, speeds[molecule]), axis=0)
                    seg_all = np.concatenate((seg_all, mol_segments), axis=0)
            except IndexError:
                self.showCautionMessageBox('It seems your data does not contain information about\nall molecules listed in the structure file.\nWe cannot show this mapping.')
                return
            except:
                return
            self.parent.data.drug_molecules_linecol = Line3DCollection(seg_all, cmap=plt.get_cmap('viridis'))
            self.parent.data.drug_molecules_linecol.set_array(speed_all)
            self.parent.data.drug_molecules_show_mapping = True
            self.parent.data.drug_molecules_mseg = seg_all
            self.parent.data.drug_molecules_mvalues = speed_all
            self.parent.data.drug_molecules_mlabel = 'speed [nm/ns]' # speed in angrtrˆm: [$\mathrm{\AA}$/$\mathrm{\mu}$s]'
            self.parent.data.drug_molecules_morder = molecule_order
            self.parent.data.extremePos()
            self.parent.data.drug_molecules_text = extremes

            # clear and redraw MainView
            self.parent.sc.ax.cla()
            # draw figure to ensure constrained layout can be computed
            self.parent.sc.fig.canvas.draw()
            self.parent.sc.fig.set_constrained_layout(True)
            self.parent.sc.updateView(
                True,
                positionMolecules=self.parent.slider_position_molecules,
                positionMembrane=self.parent.slider_position_membrane,
                data=self.parent.data
                )
            self.parent.sc.ax.set_anchor('C')
            self.parent.sc.draw()
            # update toolbar to set current "home"
            self.parent.toolbar.update()

            # Update info panel
            unit_text = '[nm/ns]'
            unit_text = 'Speed, ' + unit_text
            self.updateInfoPanelMapping(extremes, unit_text)
            self.showCautionMessageBox('Speed mapping has been calculated.', QMessageBox.Icon.Information, 'Information')

    def updateInfoPanelMapping(self, extremes, unit_text):
        """Update the Analysis, Unit, min, max info panel."""
        names = extremes.keys()
        for child in self.parent.findChildren(QLabel):
            for name in names:
                if child.accessibleName() == name+'min':
                    text = str('{:.2f}'.format(extremes[name][1]))
                    child.setText(text)
                    child.setHidden(False)
                if child.accessibleName() == name+'minframe':
                    if len(extremes[name][0]) > 1:
                        text = str(extremes[name][0][0]) + ',..'
                    else:
                        text = str(extremes[name][0][0])
                    child.setText(text)
                    child.setHidden(False)
                if child.accessibleName() == name+'max':
                    text = str('{:.2f}'.format(extremes[name][3]))
                    child.setText(text)
                    child.setHidden(False)
                if child.accessibleName() == name+'maxframe':
                    if len(extremes[name][2]) > 1:
                        text = str(extremes[name][2][0]) + ',..'
                    else:
                        text = str(extremes[name][2][0])
                    child.setText(text)
                    child.setHidden(False)
                if child.accessibleName() == 'analysis_unit_label':
                    text = unit_text
                    child.setText(text)
                    child.setHidden(False)
                    child.adjustSize()
                if child.accessibleName() == 'analysis_min_label':
                    child.setHidden(False)
                if child.accessibleName() == 'analysis_max_label':
                    child.setHidden(False)

    def belowSpeed(self):
        """Calculate everything necessary to show speed information for selected molecules/atoms below the main view."""
        if self.parent.data.top is None:
            self.showCautionMessageBox()
        else:
            single_view = True
            if self.accessibleName() == 'multiple': #FIX by Beat, removed .sender()
                single_view = False
            dlg = AnalysisSelectionDialog(
                'drug_molecules',
                atoms=self.parent.data.atoms,
                atoms_dict=self.parent.data.atoms_dict,
                distinct_mol=self.parent.data.distinct_molecules,
                distinct_mol_t=self.parent.data.distinct_molecule_types,
                drug_mol=self.parent.data.drug_molecules,
                windowtitle='Speed Selection',
                universe=self.parent.data.universe
                )
            dlg.exec()
            selected_molecules = dlg.selected_molecules
            if len(selected_molecules) == 0: # don't do anything, just return
                return
            else:
                self.parent.additional_data_view.setHidden(False) # show bottom view group
            # Setup for shared scales
            if not 'speed' in self.parent.bottom_views_scales.keys():
                self.parent.bottom_views_scales['speed'] = [0, 0] # y min, y max
            y_maxs = []
            # Collect data
            below_speeds = {}
            colors = {}
            time_steps = {}
            # calculate speed, setup data, call add view function
            for item in selected_molecules:
                selection_list = item.split()
                if len(selection_list) == 5:
                    molecule = selection_list[1] + selection_list[4]
                else: # len(selection_list) == 8
                    molecule = selection_list[1] + selection_list[4] + selection_list[7]
                analysis = analyses.speedMolecule(self.parent.data.universe, item)
                time_steps[molecule] = analysis[1]
                below_speeds[molecule] = analysis[0]
                color = self.parent.data.drug_molecules_colors[selection_list[1] + selection_list[4]]
                # n is chosen arbitrarily, feel free to change this
                if not single_view:
                    n = 5
                    while any((x==color).all() for x in colors.values()):
                        if self.parent.data.colormap in colormap.color_map_grey:
                            color = returnColors(n, self.parent.data.colormap)[random.randint(1,n-1)]
                        else:
                            color = returnColors(n, self.parent.data.colormap)[random.randint(0,n-1)]
                        n += 7
                        if n >= 256:
                            n = 3
                colors[molecule] = color
            # Labels etc. for the plot (unit_text for Qt text, unit_mpl for matplotlib plots)
            # [Angstrom / u-meter]
            # unit_text = '['.encode('utf-8') + u'\u00C5'.encode('utf-8') + '/'.encode('utf-8') + u'\u00B5'.encode('utf-8') + 'm]'.encode('utf-8')
            # unit_text = unit_text.decode('utf-8')
            unit_text = '[nm / ns]'
            unit_mpl= 'speed [nm / ns]' # 'speed [$\mathrm{\AA}$/$\mathrm{\mu}$s]'
        
            lines = []
            for mol, speed in below_speeds.items():
                line = SingleDataSet(list(range(len(time_steps[mol])+1))[1:], speed, mol, colors[mol], t=time_steps[mol])
                for index in line.y_mini:
                    line.addLabel(index, 'min')
                for index in line.y_maxi:
                    line.addLabel(index, 'max')
                if single_view:
                    data = DataSub()
                    data.setData(u=unit_text, u_mpl=unit_mpl, analysis='Speed')
                    data.slider_range = (1, len(line.t))
                    data.addLine(line)
                    data.calculateExtremes()
                    y_maxs.append(data.max[1])
                    # use unique name
                    data.id = mol
                    for name in self.parent.additional_data_view.bottom_views.keys():
                        if data.id == name:
                            data.id = data.id + '_' + str(len(self.parent.additional_data_view.bottom_views))
                    bv = BottomView(self.parent.additional_data_view, data.id, line.color, data)
                    self.parent.additional_data_view.newView(bv)
                    bv.view.linePlot2D(True, data)
                    bv.setInfoPanel()
                else:
                    lines.append(line)
            if not single_view:
                data = DataSub()
                data.setData(u=unit_text, u_mpl=unit_mpl, analysis='Speed')
                for l in lines:
                    data.addLine(l)
                data.slider_range = (1, len(lines[0].t)) # all lines have same length here
                data.calculateExtremes()
                y_maxs.append(data.max[1])
                data.id = 'multipleSpeed'
                # use unique name
                for name in self.parent.additional_data_view.bottom_views.keys():
                    if data.id == name:
                        data.id = data.id + '_' + str(len(self.parent.additional_data_view.bottom_views))
                bv = BottomView(self.parent.additional_data_view, data.id, [0, 0, 0], data)
                self.parent.additional_data_view.newView(bv)
                bv.view.linePlot2D(False, data)
                bv.setInfoPanel()
            # All speed calculations have '0' as minimum y axis value
            # update all axes from speed views:
            current_y_max = self.parent.bottom_views_scales['speed'][1]
            y_maxs.append(current_y_max)
            self.parent.bottom_views_scales['speed'][1] = np.array(y_maxs).max()+0.1
            for key, value in self.parent.additional_data_view.bottom_views.items():
                if value.data.analysis == 'Speed':
                    value.view.ax.set_ylim(0, self.parent.bottom_views_scales['speed'][1])
            self.parent.get_current_window_size()
            if self.parent.width < 980:
                self.parent.resize(720+self.parent.size_info_max,self.parent.height)
        
    def changeFrames(self):
        """Show open dialogs needed to change the view settings of the main view."""
        if self.parent.data.top is None:
            self.showCautionMessageBox()
        else:
            if self.parent.data.drug_molecules_show_mapping:
                self.resetView()
            dlg = ChangeViewDataDialog(self.parent, 1, self.parent.data.max_frame + 1) # frames start at 1
            dlg.exec()
            if dlg.from_value is None or dlg.to_value is None:
                return
            else:
                # dlg from and to values are not indices! conversion: - 1
                self.parent.data.changeSettings(dlg.from_value-1, dlg.to_value-1)
                self.parent.slider_molecules.setValue(0)
                self.parent.slider_membrane.setValue(0)
                self.parent.slider_molecules.setMinimum(self.parent.data.min_frame)
                self.parent.slider_membrane.setMinimum(self.parent.data.min_frame)
                self.parent.slider_molecules.setMaximum(self.parent.data.max_frame)
                self.parent.slider_membrane.setMaximum(self.parent.data.max_frame)
                self.parent.show_selection = True
                self.parent.selection_from = dlg.from_value
                self.parent.selection_to = dlg.to_value
                self.parent.set_frame.setMinimum(dlg.from_value)
                self.parent.set_frame.setMaximum(dlg.to_value)
                # redraw view
                self.parent.redrawMainView()
                # update linkage, vis
                self.parent.slider_from.setText(str(dlg.from_value))
                self.parent.slider_to.setText(str(dlg.to_value))
                # update overview panel
                self.parent.shown_frames.setText(
                    'shown frames: ' \
                    + '{:,}'.format(self.parent.data.min_frame + dlg.from_value).replace(',',' ') \
                    + ' - ' + '{:,}'.format(self.parent.data.max_frame + dlg.from_value).replace(',',' ')
                    )

    def resetView(self, **kwargs):
        """Reset the view to show original (first) view."""
        if self.parent.data.top is None:
            self.showCautionMessageBox('No data loaded so no reset needed.')
        else:
            if self.parent.show_selection:
                self.parent.data.resetSettings()
                self.parent.slider_molecules.setValue(0)
                self.parent.slider_membrane.setValue(0)
                self.parent.slider_molecules.setMinimum(self.parent.data.min_frame)
                self.parent.slider_membrane.setMinimum(self.parent.data.min_frame)
                self.parent.slider_molecules.setMaximum(self.parent.data.max_frame)
                self.parent.slider_membrane.setMaximum(self.parent.data.max_frame)
                self.parent.show_selection = False
                self.parent.selection_from = None
                self.parent.selection_to = None
                self.parent.set_frame.setMinimum(self.parent.data.min_frame + 1)
                self.parent.set_frame.setMaximum(self.parent.data.max_frame + 1)
                # redraw view
                self.parent.redrawMainView()
                # update linkage, vis
                self.parent.slider_from.setText('1')
                self.parent.slider_to.setText(str(self.parent.data.max_frame + 1))
                # update overview panel
                self.parent.shown_frames.setText('{:,}'.format(self.parent.data.max_frame + 1).replace(',',' ')+' frames shown')
            else:
                message = True
                for key, value in kwargs.items():
                    if key == 'message':
                        message = value
                        reset = True
                if message:
                    reset = self.showOkCancelMessageBox('Reset', 'Are you sure you want to reset the main view?')
                if self.parent.data.drug_molecules_show_mapping:
                    if reset:
                        # reset variables
                        self.parent.data.drug_molecules_show_mapping = False
                        for mol, val in self.parent.data.drug_molecules_show.items():
                            self.parent.data.drug_molecules_show[mol] = True
                        self.parent.data.drug_molecules_mseg = None
                        self.parent.data.drug_molecules_mseg_original = None
                        self.parent.data.drug_molecules_linecol = None
                        self.parent.data.drug_molecules_mcmap = 'viridis_r'
                        self.parent.data.drug_molecules_mvalues = None
                        self.parent.data.drug_molecules_moffset = {}
                        self.parent.data.drug_molecules_morder = []
                        self.parent.data.drug_molecules_mlabel = ''
                        self.parent.data.drug_molecules_text = {}
                        self.parent.data.mextreme = []
                        
                        # reset figure
                        self.parent.sc.colorbar.remove()
                        self.parent.sc.ax.cla()
                        # draw the canvas, so constrained layout can be computed
                        self.parent.sc.fig.canvas.draw()
                        self.parent.sc.fig.set_constrained_layout(True)
                        self.parent.sc.ax.set_anchor('C')
                        self.parent.sc.updateView(
                            False,
                            positionMolecules=self.parent.slider_position_molecules,
                            positionMembrane=self.parent.slider_position_membrane,
                            data=self.parent.data
                            )
                        self.parent.sc.draw()
                        # update toolbar to set current "home"
                        self.parent.toolbar.update()

                        names = self.parent.data.drug_molecules.keys()
                        for child in self.parent.findChildren(QLabel):
                            for name in names:
                                if child.accessibleName() == name+'min':
                                    child.setText('')
                                    child.setHidden(False)
                                if child.accessibleName() == name+'minframe':
                                    child.setText('')
                                    child.setHidden(False)
                                if child.accessibleName() == name+'max':
                                    child.setText('')
                                    child.setHidden(False)
                                if child.accessibleName() == name+'maxframe':
                                    child.setText('')
                                    child.setHidden(False)
                                if child.accessibleName() == 'analysis_unit_label':
                                    child.setText('Analysis, Unit')
                                    child.setHidden(False)
                                if child.accessibleName() == 'analysis_min_label':
                                    child.setHidden(False)
                                if child.accessibleName() == 'analysis_max_label':
                                    child.setHidden(False)
                    else:
                        return
                else:
                    return

    def savelegendMainView(self):
        """Safe the legend of the main view as an image file."""
        if self.parent.data.top is None:
            self.showCautionMessageBox()
        else:
            if self.parent.data.drug_molecules_show_mapping:
                # legend can be printed with printing function
                self.showOkCancelMessageBox(
                    'Save legend of main view',
                    'This view has a legend to be saved with the view.\n'
                    + 'You can save the view and all legends via\n'
                    + 'the printing icon above the view.'
                    )
                return
            else:
                # inform user, that only currently shown molecules are exported.
                all_shown = True
                any_shown = False
                for name, val in self.parent.data.drug_molecules_show.items():
                    all_shown *= val
                    any_shown += val
                if not any_shown:
                    self.showCautionMessageBox(
                        'No molecules shown as line representation.\n'
                        +'No legend to be saved.\n'
                        +'Activate the view of line-representation molecules.',
                        QMessageBox.Icon.Information
                        )
                    return
                save = True
                if not all_shown:
                    save = self.showOkCancelMessageBox(
                        'Save legend of main view',
                        'The legend shows the currently shown\nmolecules of the line representation.'
                        )
                if save:
                    # get file name
                    file_filter = 'Image File (*.png)'
                    response = QFileDialog.getOpenFileName(
                                parent=self.parent,
                                caption='Select a filename',
                                directory=os.getcwd(),
                                filter=file_filter,
                                initialFilter='Image File (*.png)'
                                )
                    filename = response[0]

                    # generate legend for line molecules
                    colors = []
                    names = []
                    for name, col in self.parent.data.drug_molecules_colors.items():
                        if self.parent.data.drug_molecules_show[name]:
                            names.append(name)
                            colors.append(col)
                        else:
                            # molecule is not shown currently
                            pass
                    f = lambda m, c: plt.plot([],[], color=c, ls='none', marker=m)[0]
                    handles = [f('_', colors[i]) for i in range(len(colors))]
                    legend = plt.legend(handles, names, loc=3, framealpha=1, frameon=False)
                    fig = legend.figure
                    fig.canvas.draw()
                    bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())

                    # save legend
                    fig.savefig(filename, dpi=600, bbox_inches=bbox)
                else:
                    # user declined to save legend
                    return

    def showOkCancelMessageBox(self, title, message):
        """Show a message box which asks the user for consent."""
        reply = QMessageBox.question(
            self.parent,
            title,
            message,
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes
            )
        if reply == QMessageBox.StandardButton.Yes:
            return True
        else:
            return False

    def showCautionMessageBox(
            self,
            text='You need to select data to be shown first.\n Open files first.',
            icon=QMessageBox.Icon.Warning, title='Caution!'
            ):
        """Show message box informing the user about a warning.

        The title of this QMessageBox is 'Caution!', the default
        text asks the user to open files first.

        Parameters:
        ------
        text (str):
            The text to be presented by this warning message box.
            Default: 'You need to select data to be shown first.
            Open files first.'
        """
        message_box = QMessageBox()
        message_box.setWindowTitle(title)
        message_box.setText(text)
        message_box.setIcon(icon)
        message_box.setStandardButtons(QMessageBox.StandardButton.Ok)
        self.buttonBox = message_box.exec()


