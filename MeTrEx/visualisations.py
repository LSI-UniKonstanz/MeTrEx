"""This is the visualisations module. 

The visualisations module contains different classes for the 
visualization of objects from the data module. 

Classes:
------
MainPlotCanvas:
    Canvas for the main plot of the MeTrEx GUI.
SubPlotCanvas:
    Canvas for different plots used in MeTrEx. The plots can be used
    with dialogs, BottomViews, SubWindows and other.
"""

__version__ = '1.0'
__author__ = 'Christiane Rohse'

from PyQt6.QtCore import pyqtSignal
import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from vis_dialogs import AnnotationDialog
import matplotlib.pyplot as plt

# adjustments of parameters for saving figures
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['savefig.bbox'] = 'tight'


class MainPlotCanvas(FigureCanvasQTAgg):
    """Show the Data.

    Parameters:
    ------
    width (int):
        Width of figure canvas in inches. Default: 0
    height (int):
        Height of the figure canvas in inches. Default: 0

    Methods:
    -----
    updateView() :
        Fill the axes of this view, draw a colorbar if selected.
    """
    def __init__(self, width=0, height=0, dpi=0): #, **kwargs
        """Set up the canvas and call updateView to show data."""
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig.set_constrained_layout(True)
        # tight layout is horrible, when pressing home button multiple 
        # times, therefore only use contrained layout
        self.ax = self.fig.add_subplot(111, projection='3d')
        # if you want to see where the figure and where the axes are
        # uncomment the following:
        #self.fig.set_facecolor('yellow')
        # set background 'invisible' / white
        self.ax.xaxis.pane.fill = False
        self.ax.yaxis.pane.fill = False
        self.ax.zaxis.pane.fill = False
        # white background, so white is 'invisible'
        self.ax.xaxis.pane.set_edgecolor('w')
        self.ax.yaxis.pane.set_edgecolor('w')
        self.ax.zaxis.pane.set_edgecolor('w')
        #self.ax.grid(False)
        self.colorbar = None
        super(MainPlotCanvas, self).__init__(self.fig)

    def updateView(self, mapping, **kwargs):
        """Update the main view according to the Data object stored in
        self.data and the slider position.
        
        """
        self.lines = {}
        self.spheres = {}
        self.leaflets = {}
        self.extreme_labels = {}
        self.line_labels = {}
        self.show_legend = False

        for key, value in kwargs.items():
            if key == 'positionMolecules':
                self.sphere_position = value
            if key == 'positionMembrane':
                self.membrane_index = value
            if key == 'data': # value = self.parent.data
                # PLOT MAPPING
                if mapping:
                    for name, val in value.drug_molecules_positions.items():
                        if value.drug_molecules_show[name]:
                            positions = value.drug_molecules_positions[name]
                            xpos = positions[:,0]
                            ypos = positions[:,1]
                            zpos = positions[:,2]
                            self.line_labels[name] = self.ax.text(xpos[0], ypos[0], zpos[0], name, None)
                            self.spheres[name] = self.ax.scatter(
                                xpos[self.sphere_position], 
                                ypos[self.sphere_position], 
                                zpos[self.sphere_position], 
                                color='fuchsia', 
                                s=50,
                                label='_nolegend_'
                                )
                            self.extreme_labels[name] = []
                            try:
                                j = 0
                                n = len(value.drug_molecules_text[name][2])
                                if n > 0:
                                    self.show_legend = True
                                while j < n:
                                    ind = value.drug_molecules_text[name][2][j]
                                    # Change here, if you want to have min and max text labels
                                    #txt = 'max: '+str(int(value.drug_molecules_text[name][3]))
                                    #label = self.ax.text(xpos[ind], ypos[ind], zpos[ind], txt, None)
                                    #label.set_fontsize('small')
                                    marker = self.ax.scatter(xpos[ind], ypos[ind], zpos[ind], marker='P',color='fuchsia', s=20, label='max')
                                    #self.extreme_labels[name].append([label, marker])
                                    self.extreme_labels[name].append([marker])
                                    j += 1
                                i=0
                                n = len(value.drug_molecules_text[name][0])
                                while i < n:
                                    ind = value.drug_molecules_text[name][0][i]
                                    # Change here, if you want to have min and max text labels
                                    #txt = 'min: '+str(int(value.drug_molecules_text[name][1]))
                                    #label = self.ax.text(xpos[ind], ypos[ind], zpos[ind], txt, None)
                                    #label.set_fontsize('small')
                                    marker = self.ax.scatter(xpos[ind], ypos[ind], zpos[ind], marker='d', color='fuchsia', s=20, label='min')
                                    #self.extreme_labels[name].append([label, marker])
                                    self.extreme_labels[name].append([marker])
                                    i += 1
                            except: 
                                # no extreme values available (e.g. map position)
                                pass
                    if self.show_legend:
                        self.legend_labels = {}
                        self.ax.legend(['max', 'min'], bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
                    else:
                        pass

                    try:
                        self.ax.set_xlim(value.mextreme[0], value.mextreme[1])
                        self.ax.set_ylim(value.mextreme[2], value.mextreme[3])
                        self.ax.set_zlim(value.mextreme[4], value.mextreme[5])
                    except:
                        # no extreme values available (e.g. map position)
                        pass
                    self.ax.set_xlabel(value.label_x)
                    self.ax.set_ylabel(value.label_y)
                    self.ax.set_zlabel(value.label_z)
                    self.linecoll = self.ax.add_collection3d(value.drug_molecules_linecol)
                    self.colorbar = self.fig.colorbar(self.linecoll, ax=self.ax, shrink=0.7, pad=0.13, location='right', label=value.drug_molecules_mlabel) # 'distance [$\mathrm{\AA}$]'
                # PLOT LINE REPRESENTATION ORIGINAL
                else:
                    for name, val in value.drug_molecules_positions.items():
                        if value.drug_molecules_show[name]:
                            try:
                                positions = value.drug_molecules_positions[name]
                                xpos = positions[:,0]
                                ypos = positions[:,1]
                                zpos = positions[:,2]
                                line = self.ax.plot(xpos, ypos, zpos, color=value.drug_molecules_colors[name], label=name)
                                text = self.ax.text(xpos[0], ypos[0], zpos[0], name, None)
                                self.lines[name] = [line, text]
                                self.spheres[name] = self.ax.scatter(
                                    xpos[self.sphere_position], 
                                    ypos[self.sphere_position], 
                                    zpos[self.sphere_position], 
                                    color='fuchsia', 
                                    s=50
                                    )
                            except: # there are atoms without selection of position (no 'C'-Atom)
                                pass
                        else:
                            pass
                # PLOT SURFACE
                # colors: selected from 'imola' colormap
                if value.lipid_molecules_show[0]:
                    if value.lipid_molecules_show_abstraction:
                        self.leaflets['leaflet0'] = self.ax.plot_surface(
                            value.lipid_molecules_surface['leaflet0'][self.membrane_index][0], 
                            value.lipid_molecules_surface['leaflet0'][self.membrane_index][1], 
                            value.lipid_molecules_surface['leaflet0'][self.membrane_index][2], 
                            color=value.leaflet_colors['leaflet0'], #[0.69114, 0.87318, 0.4097], 
                            antialiased=True, 
                            alpha=0.2
                            )
                    else:
                         self.leaflets['leaflet0'] = self.ax.plot_trisurf(
                             value.lipid_molecules_positions['leaflet0'][self.membrane_index][:,0], 
                             value.lipid_molecules_positions['leaflet0'][self.membrane_index][:,1], 
                             value.lipid_molecules_positions['leaflet0'][self.membrane_index][:,2], 
                             color=value.leaflet_colors['leaflet0'], #[0.69114, 0.87318, 0.4097], 
                             linewidth=0.2, 
                             antialiased=True, 
                             alpha=0.2, 
                             edgecolor='grey'
                             )
                if value.lipid_molecules_show[1]:
                    if value.lipid_molecules_show_abstraction:
                        self.leaflets['leaflet1'] = self.ax.plot_surface(
                            value.lipid_molecules_surface['leaflet1'][self.membrane_index][0], 
                            value.lipid_molecules_surface['leaflet1'][self.membrane_index][1], 
                            value.lipid_molecules_surface['leaflet1'][self.membrane_index][2], 
                            color=value.leaflet_colors['leaflet1'], #[0.54228, 0.74134, 0.44047], 
                            antialiased=True, 
                            alpha=0.2
                            )
                    else:
                        self.leaflets['leaflet1'] = self.ax.plot_trisurf(
                            value.lipid_molecules_positions['leaflet1'][self.membrane_index][:,0], 
                            value.lipid_molecules_positions['leaflet1'][self.membrane_index][:,1], 
                            value.lipid_molecules_positions['leaflet1'][self.membrane_index][:,2], 
                            color=value.leaflet_colors['leaflet1'], #[0.54228, 0.74134, 0.44047], 
                            linewidth=0.2, 
                            antialiased=True, 
                            alpha=0.2, 
                            edgecolor='grey'
                            )
                self.ax.set_xlabel(value.label_x)
                self.ax.set_ylabel(value.label_y)
                self.ax.set_zlabel(value.label_z)
        self.draw()


class SubPlotCanvas(FigureCanvasQTAgg):
    """Show given Data in a SubWindow / any other position.
    
    Parameters:
    -----
    vline_changed(pyqtSignal):
        Emitted signal when the position vline of a 2D line plot
        has changed.

    Methods:
    -----
    setAccessibleName():
        Set the name of the plot canvas.
    accessibleName():
        Get the name of the plot canvas.
    linePlot2D():
        Show lines in a 2D plot.
    addMinMaxLabels2DLine():
        Add labels for the minimum and maximum values of a line
        plot. Only applicable for one line shown.
    scatterPlot2D():
        Show data as scatter plot in 2D.
    scatterPlot3D():
        Show data as scatter plot in 3D.
    updatePointsets():
        Change the name of a SingleData object and recalculate the view.
    updateSpheresScatterPlot():
        Update the position of the spheres in the scatter plot. 
    """
    vline_changed = pyqtSignal(int, name='vlineChanged')

    def __init__(self, parent, width=0, height=0, dpi=0):
        """Set up an empty plotting canvas to be filled with the 
        corresponting analysis visualization.
        
        """
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = None
        super(SubPlotCanvas, self).__init__(self.fig)
        self.dim = None # 2d or 3d

    def setAccessibleName(self, name):
        """Set the name attribute."""
        self.name = name

    def accessibleName(self):
        """Return the name attribute."""
        return self.name

    def linePlot2D(self, single, data):
        """Draw a single-line line plot with a vertical line. On its
        top you can read the current x and y values optionally 
        (if single line selected). 
        
        Parameter: 
        ------
        single (bool): Display a single line (True) or 
            multiple lines(false).
        data: DataSun2D object holding all needed data and information.

        Events: 
        ------
        on press: change slider (vline) position variable
        on release: redraw slider (vline) position
        on hover: show values (optionally: if single line selected)
        on leave: hide values (optionally: if single line selected)
        """
        self.single = single # bool: True if only one line
        self.ax = self.fig.add_subplot(111)
        self.fig.set_constrained_layout(True)
        self.data = data
        self.change_vline = False
        self.annotations = {}
        self.lines = {}
        for l in self.data.lines:
            x = l.x
            y = l.y
            c = l.color
            self.lines[l.selection_string], = self.ax.plot(x, y, c=c)
            self.lines[l.selection_string].set_label(l.selection_string)
        if self.single:
            self.addMinMaxLabels2DLine()
            self.ax.axvline(self.data.lines[0].y_mini[0]+1, c='darkgrey', linestyle='dashed', linewidth=0.5) 
            self.ax.axvline(self.data.lines[0].y_maxi[0]+1, c='darkgrey', linestyle='dashed', linewidth=0.5) 

        # vline is positioned at x=1
        self.vline = self.ax.axvline(1, c='fuchsia', linewidth=2.5)

        self.ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
        self.ax.set_xmargin(0.001) # otherwise you cannot select vline and drag it to different position
        self.ax.set_ylabel(self.data.unit_matplotlib)
        #self.ax.set_xlabel('Frame')

        if self.single:
            vline_info_text = ' frame: ' + str(1) + '\n' + ' ' + str('{:.2f}'.format(self.data.lines[0].y[0]))
        else: 
            vline_info_text = ' frame: ' + str(1)
        self.vline_info = self.ax.text(0, 1.05, vline_info_text, transform=self.ax.transAxes)

        #self.annotation_test = self.ax.annotate('user annotation', xy=(0, self.extremes[1]))

        def onPress(event): 
            try:
                if self.vline.get_xdata()[0] == int(event.xdata): 
                    self.change_vline = True
                else:
                    self.change_vline = False
            except: # you clicked next to figure = nothing to do
                pass
            if event.dblclick:
                # Here you can implement anything triggered by a double click. 
                # This could be adding an anntoation.
                dlg = AnnotationDialog('', {})
                # uncomment the following line and implement the 
                # Annotation Dialog as well as the functionality to add
                # and remove annotations. Keep in mind, to check other
                # view related functionalities for changes as well.
                #dlg.exec()
                if dlg.annotations_changed:
                    # remove old annotations
                    # draw new annotations
                    # update annotation-list?
                    self.draw()
                else:
                    return

        def onRelease(event):
            if self.change_vline: 
                x=int(event.xdata) # = index + 1
                self.vline_changed.emit(x)
                self.vline.remove()
                self.vline = self.ax.axvline(x, c='fuchsia', linewidth=2.5)
                if self.single:
                    vline_info_text = ' frame: ' + str(x) + '\n' + ' ' + str('{:.2f}'.format(self.data.lines[0].y[x-1]))
                else: 
                    vline_info_text = ' frame: ' + str(x)
                self.vline_info.remove()
                x_pos = float((x - 1) / self.ax.get_xbound()[1])
                self.vline_info = self.ax.text(x_pos, 1.05, vline_info_text, transform=self.ax.transAxes)
                self.draw()
            else: 
                return

        def onLeaveAxes(event):
            try: 
                self.event_label.remove()
            except:
                pass
            self.draw()

        def onHover(event):
            try: 
                self.event_label.remove()
            except:
                pass
            if event.inaxes: 
                x = int(event.xdata) # = index + 1
                if x < self.data.lines[0].x[0]:
                    return
                if self.single:
                    label_text = str(x)+': '+ str('{:.2f}'.format(self.data.lines[0].y[x-1]))
                    self.event_label = self.ax.text(x, int(self.ax.get_ylim()[1]*0.9), label_text)
                    self.draw()
            
        connect_enter = self.fig.canvas.mpl_connect('button_press_event', onPress)
        connect_leave = self.fig.canvas.mpl_connect('button_release_event', onRelease)
        connect_hover = self.fig.canvas.mpl_connect('motion_notify_event', onHover)
        connect_leave_axes = self.fig.canvas.mpl_connect('axes_leave_event', onLeaveAxes)

    def addMinMaxLabels2DLine(self):
        """Set min and max labels for LinePlot2D plots."""
        self.label_min = self.ax.text(self.data.lines[0].y_mini[0]+1, -0.1, 'min', bbox=dict(boxstyle='square,pad=0.0', edgecolor='white'))
        self.label_min.set_fontsize('small')
        self.label_min.set_backgroundcolor('white')
        self.label_max = self.ax.text(self.data.lines[0].y_maxi[0]+1, -0.1, 'max', bbox=dict(boxstyle='square,pad=0.0', edgecolor='white'))
        self.label_max.set_fontsize('small')
        self.label_max.set_backgroundcolor('white')

    def scatterPlot2D(self, data, sphere=True, legend=True):
        """Draw a 2D scatter plot with a position indicating sphere 
        per data set.

        Parameter: 
        ------
        data (DataSub): Data object that hold all relevant 
            information.
        """
        self.dim = '2d'
        self.ax = self.fig.add_subplot(111)
        self.fig.set_constrained_layout(True)
        self.data = data
        self.s = sphere
        self.l = legend
        self.spheres = {}
        self.pointsets = {}
        self.a = 1.0
        self.legend = None

        # pointset, sphere
        for name, ds in self.data.pointsets.items():
            self.pointsets[ds.selection_string] = self.ax.scatter(ds.x, ds.y, color=ds.color, label=ds.selection_string, alpha=self.a) # c or color! 
            if sphere:
                self.spheres[ds.selection_string] = self.ax.scatter(ds.x[0], ds.y[0], color='fuchsia', s=50)
        # axis labels
        self.ax.set_xlabel(self.data.unit_mpl_x)
        self.ax.set_ylabel(self.data.unit_mpl_y)

        # legend
        if legend:
            self.legend = self.ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.0)

        # Further: picker to get positions, maybe to add annotations. 
        # Implement a method here, as a sub-method.
        # You also need to add the argument picker=True to your scatter
        # plot (e.g. ax.scatter(x, y, picker=True) and connect the 
        # picking event to the figure 
        # (e.g. self.fig.canvas.mpl_connect('pick_event', onPick)).
        # def onPick(event):
        #    # ind = index
        #    ind = event.ind 
        #    # x, y are the x and y values from the scatter plot
        #    print('selected index and values:', ind, x[ind], y[ind]) 
        #    self.ax.redraw_in_frame()

    def scatterPlot3D(self, data, sphere=True, legend=True):
        """Draws a scatter plot in 3D.
        
        Parameter:
        ------
        data (DataSub): DataSub object
        """
        self.fig.set_constrained_layout(True)
        self.s = sphere
        self.l = legend
        self.dim = '3d'
        self.ax = self.fig.add_subplot(projection='3d')
        # white background
        self.ax.xaxis.pane.fill = False
        self.ax.yaxis.pane.fill = False
        self.ax.zaxis.pane.fill = False
        self.ax.xaxis.pane.set_edgecolor('w')
        self.ax.yaxis.pane.set_edgecolor('w')
        self.ax.zaxis.pane.set_edgecolor('w')
        #self.fig.size(5, 5) # TODO: does this make sense here?! 

        self.data = data
        self.spheres = {}
        self.pointsets = {}
        self.legend = None
        self.a = 0.3

        self.ax.set_xlabel(data.unit_mpl_x)
        self.ax.set_ylabel(data.unit_mpl_y)
        self.ax.set_zlabel(data.unit_mpl_z)

        for name, ds in self.data.pointsets.items():
            self.pointsets[ds.selection_string] = self.ax.scatter(ds.x, ds.y, ds.z, color=ds.color, alpha=self.a, label=ds.selection_string)
            if sphere:
                self.spheres[ds.selection_string] = self.ax.scatter(ds.x[0], ds.y[0], ds.z[0], color='fuchsia', s=80)
        if legend:
            self.legend = self.ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.0)

    def updatePointsets(self):
        """Change the name of a SingleData object."""
        old_ps = list(self.pointsets.keys())
        for scat, reference in self.pointsets.items():
            reference.remove()
        for key in old_ps:
            self.pointsets.pop(key)
        old_spheres = list(self.spheres.keys())
        for sphere, reference in self.spheres.items():
            reference.remove()
        for key in old_spheres:
            self.spheres.pop(key)
        for name, ds in self.data.pointsets.items():
            if self.dim == '3d':
                self.pointsets[ds.selection_string] = self.ax.scatter(ds.x, ds.y, ds.z, color=ds.color, alpha=self.a, label=ds.selection_string)
                if self.s:
                    self.spheres[ds.selection_string] = self.ax.scatter(ds.x[0], ds.y[0], ds.z[0], color='fuchsia', s=80)
            else:
                self.pointsets[ds.selection_string] = self.ax.scatter(ds.x, ds.y, color=ds.color, label=ds.selection_string, alpha=self.a) # c or color! 
                if self.s:
                    self.spheres[ds.selection_string] = self.ax.scatter(ds.x[0], ds.y[0], color='fuchsia', s=50)
        if self.l:
            self.legend = self.ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.0)
        self.draw()

    def updateSpheresScatterPlot(self, pos, spheres=True):
        """Update the spheres positions of 2d and 3d scatter plots."""
        if spheres:
            if self.dim == '3d':
                for sphere, reference in self.spheres.items():
                    reference.remove()
                for name, ds in self.data.pointsets.items():
                    self.spheres[name] = self.ax.scatter(ds.x[pos], ds.y[pos], ds.z[pos], color='fuchsia', s=80)
            elif self.dim == '2d':
                for sphere, reference in self.spheres.items():
                    reference.remove()
                for name, ds in self.data.pointsets.items():
                    self.spheres[name] = self.ax.scatter(ds.x[pos], ds.y[pos], color='fuchsia', s=80)
        else: 
            return


