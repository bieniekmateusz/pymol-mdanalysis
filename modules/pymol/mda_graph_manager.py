#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-*
#-*
#Z* -------------------------------------------------------------------

"""
This class manages the graphs produced by the MDAnalysis.

Possible responsibilities to consider:
- computing graphs, such as RMSD
- storing a graph, like RMSD, this means saving the .py analysis script,
    along the data, and the graphing .py
- updating the index of graphs, understanding the mapping between the label and the graphs,
- allows to load all the graphs from the hard drive
    (this could use the .mse which is a MDAnalysis session like .pse)
- retrieve the list of graphs (which will be handy for the creation of a menu)

The Graph_Manager will help mange the paths in the system.


Lower priority improvement:
- there should be a central way to obtain the configuration, rather than hard coding it here,
    ie: the deafult paths for the storing graph data

"""

import os, sys
import json
import shutil
import numpy
import importlib
from matplotlib.widgets import SpanSelector
from . import cmd



class GraphManager():
    # fixme - should in the configuration file
    TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'plotting_templates')
    PLOTS_DIR = os.path.join(os.path.expanduser('~'), '.pymol', 'plotting', 'plots')


    @staticmethod
    def save_graph(label, graph_type, atom_group, R):
        """

        :param label: PyMOL label
        :param graph_type: a str, e.g. 'rmsd'
        :param atom_group:
        :param R:
        :return:
        """

        # check if the rmsd directory exists
        rmsd_dir = os.path.join(GraphManager.PLOTS_DIR, graph_type)
        if not os.path.isdir(rmsd_dir):
            os.makedirs(rmsd_dir)

        # save the data in a file
        universe_filename = os.path.splitext(os.path.basename(atom_group.universe.filename))[0]
        datafile_name = os.path.join(GraphManager.PLOTS_DIR, 'rmsd/%s_%s.dat' % (universe_filename, label))
        numpy.savetxt(datafile_name, R.rmsd)

        # copy the plotting file from the templates
        template_rmsd_plotter = os.path.join(GraphManager.TEMPLATE_DIR, 'rmsd.py')
        plotter_filename = os.path.join(GraphManager.PLOTS_DIR, 'rmsd/%s_%s.py' % (universe_filename, label))
        shutil.copyfile(template_rmsd_plotter, plotter_filename)

        return universe_filename


    @staticmethod
    def plot_graph(label, graph_type, universe_filename):
        # use the basic template to visualise the results
        rmsd_dir = os.path.join(GraphManager.PLOTS_DIR, 'rmsd')
        sys.path.append(rmsd_dir)
        # import the saved rmsd plotter
        plotter = importlib.import_module('%s_%s' % (universe_filename, label))
        sys.path.remove(rmsd_dir)

        # attach the interactive features to the functions
        def onclick(event):
            # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            #       ('double' if event.dblclick else 'single', event.button,
            #        event.x, event.y, event.xdata, event.ydata))

            # if the click happened outside of the graph area
            if not type(event.xdata) is numpy.float64:
                return

            # find the nearest point on x line
            closest_index = numpy.searchsorted(plotter.pymol_x_axis, (event.xdata,))
            # update the frame on the screen
            cmd.frame(closest_index)

        plotter.fig.canvas.mpl_connect('button_press_event', onclick)

        def onselect(xmin, xmax):
            # print(xmin, xmax)
            # update the data
            indmin, indmax = numpy.searchsorted(plotter.pymol_x_axis, (xmin, xmax))
            indmax = min(len(plotter.pymol_x_axis) - 1, indmax)
            # print(indmin, indmax)
            # print(plotter.pymol_y_axis[indmin:indmax])

            plotter.pymol_hist_ax.cla()
            plotter.pymol_hist_ax.hist(plotter.pymol_y_axis[indmin:indmax])

        # FIXME - Find a better way than a global variable
        # We have to ensure the object survives when the method is left.
        # Otherwise the class.methods will not be available.
        # See weak references and the Note in https://matplotlib.org/users/event_handling.html
        global G_MATPLOTLIB_CALLBACK_SPAN
        G_MATPLOTLIB_CALLBACK_SPAN = SpanSelector(plotter.pymol_plot_ax, onselect, 'horizontal', useblit=True,
                                                  rectprops=dict(alpha=0.5, facecolor='red'),
                                                  button=1)

