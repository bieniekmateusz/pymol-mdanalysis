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
- storing a graph, like RMSD, this means saving the .py analysis script,
    along the data, and the graphing .py
- retrieve the list of graphs (which will be handy for the creation of a menu), and updating the list of graphs,
    understanding the mapping between the label and the graphs,
    * How does it retrieve the graphs? It is given a set of systems, with their selections (list of Atom IDs),
- allows to load all the graphs from the hard drive
    (this could use the .mse which is a MDAnalysis session like .pse)

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
import re
import hashlib


class GraphManager():
    # fixme - should in the configuration file
    TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'plotting_templates')
    PLOTS_DIR = os.path.join(os.path.expanduser('~'), '.pymol', 'plotting')

    # these map directly to the directory names in the plotting directory
    # fixme - convert this into the ENUM type
    GRAPH_TYPES = {'RMSD', 'RMSF'}

    # use the same tree structure: RMSD-HashFilePathFileName-HashSelectedAtoms
    Graphs = {}

    @staticmethod
    def get_hash_from_filepath(full_filepath):
        """
        Extract the alphanumeric characters from the path and get the hash from the 0-9 a-z A-Z
        :param filename:
        :return:
        """

        alphanumeric_filepath = re.sub(r'\W+', '', full_filepath)

        return hashlib.sha1(alphanumeric_filepath.encode()).hexdigest()


    @staticmethod
    def get_hash_from_selection(atom_group):
        """
        Take all atom IDs and convert them into a string and then extract a hash from it.
        This way, when the label changes, the graphs are not affected.
        It is also reversible: we can find the graph having an atom_group
        :param atom_group: atom group
        :return: A hash created from the IDs
        """

        return hashlib.sha1(atom_group._pymol_used_selection.encode()).hexdigest()


    @staticmethod
    def find_graphs(atom_groups):
        """

        :param atom_groups: The atoms
        :return:
        """
        pass


    @staticmethod
    def save_graph(atom_group, rmsd_data, category):
        """

        :param label: PyMOL label
        :param category: a str, e.g. 'rmsd', 'rmsf'
        :param atom_group: used for the filename information
        :param rmsd_data: rmsd results from MDAnalysis
        :return:
        """

        # check if the rmsd directory exists
        filepath_hash = GraphManager.get_hash_from_filepath(atom_group.universe.filename)
        selection_hash = GraphManager.get_hash_from_selection(atom_group)

        datafile_dir = os.path.join(GraphManager.PLOTS_DIR, category, filepath_hash, selection_hash)
        datafile_name = os.path.join(GraphManager.PLOTS_DIR, category, filepath_hash, selection_hash, 'graph.dat')
        if not os.path.isdir(datafile_dir):
            os.makedirs(datafile_dir)

        # save the data in a file
        numpy.savetxt(datafile_name, rmsd_data)

        # copy the plotting file from the templates
        template_rmsd_plotter = os.path.join(GraphManager.TEMPLATE_DIR, '%s.py' % category)
        plotter_filename = os.path.join(GraphManager.PLOTS_DIR, category, filepath_hash, selection_hash, 'graph.py')
        shutil.copyfile(template_rmsd_plotter, plotter_filename)


    @staticmethod
    def plot_graph(atom_group, category):
        # use the basic template to visualise the results

        filepath_hash = GraphManager.get_hash_from_filepath(atom_group.universe.filename)
        selection_hash = GraphManager.get_hash_from_selection(atom_group)
        graph_dir = os.path.join(GraphManager.PLOTS_DIR, category, filepath_hash, selection_hash)

        sys.path.append(graph_dir)
        # import the saved rmsd plotter
        plotter = importlib.import_module('graph')
        sys.path.remove(graph_dir)

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

