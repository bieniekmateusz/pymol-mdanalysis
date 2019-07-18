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

from .mdanalysis_manager import MDAnalysisManager
from enum import Enum, auto


class GraphManager():
    # fixme - should in the configuration file
    TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'plotting_templates')
    PLOTS_DIR = os.path.join(os.path.expanduser('~'), '.pymol', 'plotting')

    # these map directly to the directory names in the plotting directory
    class GRAPH_TYPES(Enum):
        # we ignore the "value" and use only the "name"
        RMSD = auto(),
        RMSF = auto()


    # use the same tree structure: RMSD-HashFilePathFileName-HashSelectedAtoms
    Graphs = {}

    @staticmethod
    def _get_hash_from_filepath(full_filepath):
        """
        Extract the alphanumeric characters from the path and get the hash from the 0-9 a-z A-Z
        :param filename:
        :return:
        """

        alphanumeric_filepath = re.sub(r'\W+', '', full_filepath)

        return hashlib.sha1(alphanumeric_filepath.encode()).hexdigest()


    @staticmethod
    def _get_hash_from_selection(selection):
        """
        Take all atom IDs and convert them into a string and then extract a hash from it.
        This way, when the label changes, the graphs are not affected.
        It is also reversible: we can find the graph having an atom_group
        :param atom_group: atom group
        :return: A hash created from the IDs
        """

        return hashlib.sha1(selection.encode()).hexdigest()


    @staticmethod
    def find_graphs(systems):
        """
        Use the sha1 hashing to figure out what graphs we have for each label/selection
        :param systems: The system from MDAnalysis_Manager. They should contain both, 1) atom_groups which provides
        access to the filename and path, and 2) selection
        """
        found_graphs = {}
        for label, dict in systems.items():
            atom_group = dict['system']
            selection = dict['selection']

            # check if the system has a directory for this filepath
            filepath_hash = GraphManager._get_hash_from_filepath(atom_group.universe.filename)
            filepath_hash_dir = os.path.join(GraphManager.PLOTS_DIR, filepath_hash)
            if not os.path.isdir(filepath_hash_dir):
                continue

            # check if there is any graph for this selection
            selection_hash = GraphManager._get_hash_from_selection(selection)
            selection_hash_dir = os.path.join(filepath_hash_dir, selection_hash)
            if not os.path.isdir(selection_hash_dir):
                continue

            for graph_type in GraphManager.GRAPH_TYPES:
                # check if there is a directory
                graph_dir = os.path.join(selection_hash_dir, graph_type.name)
                if not os.path.isdir(graph_dir):
                    continue

                # there should be the graph.py and other files in the directory
                # fixme - should be constants not strings
                graph_py = os.path.join(graph_dir, 'graph.py')
                graph_dat = os.path.join(graph_dir, 'graph.dat')
                if not os.path.isfile(graph_py) or not os.path.isfile(graph_dat):
                    continue

                # add the graphs to the found list
                if label in found_graphs:
                    found_graphs[label][graph_type.name] = graph_dir
                else:
                    found_graphs[label] = {graph_type.name: graph_dir}
        return found_graphs


    @staticmethod
    def update_menu():
        # fixme - the graph_manager should not be responsible for updating Qt GUI
        # fixme - ideally, the graph manager would allow the GUI to subscribe for a notification
        #   so that each time there is a new graph, the GUI can add it to the menu
        found_graphs = GraphManager.find_graphs(MDAnalysisManager.Systems)
        # generate the new graphs in the menus directly
        menu_bar = cmd.gui.get_qtwindow().menuWidget()

        # add the menu only if it has not been added before
        import PyQt5.QtWidgets
        plots_menu = menu_bar.findChild(PyQt5.QtWidgets.QMenu, name='Plots')
        if not plots_menu and found_graphs:
            # create the Plots menu
            plots_menu = menu_bar.addMenu('Plots')
            plots_menu.setObjectName('Plots')

        # create a menu plotting action for each existing graph
        for label, graph_types in found_graphs.items():
            # create the label menu only if it does not exist
            label_menu = plots_menu.findChild(PyQt5.QtWidgets.QMenu, name=label)
            if not label_menu:
                label_menu = plots_menu.addMenu(label)
                label_menu.setObjectName(label)

            for graph_type, graph_dir in graph_types.items():
                # ignore if the graph item already exists in the menu
                if label_menu.findChild(PyQt5.QtWidgets.QAction, name=graph_type):
                    continue

                def create_load_graph(label, graph_type):
                    def load_graph():
                        GraphManager.plot_graph(label, graph_type)
                    return load_graph
                action = label_menu.addAction(graph_type, create_load_graph(label, graph_type))
                action.setObjectName(graph_type)

    @staticmethod
    def _add_menu_item(label, category):
        """
        When the user generates a new graph, this function adds it to the GUI menu
        """
        menu_bar = cmd.gui.get_qtwindow().menuWidget()

        # add the menu only if it has not been added before
        import PyQt5.QtWidgets
        plots_menu = menu_bar.findChild(PyQt5.QtWidgets.QMenu, name='Plots')
        if not plots_menu:
            # create the Plots menu
            plots_menu = menu_bar.addMenu('Plots')
            plots_menu.setObjectName('Plots')

        # create the label menu only if it does not exist
        label_menu = plots_menu.findChild(PyQt5.QtWidgets.QMenu, name=label)
        if not label_menu:
            label_menu = plots_menu.addMenu(label)
            label_menu.setObjectName(label)

        # ignore if the graph item already exists in the menu
        if label_menu.findChild(PyQt5.QtWidgets.QAction, name=category):
            # fixme - unexpected return
            return

        def create_load_graph(label, graph_type):
            def load_graph():
                GraphManager.plot_graph(label, graph_type)

            return load_graph

        action = label_menu.addAction(category, create_load_graph(label, category))
        action.setObjectName(category)


    @staticmethod
    def save_graph(label, rmsd_data, category):
        """

        :param label: PyMOL label
        :param category: a str, e.g. 'rmsd', 'rmsf'
        :param atom_group: used for the filename information
        :param rmsd_data: rmsd results from MDAnalysis
        :return:
        """

        atom_group = MDAnalysisManager.getSystem(label)
        selection = MDAnalysisManager.getSelection(label)

        # check if the rmsd directory exists
        filepath_hash = GraphManager._get_hash_from_filepath(atom_group.universe.filename)
        selection_hash = GraphManager._get_hash_from_selection(selection)

        datafile_dir = os.path.join(GraphManager.PLOTS_DIR, filepath_hash, selection_hash, category)
        datafile_name = os.path.join(datafile_dir, 'graph.dat')
        if not os.path.isdir(datafile_dir):
            os.makedirs(datafile_dir)

        # save the data in a file
        numpy.savetxt(datafile_name, rmsd_data)

        # copy the plotting file from the templates
        template_rmsd_plotter = os.path.join(GraphManager.TEMPLATE_DIR, '%s.py' % category)
        plotter_filename = os.path.join(GraphManager.PLOTS_DIR, filepath_hash, selection_hash, category, 'graph.py')
        shutil.copyfile(template_rmsd_plotter, plotter_filename)

        # update the GUI menu to give access to the graph
        GraphManager._add_menu_item(label, category)


    @staticmethod
    def plot_graph(label, category):
        # use the basic template to visualise the results

        atom_group = MDAnalysisManager.getSystem(label)
        selection = MDAnalysisManager.getSelection(label)

        filepath_hash = GraphManager._get_hash_from_filepath(atom_group.universe.filename)
        selection_hash = GraphManager._get_hash_from_selection(selection)
        graph_dir = os.path.join(GraphManager.PLOTS_DIR, filepath_hash, selection_hash, category)

        sys.path.append(graph_dir)
        # import the saved rmsd plotter
        module_name = 'graph'
        # # check if the graph already has been imported / graphed
        if module_name in sys.modules:
            # already imported so reload
            plotter = importlib.reload(sys.modules[module_name])
        else:
            plotter = importlib.import_module(module_name)
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

