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
 - FIXME: Thomas: when trajectory is not loaded from a trajectory file, but e.g. with a
    script which adds frames, or by loading a set of PDB files into the
    same object, or by copying an object within PyMOL.
 - TODO: PyMOL and MDAnalysis: create an issue on mdanalysis git hub and ask if you can add frames
    to the universe. This would come in handy in PyMOL.
 - TODO: MDAnalysis: "trajectory.filenames" should always be a list that contains the
    loaded trajectories. Let's create an issue and ensure a consistent behaviour.
 - TODO: We're in need of IDs for each selection so that the label change would not affect the simulation

"""

import os
import MDAnalysis
from enum import Enum
from . import cmd
import json


class MDAnalysisManager():
    """
    Stores all meta-data related to MDAnalysis code, like handles to the trajectories.

    TODO To Think About:
     - what if there is more trajectories which have their own callbacks?
    """

    # contain a list of loaded objects / states / names before we know how to extract them ourselves
    # This object is saved along with a session
    MDAnalysisSystems = {}

    # A list of callbacks for rendering
    # This object is saved along with a session
    callbacks = {}
    # constants
    MDA_FRAME_CHANGED_CALLBACK = "MDA_FRAME_CHANGED_CALLBACK"

    # fixme - should in the configuration file
    TEMPLATE_DIR = os.path.join(__file__.rsplit('/', maxsplit=1)[0], 'plotting_templates')
    PLOTS_DIR = os.path.expanduser('~/.pymol/plotting/plots')


    @staticmethod
    def toJSON():
        # Save all the important metadata from which the existing functionality can be recovered
        # Metadata:
        #  - labels
        #  - filenames
        #  - selection string

        metadata = {}
        for label, atom_group in MDAnalysisManager.MDAnalysisSystems.items():
            metadata['label'] = label
            metadata['selection'] = atom_group._pymol_used_selection

            filenames = {}
            filenames['top_filename'] = atom_group.universe.filename
            filenames['trajectory'] = atom_group.universe.trajectory.filename
            metadata['filenames'] = filenames


        return json.dumps(metadata)



    @staticmethod
    def fromJSON(json):
        # Reinitialise this class static methods
        pass


    @staticmethod
    def getMDAnalysisSystems():
        return MDAnalysisManager.MDAnalysisSystems


    @staticmethod
    def updateLabel(old, new):
        MDAnalysisManager.MDAnalysisSystems[new] = MDAnalysisManager.MDAnalysisSystems[old]
        del MDAnalysisManager.MDAnalysisSystems[old]


    @staticmethod
    def load(label, topology_filename):
        """
        Loads the topology file
        :param label: the PyMOL label used in the system, which the user can see and recognize
        """

        u = MDAnalysis.Universe(topology_filename)
        MDAnalysisManager.MDAnalysisSystems[label] = u.atoms

        # fixme - a hacky way to store selections for each atom group
        u.atoms._pymol_used_selection = 'all'


    @staticmethod
    def loadTraj(label, trajectory_filename):
        """
        Load the trajectory universe into the existing label.
        fixme: How would this work if the trajectory was the topology? ie the .pdb file.
        :param label: The name of the
        :param trajectory_filename:
        :return:
        """

        # get the universe for the label
        atom_group = MDAnalysisManager.MDAnalysisSystems[label]

        # load the topology with its trajectory
        atom_group.universe.load_new(trajectory_filename)

        # set up frame slider (PyMOL movie panel)
        # fixme - what if there are two separate simulations? separate sliders? focus?
        cmd.mset('1x{}'.format(atom_group.universe.trajectory.n_frames))

        MDAnalysisManager.renderWithMDAnalysis(label)


    @staticmethod
    def renderWithMDAnalysis(label):

        def fetch_frame_coordinates(frame):
            '''
            Updates coordinates to the selected frame in each universe.
            fixme - should update them only for the selected universe?
            :param frame: 1-based frame index
            '''
            for atom_group_label in MDAnalysisManager.MDAnalysisSystems.keys():
                atom_group = MDAnalysisManager.MDAnalysisSystems[atom_group_label]
                # MDAnalysis: switch to the requested frame
                atom_group.universe.trajectory[int(frame) - 1]

                # 0 is for append, 1 we guess means replacing
                replace = 1
                cmd.load_coordset(atom_group.positions, atom_group_label, replace)

        # MDAnalysis universe
        atom_group = MDAnalysisManager.MDAnalysisSystems[label]

        # This should be the default
        cmd.set('retain_order', 1, label)

        # reduce PyMOL's logging (optional)
        cmd.feedback('disable', 'executive', 'actions')

        # set the per-frame call to update coordinates in state 1 ("in place")
        MDAnalysisManager.callbacks[MDAnalysisManager.MDA_FRAME_CHANGED_CALLBACK] = fetch_frame_coordinates
        for frame in range(1, atom_group.universe.trajectory.n_frames + 1):
            cmd.mdo(frame, '{}.{}'.format(MDAnalysisManager.MDA_FRAME_CHANGED_CALLBACK, frame))


# fixme - there is a better place to initialise this, during the installation?
if not os.path.isdir(MDAnalysisManager.PLOTS_DIR):
    os.makedirs(MDAnalysisManager.PLOTS_DIR)