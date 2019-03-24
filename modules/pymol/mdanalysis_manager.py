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

"""

import MDAnalysis
from enum import Enum
from . import cmd


class MEMORY_MODE(Enum):
    """
    Decides whether PYMOL should load/store the trajectory data.
    MDAnalysis does not load the trajectory into memory and fetches the frames it needs at any time.

    PYMOL mode is added for testing purposes. However, in the long term MDAnalysis should take over this function.

    If PYMOL is used, MDAnalysis will still load the metadata, but it won't be used for rendering the molecules.
    """
    MDANALYSIS = 1  # MDAnalysis takes over loading of trajectories. PyMOL stores only metadata.
    PYMOL = 2       # PyMOL reads the trajectories into RAM using its own settings.

class MDAnalysisManager():
    """
    Stores all data related to MDAnalysis code. This most often means
    storing the handles to the trajectories.

    TODO To Think About:
     - what if there is more trajectories which have their own callbacks?
    """

    # contain a list of loaded objects / states / names before we know how to extract them ourselves
    MDAnalysisSystems = {}
    MODE = MEMORY_MODE.MDANALYSIS

    # A list of callbacks that are related to MDAnalysis
    callbacks = {}

    # constants
    MDA_FRAME_CHANGED_CALLBACK = "MDA_FRAME_CHANGED_CALLBACK"


    @staticmethod
    def getMDAnalysisSystems():
        return MDAnalysisManager.MDAnalysisSystems


    @staticmethod
    def load(label, topology_filename):
        """
        Loads the topology file
        :param label: the PyMOL label used in the system, which the user can see and recognize
        """

        u = MDAnalysis.Universe(topology_filename)
        MDAnalysisManager.MDAnalysisSystems[label] = u


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
        u = MDAnalysisManager.MDAnalysisSystems[label]

        # load the topology with its trajectory
        u.load_new(trajectory_filename)

        # set up frame slider (PyMOL movie panel)
        cmd.mset('1x{}'.format(u.trajectory.n_frames))

        if MDAnalysisManager.MODE == MEMORY_MODE.MDANALYSIS:
            MDAnalysisManager.renderWithMDAnalysis(label)


    @staticmethod
    def renderWithMDAnalysis(label):

        def fetch_frame_coordinates(frame, label):
            '''
            Updates coordinates in state 1 from universe frame
            :param frame: 1-based frame index
            '''
            universe = MDAnalysisManager.MDAnalysisSystems[label]
            cmd.load_coordset(universe.trajectory[int(frame) - 1].positions, label, 1)

        # MDAnalysis universe
        universe = MDAnalysisManager.MDAnalysisSystems[label]

        # This should be the default
        cmd.set('retain_order', 1, label)

        # reduce PyMOL's logging (optional)
        cmd.feedback('disable', 'executive', 'actions')

        # set the per-frame call to update coordinates in state 1 ("in place")
        MDAnalysisManager.callbacks[MDAnalysisManager.MDA_FRAME_CHANGED_CALLBACK] = [fetch_frame_coordinates, label]
        for frame in range(1, universe.trajectory.n_frames + 1):
            cmd.mdo(frame, '{}.{}'.format(MDAnalysisManager.MDA_FRAME_CHANGED_CALLBACK, frame))