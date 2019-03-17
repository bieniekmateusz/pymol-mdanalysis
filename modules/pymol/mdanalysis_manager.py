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


class MDACallbacks:
    """
    A list of callbacks that are related to MDAnalysis
    """

    callbacks = {}


class MDAnalysisManager():
    """
    Stores all data related to MDAnalysis code. This most often means
    storing the handles to the trajectories.

    TODO To Think About:
     - what if there is more trajectories which have their own callbacks?
    """

    # contain a list of loaded objects / states / names before we know how to extract them ourselves
    MDAnalysisSystems = {}



    class MEMORY_MODE(Enum):
        PYMOL = 1   # PyMOL reads the trajectories into RAM using its own settings
        MDANALYSIS = 2  # MDAnalysis takes over loading of trajectories. PyMOL stores only metadata.
        PYMOL_MDANALYSIS = 3    # 1 and 2 together

    MODE = MEMORY_MODE.MDANALYSIS

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
        topology_file = u.filename

        # load the topology with its trajectory
        u_withTraj = MDAnalysis.Universe(topology_file, trajectory_filename)
        MDAnalysisManager.MDAnalysisSystems[label] = u_withTraj

        # FIXME do i need to know whether it's mdanalysis only?
        MDAnalysisManager.loadIntoPyMOL(label)


    @staticmethod
    def loadIntoPyMOL(label):

        universe = MDAnalysisManager.MDAnalysisSystems[label]
        def get_mdanalysis_load_frame(label):
            def mdanalysis_load_frame(frame):
                '''
                Updates coordinates in state 1 from universe frame

                :param frame: 1-based frame index
                '''
                universe = MDAnalysisManager.MDAnalysisSystems[label]
                cmd.load_coordset(universe.trajectory[int(frame) - 1].positions, label, 1)
            return mdanalysis_load_frame

        # cmd.load(gro)

        # This should be the default
        cmd.set('retain_order', 1, label)

        # load trajectory into MDAnalysis universe
        # universe = MDAnalysis.Universe(gro, xtc)

        # set up frame slider (PyMOL movie panel)
        cmd.mset('1x{}'.format(universe.trajectory.n_frames))

        # per-frame callback to update coordinates in state 1
        MDACallbacks.callbacks['MDA_FRAME_CHANGED_CALLBACK'] = get_mdanalysis_load_frame(label)
        for frame in range(1, universe.trajectory.n_frames + 1):
            cmd.mdo(frame, 'MDA_FRAME_CHANGED_CALLBACK.{}'.format(frame))

        # reduce PyMOL's logging (optional)
        cmd.feedback('disable', 'executive', 'actions')

        # init # TODO prettify / explain
        get_mdanalysis_load_frame(label)(1)
        cmd.zoom(animate=1)
        # cmd.mplay()


