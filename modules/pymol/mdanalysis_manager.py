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


class MDAnalysisManager():
    """
    Stores all data related to MDAnalysis code. This most often means
    storing the handles to the trajectories.
    """

    # contain a list of loaded objects / states / names before we know how to extract them ourselves
    MDAnalysisSystems = {}

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

