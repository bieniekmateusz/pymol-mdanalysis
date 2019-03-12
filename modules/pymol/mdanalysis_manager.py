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


import MDAnalysis

"""
A decorator for the PyMOL "load_traj" function. It loads the corresponding 
MDAnalysis universe for the loaded trajectory in PyMOL. 
"""
def mdaLoadTraj(pymolLoadTraj):
    def decorator(*args, **kwargs):
        # args and kwargs are what was called

        # I suspect that the reply holds information on any error
        # that would occur during loading the trajectory
        reply, metadata = pymolLoadTraj(*args, **kwargs)
        label = metadata['label']
        trajectory_filename = metadata['trajectory_filename']

        # Load the MDAnalysis part
        MDAnalysisManager.loadTraj(label, trajectory_filename)

        return reply
    return decorator


"""
A decorator for the PyMOL load() function. It loads the corresponding 
MDAnalysis universe for the loaded file in PyMOL.
"""
def mdaLoad(pymolLoad):
    def decorator(*args, **kwargs):
        # args and kwargs are what was called

        # I suspect that the reply holds information on any error
        # that would occur during loading the trajectory
        reply, metadata = pymolLoad(*args, **kwargs)
        label = metadata['oname']
        topology_filename = metadata['finfo']

        # Load into MDAnalysis
        MDAnalysisManager.load(label, topology_filename)
        return reply
    return decorator


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
    def loadTraj(label, trajectory):
        """
        Load the trajectory universe into the existing label.
        fixme: How would this work if the trajectory was the topology? ie the .pdb file.
        :param label: The name of the
        :param trajectory:
        :return:
        """

        # get the universe for the label
        u = MDAnalysisManager.MDAnalysisSystems[label]
        topology_file = u.filename

        # load the topology with its trajectory
        u_withTraj = MDAnalysis.Universe(topology_file, trajectory)
        MDAnalysisManager.MDAnalysisSystems[label] = u_withTraj

