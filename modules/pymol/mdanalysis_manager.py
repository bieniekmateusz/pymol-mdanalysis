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

        # Load the MDAnalysis part
        mdsystems = MDAnalysisManager.getMDAnalysisSystems()

        # if the object name is the same as the previously loaded one,
        # then reload MDAnalysis together with the trajectory
        if mdsystems[metadata['oname']] != None:
            current = mdsystems[metadata['oname']]
            current['trajectory'] = metadata['fname']
            current['universe'] = MDAnalysis.Universe(current['topology'], current['trajectory'])

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

        # Load the MDAnalysis part
        mdsystems = MDAnalysisManager.getMDAnalysisSystems()

        # get the name of the object, and the file, and use it to set up the MDAnalysis object
        # contain a list of loaded objects / states / names before we know how to extract them ourselves
        mdsystems[metadata['oname']] = {
            'mdanalysis_universe': MDAnalysis.Universe(metadata['finfo']),
            'topology': metadata['finfo']
        }

        return reply
    return decorator


class MDAnalysisManager():
    """
    Stores all data related to MDAnalysis code. This most often means
    storing the handles to the trajectories.
    """

    MDAnalysisSystems = {}

    @staticmethod
    def getMDAnalysisSystems():
        return MDAnalysisManager.MDAnalysisSystems

