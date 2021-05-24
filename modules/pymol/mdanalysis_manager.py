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
import hashlib
import json

import sqlite3
import MDAnalysis

from . import cmd
from . import CmdException

PYMOL_STORAGE=os.path.join(os.path.expanduser('~'), '.pymol')
# create a ~/.pymol directory if it does not exist
if not os.path.isdir(PYMOL_STORAGE):
    os.makedirs(PYMOL_STORAGE)

class SelectionHistoryManager():
    # remembers the created labels with the atomic indices for each file (like .bash_history)
    SELECTION_HISTORY = os.path.join(PYMOL_STORAGE, 'mda_selection_history.db')

    con = sqlite3.connect(SELECTION_HISTORY)

    # connect to the database
    cursor = con.cursor()

    # create a table for selections if it does not yet exist
    cursor.execute(
        """CREATE TABLE IF NOT EXISTS selections(
            filehash text,
            label text,
            atom_ids text,
            PRIMARY KEY (filehash, label)
            )"""
    )
    con.commit()

    @staticmethod
    def addSelection(filepath, label, atom_ids):
        """
        Adds the selection to the database. If there is already such selection for this
        filename, overwrite the previous selection.
        :param filepath:
        :param label:
        :param atom_ids:
        :return:
        """
        hash = SelectionHistoryManager._getHash(filepath)
        SelectionHistoryManager.cursor.execute(
            """INSERT OR REPLACE INTO selections(filehash, label, atom_ids)
            VALUES(?, ?, ?)
            """, (hash, label, atom_ids)
        )
        SelectionHistoryManager.con.commit()

    @staticmethod
    def getSelectionsLabels(filepath):
        """
        Retrieves the selection labels associated with a given file.
        :param filepath: universe topology filename
        :return: list of selction labels
        """
        hash = SelectionHistoryManager._getHash(filepath)
        SelectionHistoryManager.cursor.execute("SELECT label FROM selections "
                                               "WHERE filehash=?", (hash, ))
        labels = SelectionHistoryManager.cursor.fetchall()
        return [label[0] for label in labels]

    @staticmethod
    def getMdaSelection(filepath, label):
        hash = SelectionHistoryManager._getHash(filepath)
        SelectionHistoryManager.cursor.execute("SELECT atom_ids FROM selections "
                                               "WHERE filehash=? and label=?",
                                               (hash, label))
        atom_ids = SelectionHistoryManager.cursor.fetchone()
        return 'bynum ' + atom_ids[0]


    @staticmethod
    def deleteMdaSelection(filepath, label):
        hash = SelectionHistoryManager._getHash(filepath)
        SelectionHistoryManager.cursor.execute("DELETE FROM selections "
                                               "WHERE filehash=? and label=?",
                                               (hash, label))
        SelectionHistoryManager.con.commit()
        return


    @staticmethod
    def _getHash(filepath):
        return hashlib.md5(open(filepath).read().encode('utf-8')).hexdigest()


class MDAnalysisManager():
    """
    Stores all meta-data related to MDAnalysis code, like handles to the trajectories.

    TODO To Think About:
     - what if there is more trajectories which have their own callbacks?
    """

    # contain a list of loaded objects / states / names before we know how to extract them ourselves
    # This object is saved along with a session
    Systems = {}

    # A list of callbacks for rendering
    # This object is saved along with a session
    callbacks = {}
    # constants
    MDA_FRAME_CHANGED_CALLBACK = "MDA_FRAME_CHANGED_CALLBACK"

    # fixme - delete - it's been moved to graph_manager
    TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'plotting_templates')
    PLOTS_DIR = os.path.join(PYMOL_STORAGE, 'plotting')

    SESSION = None
    SESSION_PATH = None


    @staticmethod
    def toJSON():
        """

        :return: string json
        """
        # Save all the important metadata from which the existing functionality can be recovered
        # Metadata:
        #  - labels
        #  - filenames
        #  - selection string

        metadata = {'SESSION': MDAnalysisManager.SESSION,
                    'SESSION_PATH': MDAnalysisManager.SESSION_PATH}
        # list of labels where each points to metadata
        metadata['labels'] = {}
        labels = metadata['labels']
        for label, item in MDAnalysisManager.Systems.items():
            atom_group = item['system']
            selection = item['selection']

            labels[label] = {}
            labels[label]['selection'] = selection

            filenames = {}
            filenames['top_filename'] = atom_group.universe.filename
            # fixme - note that when multiple trajectories are loaded, it becomes trajectory.filenames
            filenames['trajectory'] = atom_group.universe.trajectory.filename
            labels[label]['filenames'] = filenames

        return json.dumps(metadata)



    @staticmethod
    def fromJSON(json_metadata):
        """
        Initialise the MDAnalysis labels and metadata
        :param json_metadata:
        :return:
        """
        # Reinitialise this class static methods
        metadata = json.loads(json_metadata)
        MDAnalysisManager.SESSION = metadata['SESSION']
        MDAnalysisManager.SESSION_PATH = metadata['SESSION_PATH']
        for label, meta in metadata['labels'].items():
            MDAnalysisManager.load(label, meta['filenames']['top_filename'])
            MDAnalysisManager.select(label, label, meta['selection'])

            # load the trajectory for the filename (if it's there?)
            if 'trajectory' in meta['filenames']:
                MDAnalysisManager.loadTraj(label, meta['filenames']['trajectory'])
                MDAnalysisManager.renderWithMDAnalysis(label)
            # fixme - handle the multiple trajectory filenames

        print('reloaded')



    @staticmethod
    def getSystems():
        return MDAnalysisManager.Systems


    @staticmethod
    def getSystem(label):
        return MDAnalysisManager.Systems[label]['system']


    @staticmethod
    def getSelection(label):
        return MDAnalysisManager.Systems[label]['selection']


    @staticmethod
    def select(label, new_label, selection):
        atom_group = MDAnalysisManager.Systems[label]['system'].select_atoms(selection)
        MDAnalysisManager.Systems[new_label] = {'system': atom_group, 'selection': selection}

        # Compatibility with PyMOL
        # convert to indexes for PyMOL internals
        # fixme - this should be moved to some general helper function module
        def get_consecutive_index_ranges(atom_ids):
            # This is the accepted PyMOL format: "index 2100-2200 + index 2300-2400"
            import more_itertools as mit
            grouped = [list(group) for group in mit.consecutive_groups(atom_ids)]
            return ' + '.join(['index %d-%d' % (g[0], g[-1]) for g in grouped])

        pymol_selection = get_consecutive_index_ranges(atom_group.ids)
        try:
            cmd.select(new_label, pymol_selection)
        except Exception as exception:
            print('MDAnalysis encouranted an issue with recreating the selections')
            if selection != 'all':
                raise exception

        # update the selection history
        # selections that are too large are not stored / remembered
        SelectionHistoryManager.addSelection(atom_group.universe.filename, new_label, ' '.join(atom_group.ids.astype(str)))
        # if len(atom_group) < 5000:
        #     # update the selection history file
        #     with open(MDAnalysisManager.SELECTION_HISTORY, 'a') as FHIST:
        #         # store filepath, label, atomic positions
        #         output = json.dumps({'fp' : atom_group.universe.filename,
        #                   'l' : new_label,
        #                   'a' : ' '.join(atom_group.ids.astype(str))})
        #         FHIST.write(output + '\n')
        #         pass




    @staticmethod
    def newLabel(label, atom_group, selection):
        MDAnalysisManager.Systems[label] = {'system': atom_group, 'selection': selection}


    @staticmethod
    def updateLabel(old, new):
        MDAnalysisManager.Systems[new] = MDAnalysisManager.Systems[old]
        del MDAnalysisManager.Systems[old]


    @staticmethod
    def exists(label):
        if label in MDAnalysisManager.Systems:
            return True

        return False


    @staticmethod
    def load(label, topology_filename):
        """
        Loads the topology file
        :param label: the PyMOL label used in the system, which the user can see and recognize
        """

        u = MDAnalysis.Universe(topology_filename)
        # selection is the string used for selection of the atom group
        MDAnalysisManager.Systems[label] = {'system': u.atoms, 'selection': 'all'}


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
        atom_group = MDAnalysisManager.Systems[label]['system']

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
            for atom_group_label in MDAnalysisManager.Systems.keys():
                atom_group = MDAnalysisManager.Systems[atom_group_label]['system']
                # MDAnalysis: switch to the requested frame
                atom_group.universe.trajectory[int(frame) - 1]

                # 0 is for append, 1 we guess means replacing
                replace = 1
                try:
                    cmd.load_coordset(atom_group.positions, atom_group_label, replace)
                except CmdException as exception:
                    # check if this is the LoadCoords-Error. PyMOL cannot apply load_coordset to a selection.
                    # A full object is needed (ie "Extract to object" feature)
                    # unfortunately the exception here cannot be recognized, because it contains no specific information
                    # which is why I squash all errors
                    # fixme - switch to using always objects instead of selections?
                    print('Swallowed an exception - hopefully only LoadCoords-Error exception.')
                    pass

        # MDAnalysis universe
        atom_group = MDAnalysisManager.Systems[label]['system']

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