"""
Python functions to load and display refinement data from Jana2006

Version 0.1 12/8/19

Dan Porter
August 2019
"""

import os
from .jana_functions import readm40
from .jana_functions import refine
from .jana_functions import refinementtable


class Refine:
    def __init__(self, filename):
        (dirName, filetitle) = os.path.split(filename)
        (fname, Ext) = os.path.splitext(filetitle)

        self.filename = filename
        self.title = fname

        crys, error = readm40(os.path.join(dirName, fname + '.m40'))
        self.crystal_dict = crys
        self.error_dict = error

    def update(self, notes=''):
        refine(self.filename, notes)

    def create_table(self):
        refinementtable(self.filename)
