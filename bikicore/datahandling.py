"""Class for the solver object used throughout the program.
"""

import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import bikipy.bikicore.model as bkcm
import bikipy.bikicore.components as bkcc
import bikipy.bikicore.datahandling as bkcd
from traits.api import HasTraits, Int, Str, Instance, This, List


# Solver class
class Solver(HasTraits):
    
    # Initalize traits
    model = Instance(bkcm.Model)
    experiment_list = List(Instance(bkcd.Experiment))
    
    def __init__(self, model, *args, **kwargs):
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery
        self.model = model
        self.ID = uuid.uuid4()
        self.experiment_list = []
        
    