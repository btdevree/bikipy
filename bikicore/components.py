"""Classes representing physical or conceptual objects needed to describe the
biochemical system, the model, and the experiments performed. 

"""

import uuid
from traits.api import HasTraits, Str, List, Int, Instance
from bikipy.bikicore.exceptions import ComponentNotValidError

# Define classes for the different components of the biochemical system
class Drug(HasTraits):
    
    # Traits initialization
    name = Str('adrenaline')
    symbol = Str('A')
    ID = Instance(uuid.UUID) 
    
class Protein(HasTraits):
    
    # Traits initialization
    name = Str('beta adrenergic receptor')
    symbol = Str('R')
    conformation_names = List(Str(['inactive', 'active'])) 
    conformation_symbols = List(Str (['', '*']))
    ID = Instance(uuid.UUID) 
    
    # Check if the trait values are valid. Use after user editing of the object, for example.
    def check_protein_traits(self):
        
        # Check if the number of conformation names and symbols are the same
        if len(self.conformation_names) != len(self.conformation_symbols):
            raise ComponentNotValidError('The number of protein conformation names and symbols must be the same')

# Define class for the different states - nodes on network graph
class State(HasTraits):
    
    # Traits initialization
    name = Str('1')
    number = Int(1)
    ID = Instance(uuid.UUID)

# Define classes for the different types of state transitions - edges on network graph
class StateTransition(HasTraits):    
    #Superclass for all transition objects
    
    # Traits initialization
    number = Int(1)
    ID = Instance(uuid.UUID)
    
    
    
class ConformationalChange(StateTransition):
    
    # Traits initialization
    name = Str('activation')
    
class Association(StateTransition):
    
    # Traits initialization
    name = Str('association')
    
class Dissociation(StateTransition):
    
    # Traits initialization
    name = Str('dissociation')
    
    
    