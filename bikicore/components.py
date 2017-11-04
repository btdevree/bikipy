"""Classes representing physical or conceptual objects needed to describe the
biochemical system, the model, and the experiments performed. 

"""

from traits.api import HasTraits, Str, ListStr, Int
from bikipy.bikicore.exceptions import ComponentNotValidError

# Define classes for the different components of the biochemical system
class Drug(HasTraits):
    
    # Traits initialization
    name = Str('adrenaline')
    symbol = Str('A')
    
class Protein(HasTraits):
    
    # Traits initialization
    name = Str('beta adrenergic receptor')
    symbol = Str('R')
    conformation_names = ListStr(['inactive', 'active']) 
    conformation_symbols = ListStr (['', '*'])
    
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
    IDnumber = Int(1)

# Define classes for the different types of state transitions - edges on network graph
class StateTransition(HasTraits):    
    #Superclass for all transition objects
    
    # Traits initialization
    number = Int(1)
    IDnumber = Int(1)
    
class ConformationalChange(StateTransition):
    
    # Traits initialization
    name = Str('activation')
    
class Association(StateTransition):
    
    # Traits initialization
    name = Str('association')
    
class Dissociation(StateTransition):
    
    # Traits initialization
    name = Str('dissociation')
    
    
    