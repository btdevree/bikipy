"""Classes representing physical or conceptual objects needed to describe the
biochemical system, the model, and the experiments performed. 

"""

from traits.api import HasTraits, Str, ListStr

# ------ Classes ------
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
