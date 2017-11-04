"""Objects and logic for TraitsUI based GUI to create biochemical models. 

"""

from traits.api import HasTraits, Button
import bikipy.bikicore.components as bkcc

class ModelCreator(HasTraits):
    """Base object that contains the components needed to build up or edit a 
    model
    """  
    newbutton = Button("Create New Model")
    
# Create the ModelCreator object:
mc = ModelCreator()

# Start up ModelCreator object if script is called by python intrepreter
if __name__ == "__main__":
    mc.configure_traits()