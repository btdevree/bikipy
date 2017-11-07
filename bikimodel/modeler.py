"""Objects and logic for TraitsUI based GUI to create biochemical models. 

"""

from traits.api import HasTraits, Instance, List
from traitsui.api import View, Menu, MenuBar, Action, Handler, NoButtons
from traitsui.item import Item
import bikipy.bikicore.components as bkcc
import bikipy.bikicore.model as bkcm

class ModelCreator(HasTraits):
    """Base GUI object that contains the components needed to build up or edit a 
    model
    """  
    
    #Initialize traits
    current_model = Instance(bkcm.Model)
    model_list = List(Instance(bkcm.Model))    
        
    #Menu Bar methods
    def new_model(self):
        print(self.current_model)
        print(bkcm.create_new_model(10))
        new_model = bkcm.create_new_model(3)
        self.current_model = new_model
        print(new_model)
        print(self.current_model)
        self.model_list = self.model_list.append(self.current_model)
    
    new_menu_action = Action(
        name='New Model',
        action='new_model',
        tooltip='Create a new model')      
    
    # Method to configure the view
    def config_model_view(self):
        menubar_config = MenuBar(
                            Menu(self.new_menu_action, name = 'New Model'))
        
        view_config = View(
            Item(label="Lots of stuff should go here"),
            menubar = menubar_config,
            buttons=NoButtons
            )
        
        return(view_config)
        

    traits_view = config_model_view
    
    

# Create the ModelCreator object:
mc = ModelCreator()

# Start up ModelCreator object if script is called by python intrepreter
if __name__ == "__main__":
    mc.configure_traits()