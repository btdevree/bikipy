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
        new_model = bkcm.create_new_model('new', self.model_list)
        self.model_list.extend([self.new_model])
        self.select_model(new_model)
    
    def copy_model(self):
        parent_model = UserSelectModel(self.model_list)
        new_model = bkcm.create_new_model('copy', self.model_list)
        self.model_list.extend([self.new_model])
        self.select_model(new_model)
    
    def select_model(self):
        self.current_model = new_model

    
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
        
    #Setup GUI window options
    traits_view = config_model_view

# Start up ModelCreator object if script is called by python intrepreter
if __name__ == "__main__":
          
    # Create the ModelCreator object and start the traits GUI
    mc = ModelCreator()    
    mc.configure_traits()