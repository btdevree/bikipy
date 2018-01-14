"""Objects and logic for TraitsUI based GUI to create biochemical models. 

"""

from traits.api import HasTraits, Instance, List
from traitsui.api import View, Menu, MenuBar, Action, Handler, NoButtons, Group
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
        self.model_list.append(new_model)
        self.current_model = new_model
    
    def copy_model(self):
        parent_model = UserSelectModel(self.model_list)
        new_model = bkcm.create_new_model('copy', self.model_list, parent_model)
        self.model_list.append(new_model)
        self.current_model = new_model
    
    def select_model(self, selected_model):
        selected_model = UserSelectModel(self.model_list)
        self.current_model = selected_model
   
    # Method to configure the view
    def config_model_view(self):
        menubar_config = MenuBar(
                            Menu(Action(name='New Model',
                                        action='new_model',
                                        tooltip='Create a new model from scratch'), 
                                 Action(name='Copy Model',
                                        action='copy_model',
                                        tooltip='Create a new model by copying an existing one'),
                                            name = 'Create Model'))
        
        view_config = View(
            Group(
                Item(name = 'model_list', style = 'custom'),
                show_border = True),
            Item(label="Lots of stuff should go here"),
            menubar = menubar_config,
            buttons=NoButtons,
            title = 'BiKiPy Modeler'
            )
        
        return(view_config)
        
    #Setup GUI window options
    traits_view = config_model_view

# Start up ModelCreator object if script is called by python intrepreter
if __name__ == "__main__":
          
    # Create the ModelCreator object and start the traits GUI
    mc = ModelCreator()    
    mc.configure_traits()