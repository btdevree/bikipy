"""Class for the model objects used throughout the program.
"""

import uuid
import bikipy.bikicore.components as bkcc
from traits.api import HasTraits, Int, Str, Instance, This, List

#Model class
class Model(HasTraits):
    
    #Initalize traits
    number = Int
    name = Str
    parent_model = Instance(This)
    ID = Instance(uuid.UUID)
    drug_list = List(bkcc.Drug)
    protein_list = List(bkcc.Protein)
    #compartment_list = List(bkcc.Compartment) #To be implemented in future
    rule_list = List(bkcc.Rule)
    
    def __init__(self, number, name, parent_model, *args, **kwargs):
        super().__init__(*args, **kwargs) #Make sure to call the HasTraits initialization machinery
        self.number = number
        self.name = name
        self.parent_model
        self.ID = uuid.uuid4()
    
    def _copy_from(model_to_copy):
        #Code to duplicate the structure of a model with new objects
        pass           
          
#Model creation method
def create_new_model(new_model_type, model_list, model_to_copy = None):
    new_number = _find_next_model_number(model_list)
    
    if new_model_type == 'new':
        new_model = Model(new_number, 'New Model', None)
    else:
        new_name = [model_to_copy.name, '-copy']
        if new_model_type == 'copy':
            new_model = Model(new_number, new_name, None)
        elif new_model_type == 'new_child':
            new_model = Model(new_number, new_name, model_to_copy)
        elif new_model_type == 'copy_child':
            new_model = Model(new_number, new_name, model_to_copy.parent_model)
        new_model._copy_from(model_to_copy)
    return new_model
    
def _find_next_model_number(model_list):
    found_numbers = {model.number for model in model_list}
    current_int = 1
    while ({current_int} & found_numbers) != set():
       current_int += 1
    return current_int
    
        
    