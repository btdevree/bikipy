"""Class for the model objects used throughout the program.
"""

import uuid
import bikipy.bikicore.components as bkcc
from traits.api import HasTraits, Int, Str, Instance, This, List

# Model class
class Model(HasTraits):
    
    # Initalize traits
    number = Int
    name = Str
    parent_model = Instance(This)
    ID = Instance(uuid.UUID)
    drug_list = List(Instance(bkcc.Drug))
    protein_list = List(Instance(bkcc.Protein))
    # compartment_list = List(bkcc.Compartment) #To be implemented in future
    rule_list = List(Instance(bkcc.Rule))
    
    def __init__(self, number, name, parent_model, *args, **kwargs):
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery
        self.number = number
        self.name = name
        self.parent_model
        self.ID = uuid.uuid4()
        self.drug_list = []
        self.protein_list = []
        self.rule_list = []
    
    def _copy_from(self, model_to_copy):
        # Code to duplicate the structure of a model with new objects
        pass           
    
    def generate_network(self):
        # Create a network graph of the model by using the list of rules. 
        
        # Create a new Network object and work with main_graph for the next block
        self.network = bkcc.Network()
      
        # Add singleton states to graph
        for current_component in self.drug_list:
            new_state = bkcc.State()
            new_state.required_drug_list = [current_component]
            new_state.required_protein_list = []
            new_state.req_protein_conf_lists = [[]]
            self.network.main_graph.add_node(new_state)
        for current_component in self.protein_list:
            new_state = bkcc.State()
            new_state.required_drug_list = []
            new_state.required_protein_list = [current_component]
            new_state.req_protein_conf_lists = [list(range(len(current_component.conformation_names)))] # Allow all conformations
            self.network.main_graph.add_node(new_state)
        
        # Apply each rule
        for current_rule in self.rule_list:
            
            # Break into basic relationships
            if current_rule.rule == ' associates with ':
                
                #fill in states directly for now
                states = [x for x in self.network.main_graph.__iter__()]
                state1 = states[0]
                state2 = states[1]
                self._create_association(state1, state2)
            # Make association steps
        
    def _create_association(self, state1, state2):
        # Helper function to connect two states into an association relationship
        # Creates a new associated state if one cannot be found already
        
        # Look for an associated state
        associated_drug_list = [*state1.required_drug_list, *state2.required_drug_list]
        associated_protein_list = [*state1.required_protein_list, *state2.required_protein_list]
        associated_conf_list = [*state1.req_protein_conf_lists, *state2.req_protein_conf_lists]
        #DO LATER, search for a state with these components, assign to state12
        
        # Make a new associated state
        state12 = bkcc.State()
        state12.required_drug_list = associated_drug_list
        state12.required_protein_list = associated_protein_list
        state12.req_protein_conf_lists = associated_conf_list
        
        # Connect states 
        self.network.main_graph.add_edge(state1, state12)# (NeworkX silently adds a new node if state12 does not already exist)
        self.network.main_graph.add_edge(state2, state12)
        
        #new_states = [bkcc.State(x) for x in [self.drug_list, self.protein_list]]


    
# Model creation method
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

# Model creation helper method
def _find_next_model_number(model_list):
    found_numbers = {model.number for model in model_list}
    current_int = 1
    while ({current_int} & found_numbers) != set():
       current_int += 1
    return current_int


        
    