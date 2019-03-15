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
        # Create a new network graph of the model by using the list of rules. 
        
        # Create a new Network object
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
        
        # After creating the singleton graph, build the network with the given rules
        self.apply_rules_to_network()
        
    def apply_rules_to_network(self, graph = None):
        # Apply the model's rules to an existing graph

        # If default, work on the main graph
        if graph == None:
            graph = self.network.main_graph
        
        # Apply each rule to the graph
        for current_rule in self.rule_list:
            
            # Find states that fit the rule discription
            matching_subject_states, matching_object_states = self._find_states_that_match_rule(current_rule)

            # Break into basic relationships
            if current_rule.rule == ' associates with ':
                
                # We need to count the number of components that the rule implies, or else we create unintended oligomerization
                minimum_components = [*current_rule.rule_subject, *current_rule.rule_object]
                component_count_dict = self._count_components(minimum_components)
                
                self._create_association(graph, component_count_dict, matching_subject_states, matching_object_states)
                
            elif current_rule.rule == ' dissociates_from ':
                pass
    
    def _find_states_that_match_rule(self, rule):
        # Helper function that looks through a graph and returns lists of subject 
        # and object states that include a rule's required components
            
        # Look through each component that is needed for the rule
        def state_and_rule_match(current_state, rule, what_to_search): 
            # Nested function allows reuse for subjects and objects
            
            # Pick the right components for rule subjects or objects
            if what_to_search == 'subject':
                rule_components = rule.rule_subject
                rule_conformations = rule.subject_conf
            elif what_to_search == 'object':
                rule_components = rule.rule_object
                rule_conformations = rule.object_conf
            
            # We will create list of list of booleans that ask if a component in the current state under consiteration fufills contains the 
            current_tests = []
            for current_rule_component, current_rule_conf in zip(rule_components, rule_conformations):
            
                # Check if any drugs in the current state match
                if isinstance(current_rule_component, bkcc.Drug):
                    for current_drug in current_state.required_drug_list:
                        if current_drug == current_rule_component:
                            found_drug = True
                            break
                    else:
                        found_drug = False # No matching states
                else:
                    found_drug = False # Not a drug
                
                # Check if any proteins in the current state match
                if isinstance(current_rule_component, bkcc.Protein):
                    for current_protein, current_conf in zip(current_state.required_protein_list, current_state.req_protein_conf_lists):
                        if current_protein == current_rule_component:
                            if current_conf == current_rule_conf: # Exact matching conformation
                                found_protein = True
                                break
                            elif current_rule_conf == []: # Rule allows any conformation configuration
                                found_protein = True
                                break
                    else:
                        found_protein = False # No matching states
                else:
                    found_protein = False # Not a protein
                    
                # If either type of component works for the rule requirement, record a match
                current_tests.append(any([found_drug, found_protein]))
            
            # If the state matches all the required components, return True
            return all(current_tests)
          
        # Iterate through all the graph's states and add them to the lists if they match
        matching_subject_states = []
        matching_object_states = []
        for current_state in self.network.main_graph.__iter__():
            if state_and_rule_match(current_state, rule, 'subject'):
                matching_subject_states.append(current_state)
            if state_and_rule_match(current_state, rule, 'object'):
                matching_object_states.append(current_state)
            
        # Give the lists back to the calling method
        return matching_subject_states, matching_object_states
   
    def _state_match_to_component_lists(self, test_state, component_list, conformation_list, match = 'exact'):
        # Helper function that asks if a given state contains the components given in the lists
        # Returns True or False
        
        # Use modes:
        #   match = 'exact': returns true if and only if the state contains all the required components, and no extra components
        #   match = 'minimal': returns true if the state contains all the reqired components, but may have additional ones as well
        
        # Create a copy of the component and conformation lists to consume
        remaining_components = component_list
        remaining_conformations = conformation_list
        
        # Get the components and conformations that the state requires
        state_components, state_conformations = test_state.generate_requirement_lists()
                
        # We will create list of list of booleans that ask if a component in the current state under consiteration fufills contains the 
        test_list = []
        for current_rule_component, current_rule_conf in zip(state_components, state_conformations):
        
             for current_protein, current_conf in zip(current_state.required_protein_list, current_state.req_protein_conf_lists):
                    if current_protein == current_rule_component:
                        if current_conf == current_rule_conf: # Exact matching conformation
                            found_protein = True
                            break
                        elif current_rule_conf == []: # Rule allows any conformation configuration
                            found_protein = True
                            break
                else:
                    found_drug = False # No matching states
            else:
                found_drug = False # Not a drug
            
            # Check if any proteins in the current state match
            if isinstance(current_rule_component, bkcc.Protein):
            
                else:
                    found_protein = False # No matching states
            else:
                found_protein = False # Not a protein
                
            # If either type of component works for the rule requirement, record a match
            current_tests.append(any([found_drug, found_protein]))
        
        # If the state matches all the required components, return True
        return all(current_tests)
    
    def _count_components(self, components_list):
        # Counts the number of each component in the given list
        # Returns a dictionary with components as keys and a count as an integer value
        
        # Create and empty dictionary and loop through each component in the list
        count_dict = {}
        for component in components_list:
            
            # Try to get the value with the component as key
            try: 
                number_found = count_dict[component]
            # If the component is not alreay a key, add it and set the count to 1
            except KeyError:
                count_dict[component] = 1
            # If the component exists as a key, increment the count by one
            else:
                count_dict[component] = number_found + 1
                
        # Return the counting dictionary
        return count_dict

    def _create_association(self, graph, rule_count_dict, subject_states, object_states):
        # Function to connect states into an association relationships on the given graph
        # Creates a new associated state if one cannot be found to connect
        
        # Loop through each subject and object state
        for current_sub_state in subject_states:
            for current_obj_state in object_states:
                       
                # Look for an already associated state
                associated_drug_list = [*current_sub_state.required_drug_list, *current_obj_state.required_drug_list]
                associated_protein_list = [*current_sub_state.required_protein_list, *current_obj_state.required_protein_list]
                associated_protein_conf_list = [*current_sub_state.req_protein_conf_lists, *current_obj_state.req_protein_conf_lists]
                #######################################################3
                # Make a new associated state
                state12 = bkcc.State()
                state12.required_drug_list = associated_drug_list
                state12.required_protein_list = associated_protein_list
                state12.req_protein_conf_lists = associated_protein_conf_list
                
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


        
    