"""Class for the model objects used throughout the program.
"""

import uuid
import networkx as nx
import bikipy.bikicore.components as bkcc
from collections import Counter
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

            # Each type of rule needs a different treatment
            
            # Irreversable association
            if current_rule.rule == ' associates with ':
                
                # Find states that fit the rule discription
                matching_subject_states = self._find_states_that_match_rule(current_rule, 'subject')
                matching_object_states = self._find_states_that_match_rule(current_rule, 'object')
                
                # We need to count the number of components that the rule implies, or else we create unintended oligomerization
                subject_components, subject_conformations = current_rule.generate_component_list('subject')
                object_components, object_conformations = current_rule.generate_component_list('object')
                minimum_components, minimum_conformations = current_rule.generate_component_list('both')
                subject_counter = Counter(subject_components)
                object_counter = Counter(object_components)
                associated_counter = Counter(minimum_components)
                subject_difference = associated_counter - subject_counter
                object_difference = associated_counter - object_counter
                
                # Associate any valid associated states between the matching subject and object states by calling the basic association helper method
                self._create_association(graph, component_counter, matching_subject_states, matching_object_states)
                
            elif current_rule.rule == ' dissociates_from ':
                pass
    
    def _find_states_that_match_rule(self, rule, what_to_find):
        # Helper function that looks through a graph and returns lists states that include a rule's required components
        # what_to_find = 'subject', 'object', or 'both'
        
        # Get the components and conformations that we are looking for
        # Note that any proteins with a [] conformation will always be last in the list
        rule_components, rule_conformations = rule.generate_component_list(what_to_find)
        
        # Iterate through all the graph's states and add them to the lists if they match
        matching_states = []
        for current_state in self.network.main_graph.__iter__():
            if self._state_match_to_component_lists(current_state, rule_components, rule_conformations, 'minimal'):
                matching_states.append(current_state)

        # Give the list back to the calling method
        return matching_states
   
    def _state_match_to_component_lists(self, query_state, reference_component_list, reference_conformation_list, match = 'exact'):
        # Helper function that asks if a given query state contains the components given in the reference lists
        # Returns True or False
        # Use modes:
        #   match = 'exact': returns true if and only if the query state contains all the reference components, and no extra components
        #   match = 'minimal': returns true if the query state contains all the reference components, but may have additional ones as well
        
        # Get the components and conformations that the state requires
        remaining_components, remaining_conformations = query_state.generate_component_list()

        # Loop through all remaining reference components for each query component
        for current_ref_comp, current_ref_conf in zip(reference_component_list, reference_conformation_list): 
            for i, (current_query_comp, current_query_conf) in enumerate(zip(remaining_components, remaining_conformations)):
                
                # Ask if the query and reference match
                # If there are exact conformation matches, they will be consumed first before [] conformations
                if current_query_comp == current_ref_comp and (current_query_conf == current_ref_conf or current_ref_conf == []):
                    
                    #Consume the reference match
                    del remaining_components[i]
                    del remaining_conformations[i]
                    break 
            
            else: #no break
                # If we arrive here, there was not a matching query component/conformation for a reference component/conformation
                return False # Short-circut answer
            
        # If we arrive here, matches for all query and references components/conformations were found
        # If mode is 'minimal', this means there's a match
        if match == 'minimal':
            return True
        # Otherwise, we need to make sure we completely consumed the reference lists
        elif match == 'exact':
            if remaining_components == [] and remaining_conformations == []:
                return True
            else:
                return False
        # No other match modes supported, raise error
        else:
            raise ValueError('Function _state_match_to_component_lists received an incorrect argument for parameter "match"')
    
    def _signature_match_association(self, query_signature, reference_signature):
        # Helper function that asks if a given query signature matches (any of) the reference signatures in an association context 
        # Returns True or False
        # "reference_signature" can be a single signature or a list of signatures (from "conformation included" signatures computed from rule)
        
        # Make sure we have a list of references
        if not isinstance(reference_signature, list):
            reference_signature = [reference_signature]
        
        # This allows us to see the pattern of component mixing between the states
        def calc_sig_diff(query, ref):
            pass
        
        # Work through each reference signature until we find one that works
        for refsig in reference_signature:
            
            # Need to calculate the differences between third and subject/object states
            subject_dif_ref = 
        
    
    
    
    def _count_components(self, components_list):
        # Counts the number of each component in the given list
        # Returns a dictionary with components as keys and a count as an integer value
        
        # NOTE: I think this can also be done with Python's collections.counter.
        
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

    def _create_association(self, graph, associated_counter, subject_states, object_states):
        # Function to connect lists of states into an association relationships on the given graph
        # Makes all valid combinations of the subject and object states
        # Creates a new associated state if one cannot be found in existing graph
        
        # Loop through each subject and object state
        for current_sub_state in subject_states:
            for current_obj_state in object_states:
                       
                # Create sorted lists of components/conformations for a possible associated state
                sub_comp, sub_conf = current_sub_state.generate_component_list()
                obj_comp, obj_conf = current_obj_state.generate_component_list()
                associated_component_list = [*sub_comp, *obj_comp]
                associated_conformation_list = [*sub_conf, *obj_conf]
                
                # See if this possible state already exists in the graph
                for current_state in graph.__iter__():
                    if self._state_match_to_component_lists(current_state, associated_component_list, associated_conformation_list, match = 'exact'):
                        # If the state already exists, use it
                        associated_state = current_state
                        break
                
                # If no valid associated state alreay exists, create a new one
                else: # no break
                    associated_state = bkcc.State()
                    associated_state.add_component_list(associated_component_list, associated_conformation_list)
                    
                # Add edges to the associated state from the object and subject states (NetworkX will add any states that don't alreay exist in the graph)
                graph.add_edge(current_sub_state, associated_state)
                graph.add_edge(current_obj_state, associated_state)
                    
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


        
    