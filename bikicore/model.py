"""Class for the model objects used throughout the program.
"""

import uuid
import itertools
import networkx as nx
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
            new_state.req_protein_conf_lists = []
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
                
                # Get a list of accecptable signatures for the rule and read which type of signature we need
                reference_signatures = current_rule.generate_signature_list()
                
                # Find states that fit the rule discription
                matching_subject_states = self._find_states_that_match_rule(current_rule, 'subject')
                matching_object_states = self._find_states_that_match_rule(current_rule, 'object')
                
                # Find the possible pairings of subject and object states that create valid signatures
                possible_state_tuple_list = self._find_association_pairs(reference_signatures, matching_subject_states, matching_object_states)
                
                # Test if the a pair of states could create the implied internal structure required by the rule
                valid_state_tuple_list, valid_link_list = self._find_association_internal_link(current_rule, possible_state_tuple_list)
                # Associate any valid pairs of states 
                for current_state_tuple, current_link_tuple in zip(valid_state_tuple_list, valid_link_list):
                    self._create_association(graph, *current_state_tuple, current_link_tuple)
                
            # Irreversable disassociation
            elif current_rule.rule == ' dissociates_from ':
                
                # Get a list of accecptable signatures for the rule and read which type of signature we need
                reference_signatures = current_rule.generate_signature_list()
                
                # Find states that fit the rule discription
                matching_subject_states = self._find_states_that_match_rule(current_rule, 'subject')
                matching_object_states = self._find_states_that_match_rule(current_rule, 'object')
                
                # Find the possible pairings of subject and object states that create valid signatures
                possible_state_tuple_list = self._find_association_pairs(reference_signatures, matching_subject_states, matching_object_states)
                
                # Test if the a pair of states could create the implied internal structure required by the rule
                valid_state_tuple_list, valid_link_list = self._find_association_internal_link(current_rule, possible_state_tuple_list)
                # Associate any valid pairs of states 
                for current_state_tuple, current_link_tuple in zip(valid_state_tuple_list, valid_link_list):
                    self._create_association(graph, *current_state_tuple, current_link_tuple)
            
            elif current_rule.rule == ' reversibly associates with ':
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
            if self._state_match_to_component_lists(current_state, rule_components, rule_conformations, [], 'minimal'): # We do not worry about links here
                matching_states.append(current_state)

        # Give the list back to the calling method
        return matching_states
   
    def _state_match_to_component_lists(self, query_state, reference_component_list, reference_conformation_list, reference_link_list, match = 'exact'):
        # Helper function that asks if a given query state contains the components given in the reference lists
        # Returns True or False
        # Use modes:
        #   match = 'exact': returns true if and only if the query state contains all the reference components, and no extra components
        #   match = 'minimal': returns true if the query state contains all the reference components, but may have additional ones as well
        
        # NOTE: The link checking won't work properly with "any" ([]) conformation - Try to fix design in the future
        
        # Get the link list that the state requires
        remaining_links = query_state.internal_links.copy()
        
        # Loop through all remaining links for each query link
        for current_link_element1, current_link_element2 in reference_link_list: 
            for i, current_query_links in enumerate(remaining_links):
                
                # Ask if the query and reference in the normal order match
                found_match = False
                if self._compare_components_linked_tuple_element(current_query_links, (current_link_element1, current_link_element2), reference_component_list, reference_conformation_list):
                    found_match = True
                # Ask if the query and reference in the reversed order match
                elif self._compare_components_linked_tuple_element(current_query_links, (current_link_element2, current_link_element1), reference_component_list, reference_conformation_list):
                    found_match = True
                
                #Consume the reference match
                if found_match:
                    del remaining_links[i]
                    break 
            
            else: #no break
                # If we arrive here, there was not a matching link
                return False # Short-circut answer

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
            
        # If we arrive here, matches for all query and references components/conformations/links were found
        # If mode is 'minimal', this means there's a match
        if match == 'minimal':
            return True
        
        # Otherwise, we need to make sure we completely consumed the reference lists
        elif match == 'exact':
            if remaining_components == [] and remaining_conformations == [] and remaining_links == []:
                return True
            else:
                return False
        
        # No other match modes supported, raise error
        else:
            raise ValueError('Function _state_match_to_component_lists received an incorrect argument for parameter "match"')
   
    def _find_association_pairs(self, reference_signatures, matching_subject_states, matching_object_states):
        # Function that returns a list of 2-tuples containing a valid subject and object state pair for an association reaction
        
        # Get the type of signatures required
        count_type = reference_signatures[0].count_type
        
        # Loop through all the possible combinations of subject and object states and store any valid pairs
        valid_pairs = []
        for current_tuple in itertools.product(matching_subject_states, matching_object_states):

            # Create a signature from the current pair
            test_signature = bkcc.CountingSignature(count_type, *current_tuple)
            sub_comp, sub_conf = current_tuple[0].generate_component_list()
            obj_comp, obj_conf = current_tuple[1].generate_component_list()
            test_signature.count_for_third_state(sub_comp + obj_comp, sub_conf + obj_conf)
            
            # See if the test signature matches with any of the reference ones
            if self._signature_match_association(test_signature, reference_signatures):
                valid_pairs.append(current_tuple)
        
        # Give the list of tuples back
        return valid_pairs
    
    def _signature_match_association(self, query_signature, reference_signatures):
        # Helper function that asks if a given query signature matches (any of) the reference signatures in an association context 
        # Returns True or False
        # "reference_signature" can be a single signature or a list of signatures (from "conformation included" signatures computed from rule)
        
        # Make sure we have a list of references
        if not isinstance(reference_signatures, list):
            reference_signatures = [reference_signatures]
        
        # Work through each reference signature until we find one that works
        for refsig in reference_signatures:

            # These are the keys we are looking for, anything else is just extra that we ignore
            refkeys = refsig.third_state_count.keys()
            
            # Check the subject part of the signature
            # The differences between the third and object states should be equal to the reference subject
            query_diff_sub = query_signature.third_state_count - query_signature.object_count
            
            # If we can find a full set of the reference associated state components in a query, we must delete them (from, i.e., a previous dimerization, etc.)
            if all([query_diff_sub[key] >= refsig.third_state_count[key] for key in refkeys]):
                query_diff_sub = query_diff_sub - refsig.third_state_count
            
            # We need to have at least the number of components that the reference state has, but we can have more
            if not all([query_diff_sub[key] >= refsig.subject_count[key] for key in refkeys]):
                continue # Short-circut the comparison, this one doesn't match the pattern
            
            # Check the object part of the signature
            # The differences between the third and subject states should be equal to the reference object
            query_diff_obj = query_signature.third_state_count - query_signature.subject_count
            
            # If we can find a full set of the reference associated state components in a query, we must delete them (from, i.e., a previous dimerization, etc.)
            if all([query_diff_obj[key] >= refsig.third_state_count[key] for key in refkeys]):
                query_diff_obj = query_diff_obj - refsig.third_state_count
            
            # We need to have at least the number of components that the reference state has, but we can have more
            if not all([query_diff_obj[key] >= refsig.object_count[key] for key in refkeys]):
                continue # Short-circut the comparison, this one doesn't match the pattern
            
            # If we arrive here, we have a match and we can short-circut the overall result
            return True
        
        # If we arrive here, nothing matched after exhausting the list of reference signatures
        else:
            return False

    def _create_association(self, graph, subject_state, object_state, new_link):
        # Function to connect two states into an association relationship on the given graph
        # Creates a new associated state if one cannot be found in existing graph
        # New link must be given for components in the assocated state 
                     
        # Create lists of components/conformations for a possible associated state
        sub_comp, sub_conf = subject_state.generate_component_list()
        obj_comp, obj_conf = object_state.generate_component_list()
        associated_component_list = sub_comp + obj_comp
        associated_conformation_list = sub_conf + obj_conf
        associated_old_links = self._combine_internal_link_lists(subject_state, object_state)
        associated_links = associated_old_links + [new_link]
                
        # See if this possible state already exists in the graph
        for current_state in graph.__iter__():
            if self._state_match_to_component_lists(current_state, associated_component_list, associated_conformation_list, associated_links, match = 'exact'):
                # If the state already exists, use it
                associated_state = current_state
                break
        
        # If no valid associated state alreay exists, create a new one
        else: # no break
            associated_state = bkcc.State()
            associated_state.add_component_list(associated_component_list, associated_conformation_list)
            associated_state.internal_links = associated_links
            
        # Add edges to the associated state from the object and subject states (NetworkX will add any states that don't alreay exist in the graph)
        graph.add_edge(subject_state, associated_state)
        graph.add_edge(object_state, associated_state)
            
    def _translate_link_tuple(self, link_tuple, translate_dict):
        # Generator expression to translate complex tuple trees - only run when consumed by tuple later on
        if isinstance(link_tuple, int): # Need to handle a single integer element as well
            children = translate_dict[link_tuple]
            return children
        else:
            children = (self._translate_link_tuple(link_element, translate_dict) if isinstance(link_element, tuple)
                        else translate_dict[link_element] for link_element in link_tuple)
            return tuple(children)
    
    def _combine_component_list(self, list1, list2):
        # Adds component and conformation lists 
         
        # Add drugs and proteins seperately
        drug_list = []
        protein_list = []
        drug_list.extend([x for x in list1 + list2 if isinstance(x, bkcc.Drug)])
        protein_list.extend([x for x in list1 + list2 if isinstance(x, bkcc.Protein)])
        return drug_list + protein_list
    
    def _compare_components_linked_tuple_element(self, query_tuple_element, reference_tuple_element, component_list, conformation_list = []):
        # Generator expression to compare complex tuple trees by referenced component
        
        # NOTE: conformation testing would not work as intended with the "any" ([]) conformation - fix in future redesign
        
        # Recursive testing if we have matching tuples instead of a single index value
        if  isinstance(query_tuple_element, tuple) and isinstance(reference_tuple_element, tuple) and len(query_tuple_element) == len(reference_tuple_element):
            if conformation_list == []: # No conformations
                return all(self._compare_components_linked_tuple_element(*elements, component_list) for elements in zip(query_tuple_element, reference_tuple_element))
            else:
                return all(self._compare_components_linked_tuple_element(*elements, component_list, conformation_list) for elements in zip(query_tuple_element, reference_tuple_element))
        
        # If we have a single index value, test if they match
        elif isinstance(query_tuple_element, int) and isinstance(reference_tuple_element, int):
            if conformation_list == []: # No conformations
                return component_list[query_tuple_element] == component_list[reference_tuple_element]
            else:
                return component_list[query_tuple_element] == component_list[reference_tuple_element] and conformation_list[query_tuple_element] == conformation_list[reference_tuple_element]
        
        # If neither are true, we are by definition a mismatch
        else:
            return False
        
    def _combine_internal_link_lists(self, state_1, state_2, return_translation_dicts = False):
        # Function to translate the link list from state 2 into that from state 1, assuming association of the two states
        
        # How many of each type of component are in each state
        num_1_drug = len(state_1.required_drug_list)
        num_1_protein = len(state_1.required_protein_list)
        num_2_drug = len(state_2.required_drug_list)
        num_2_protein = len(state_2.required_protein_list)
#        num_12_drug = num_1_drug + num_2_drug
#        num_12_protein = num_1_protein + num_2_protein
        
        # Make 1 to 12 translation dictionary
        translate_1_to_12 = {}
        for old_index in range(0, num_1_drug):
            translate_1_to_12[old_index] = old_index
        for old_index in range(num_1_drug, num_1_drug + num_1_protein):
            translate_1_to_12[old_index] = old_index + num_2_drug
            
        # Make 2 to 12 translation dictionary
        translate_2_to_12 = {}
        for old_index in range(0, num_2_drug):
            translate_2_to_12[old_index] = old_index + num_1_drug
        for old_index in range(num_2_drug, num_2_drug + num_2_protein):
            translate_2_to_12[old_index] = old_index + num_1_drug + num_1_protein
        
        # Translate the link lists
        translated_link_1_list = [self._translate_link_tuple(link, translate_1_to_12) for link in state_1.internal_links]
        translated_link_2_list = [self._translate_link_tuple(link, translate_2_to_12) for link in state_2.internal_links]
        translated_list = translated_link_1_list + translated_link_2_list
        
        # Return the requested information
        if return_translation_dicts:
            # Make 12 to 1 and 12 to 2 dictionaries
            translate_12_to_1 = {value: key for key, value in translate_1_to_12.items()}
            translate_12_to_2 = {value: key for key, value in translate_2_to_12.items()}
            
            return translated_list, [translate_1_to_12, translate_2_to_12, translate_12_to_1, translate_12_to_2]
        else:
            return translated_list

    def _find_association_internal_link(self, rule, state_pairs):
        # Function to check if the pairs of states given can be associated to make the new internal link structure implied by the given rule
        
        # Check each pair of states
        valid_state_tuples = []
        valid_link_tuples = []
        for state_pair in state_pairs:
            
            # We may have more than one component (or set of components) in each state that matches the rule, so find any that could apply 
            matching_subject_indices = [[] for x in range(len(rule.rule_subject))]
            for rule_comp_index, subject_matched_list in enumerate(matching_subject_indices):
                for state_comp_index, (current_comp, current_conf) in enumerate(zip(*state_pair[0].generate_component_list())):
                    
                    # Record if the state and rule components match
                    if current_comp == rule.rule_subject[rule_comp_index] and (current_conf == rule.subject_conf[rule_comp_index] or rule.subject_conf[rule_comp_index] == []):
                        subject_matched_list.append(state_comp_index)
                        
            # We may have more than one component in each state that matches the rule, so find any that could apply 
            matching_object_indices = [[] for x in range(len(rule.rule_object))]
            for rule_comp_index, object_matched_list in enumerate(matching_object_indices):
                for state_comp_index, (current_comp, current_conf) in enumerate(zip(*state_pair[1].generate_component_list())):
        
                    # Record if the state and rule components match
                    if current_comp == rule.rule_object[rule_comp_index] and (current_conf == rule.object_conf[rule_comp_index] or rule.object_conf[rule_comp_index] == []):
                        object_matched_list.append(state_comp_index)
                        
            # Make all possible combinations of links with the found indices
            subject_combos = itertools.product(*matching_subject_indices)
            object_combos = itertools.product(*matching_object_indices)                    
            
            # Make the associated link index list, translation dictionaries, and component list
            associated_old_link_list, translation_dicts = self._combine_internal_link_lists(state_pair[0], state_pair[1], True)
            associated_component_list = self._combine_component_list(state_pair[0].generate_component_list(True), state_pair[1].generate_component_list(True))
            
            # Now ask if we previously made this exact same link, and if so, reject the pair
            valid_first_elements = []
            valid_second_elements = []

            #for newlink in itertools.product(subject_combos, object_combos):
            for new_link_first_element, new_link_second_element in itertools.product(subject_combos, object_combos):
 
                # If the element is a singleton tuple, convert it to a number
                if len(new_link_first_element) == 1:
                    new_link_first_element = new_link_first_element[0]
                if len(new_link_second_element) == 1:
                    new_link_second_element = new_link_second_element[0]
                    
                
                # Translate into associated state indices
                new_associated_link_first_element = self._translate_link_tuple(new_link_first_element, translation_dicts[0])
                new_associated_link_second_element = self._translate_link_tuple(new_link_second_element, translation_dicts[1])
               
                # Check if the first part of the new link is already linked with another component copy of the second part 
                # Find anything linked to the first element
                old_link_element_to_first_element = []
                for old_link in associated_old_link_list:
                    if old_link[0] == new_associated_link_first_element:
                        old_link_element_to_first_element.append(old_link[1])
                    if old_link[1] == new_associated_link_first_element:
                        old_link_element_to_first_element.append(old_link[0])
                
                # Check if anything has the same component structure as the new peice
                old_link_found = []
                for old_link_element in old_link_element_to_first_element:
                    old_link_found.append(self._compare_components_linked_tuple_element(old_link_element, new_associated_link_second_element, associated_component_list))
                        
                # Find anything linked to the second element
                old_link_element_to_second_element = []
                for old_link_element in associated_old_link_list:
                    if old_link[0] == new_associated_link_second_element:
                        old_link_element_to_second_element.append(old_link[1])
                    if old_link[1] == new_associated_link_second_element:
                        old_link_element_to_second_element.append(old_link[0])
                
                # Check if anything has the same component structure as the new peice
                for old_link_element in old_link_element_to_second_element:
                    old_link_found.append(self._compare_components_linked_tuple_element(old_link_element, new_associated_link_first_element, associated_component_list))
            
                # We should be able to use any links that haven't been formed already
                if not any(old_link_found):
                    valid_first_elements.append(new_associated_link_first_element)
                    valid_second_elements.append(new_associated_link_second_element)
            
            # Any valid link element combos should be made into a new link
            for first_element, second_element in zip(valid_first_elements, valid_second_elements):
                valid_state_tuples.append(state_pair)
                valid_link_tuples.append((first_element, second_element))

        # Send back the list of things that should work
        return valid_state_tuples, valid_link_tuples
    
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


        
    