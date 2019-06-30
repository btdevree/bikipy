"""Class for the model objects used throughout the program.
"""

import uuid
import itertools
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
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
    
    def generate_network(self, max_cycles=20):
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
            for current_conformation in range(len(current_component.conformation_names)):
                new_state = bkcc.State()
                new_state.required_drug_list = []
                new_state.required_protein_list = [current_component]
                new_state.req_protein_conf_lists = [[current_conformation]]
                self.network.main_graph.add_node(new_state)
        
        # After creating the singleton graph, reapply the network rules until no more changes happen
        old_network = self.network.main_graph.copy()
        current_cycle_number = 0
        while current_cycle_number <= max_cycles:
            self.apply_rules_to_network()
#            self._main_graph_dump('Graph_0')
#            print('Graph size = ', self.network.main_graph.number_of_nodes())
            # Check if graph changed
            if nx.algorithms.isomorphism.is_isomorphic(old_network, self.network.main_graph):
#                self._main_graph_dump('Final_graph')
#                print('Graph size = ', self.network.main_graph.number_of_nodes())
                break
            else:
                old_network = self.network.main_graph.copy()
                current_cycle_number += 1
#                self._main_graph_dump('Graph_' + str(current_cycle_number))
#                print('Graph size = ', self.network.main_graph.number_of_nodes())
        
    def apply_rules_to_network(self, graph = None):
        # Apply the model's rules to an existing graph

        # If default, work on the main graph
        if graph == None:
            graph = self.network.main_graph
        
        # Apply each rule to the graph
        for current_rule in self.rule_list:

            # Each type of rule needs a different treatment
            
            # Association
            if current_rule.rule == ' associates with ' or current_rule.rule == ' reversibly associates with ' \
                    or current_rule.rule == ' associates and dissociates in rapid equlibrium with ':
                
                # Get a list of accecptable signatures for the rule and read which type of signature we need
                reference_signatures = current_rule.generate_signature_list()
                
                # Find states that fit the rule description
                matching_subject_states = self._find_states_that_match_rule(current_rule, 'subject')
                matching_object_states = self._find_states_that_match_rule(current_rule, 'object')
                
                # Find the possible pairings of subject and object states that create valid signatures
                possible_state_tuple_list = self._find_association_pairs(reference_signatures, matching_subject_states, matching_object_states)
                
                # Test if the a pair of states could create the implied internal structure required by the rule
                valid_state_tuple_list, valid_link_list = self._find_association_internal_link(current_rule, possible_state_tuple_list)
                
                # Associate any valid pairs of states 
                for current_state_tuple, current_link_tuple in zip(valid_state_tuple_list, valid_link_list):
                    if current_rule.rule == ' associates with ':
                        self._create_association(graph, *current_state_tuple, current_link_tuple)
                    elif current_rule.rule == ' reversibly associates with ':
                        self._create_association(graph, *current_state_tuple, current_link_tuple, reversible = True)
                    elif current_rule.rule == ' associates and dissociates in rapid equlibrium with ':
                        self._create_association(graph, *current_state_tuple, current_link_tuple, reversible = True)
         
            # Dissociation
            elif current_rule.rule == ' dissociates from ' or current_rule.rule == ' reversibly dissociates from ' \
                    or current_rule.rule == ' dissociates and reassociates in rapid equlibrium from ':

                # Get a list of accecptable signatures for the rule and read which type of signature we need
                reference_signatures = current_rule.generate_signature_list()
                
                # Find states that fit the rule description
                matching_object_states = self._find_states_that_match_rule(current_rule, 'object')
                
                # Find the possible pairings of subject and object states that create valid signatures
                possible_state_split_list = self._find_dissociation_pairs(reference_signatures, matching_object_states)
                
                # Test if the a pair of states could create the implied internal structure required by the rule
                valid_state_split_list, valid_link_lists = self._find_dissociation_internal_link(current_rule, possible_state_split_list)
                
                # Associate any valid pairs of states 
                for current_state_split_tuple, current_link_tuple in zip(valid_state_split_list, valid_link_lists):
                    if current_rule.rule == ' dissociates from ':
                        self._create_dissociation(graph, *current_state_split_tuple, *current_link_tuple)
                    elif current_rule.rule == ' reversibly dissociates from ':
                        self._create_dissociation(graph, *current_state_split_tuple, *current_link_tuple, reversible = True)
                    elif current_rule.rule == ' dissociates and reassociates in rapid equlibrium from ':
                        self._create_dissociation(graph, *current_state_split_tuple, *current_link_tuple, reversible = True)
            
            # Conformational changes and reactions
            elif current_rule.rule == ' converts to ' or current_rule.rule == ' reversibly converts to ' \
                    or current_rule.rule == ' converts in rapid equlibrium to ':
                
                # Get a list of acceptable signatures for the rule and read which type of signature we need
                # Add converstion signatures
                reference_signatures = current_rule.generate_signature_list()
                
                # Find states that fit the rule description
                matching_subject_states = self._find_states_that_match_rule(current_rule, 'subject')
                
                # Find all possible conversion reactions with the matching states, returns the components involved and what they change to, but not the internal structure
                possible_conversion_tuples = self._find_conversion_pairs(current_rule, reference_signatures, matching_subject_states)

                # Validate links for possible conversion reactions
                valid_conversion_tuples = self._find_conversion_internal_link(current_rule, possible_conversion_tuples)
                
                # Make the conversion
                for convert_tuple in valid_conversion_tuples:
                    if current_rule.rule == ' converts to ':
                        self._create_conversion(graph, *convert_tuple)
                    if current_rule.rule == ' reversibly converts to ':
                        self._create_conversion(graph, *convert_tuple, reversible = True) 
                    if current_rule.rule == ' converts in rapid equlibrium to ':
                        self._create_conversion(graph, *convert_tuple, reversible = True)
            
            # Competition rule
            elif current_rule.rule == ' is competitive with ':   
                
                # Find states that fit the rule description
                matching_subject_states = self._find_states_that_match_rule(current_rule, 'both')
               
                # Validate links for matching states
                states_to_remove = self._find_competition_internal_link(current_rule, matching_subject_states)
                
                # Make the conversion
                for current_state in states_to_remove:
                    self._remove_state(graph, states_to_remove)
            
            else:
                raise ValueError("Rule not recognized")
    
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
        query_component_list, query_conformation_list = query_state.generate_component_list()
        
        # Run the comparison
        component_boolean = self._compare_component_lists(query_component_list, query_conformation_list, reference_component_list, reference_conformation_list, match)
        
        # Need to factor in the link matching as well
        if match == 'minimal':
            return component_boolean
        elif match == 'exact':
            if remaining_links == []:
                return component_boolean
            else:
                return False
        else: # No other match modes supported, raise error
            raise ValueError('Function _state_match_to_component_lists received an incorrect argument for parameter "match"')
        
    def _compare_component_lists(self, query_component_list, query_conformation_list, reference_component_list, reference_conformation_list, match = 'exact'):
        # Function to compare two sets of component and conformation lists, with no regard to order and without replacement    
        
        # Create copies of query lists for consumption
        remaining_components = query_component_list.copy()
        remaining_conformations = query_conformation_list.copy()
        
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
            if remaining_components == [] and remaining_conformations == []:
                return True
            else:
                return False
        
        # No other match modes supported, raise error
        else:
            raise ValueError('Function _compare_component_lists received an incorrect argument for parameter "match"')
    
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
      
    def _find_dissociation_pairs(self, reference_signatures, matching_object_states):
        # Function that returns a list of 2-tuples containing a tuple of indices to split off as subject and a object state for a dissociation reaction
        
        # Get the type of signatures required
        count_type = reference_signatures[0].count_type
        
        # Loop through all the found object states and store any valid splitting of components
        valid_pairs = []
        for current_obj in matching_object_states:
            obj_comp_list, obj_conf_list = current_obj.generate_component_list()
            
            # Need to split the object state into all possible subject & third states
            split_indices = []
            for num_indices in range(1, len(obj_comp_list)): # Need to repeat for single, double, etc. allowed number of conformations
                split_indices.extend(itertools.combinations(range(0, len(obj_comp_list)), num_indices))
            
            # Create a signature from each possible splitting
            for current_split in split_indices:
                test_signature = bkcc.CountingSignature(count_type, object_state = current_obj)
                test_sub_comp = [obj_comp_list[x] for x in current_split]
                test_sub_conf = [obj_conf_list[x] for x in current_split]
                test_signature.count_for_subject(test_sub_comp, test_sub_conf)
                test_signature.third_state_count = test_signature.object_count - test_signature.subject_count
                
                # See if the test signature matches with any of the reference ones
                if self._signature_match_dissociation(test_signature, reference_signatures):
                    valid_pairs.append((current_obj, current_split))
        
        # Give the list of tuples back
        return valid_pairs
          
    def _signature_match_dissociation(self, query_signature, reference_signatures):
        # Helper function that asks if a given query signature matches (any of) the reference signatures in dissociation context 
        # Returns True or False
        # "reference_signature" can be a single signature or a list of signatures (from "conformation included" signatures computed from rule)
        
        # Make sure we have a list of references
        if not isinstance(reference_signatures, list):
            reference_signatures = [reference_signatures]
        
        # Work through each reference signature until we find one that works
        for refsig in reference_signatures:

            # These are the keys we are looking for, anything else is just extra that we ignore
            refkeys = refsig.object_count.keys()
            
            # Check if the signature matches the rule, and if so, return true
            sub_test = [query_signature.subject_count[key] == refsig.subject_count[key] for key in refkeys]
            obj_test = [query_signature.object_count[key] >= refsig.object_count[key] for key in refkeys]
            third_test = [query_signature.third_state_count[key] >= refsig.third_state_count[key] for key in refkeys]
            if all(sub_test + obj_test + third_test):
                return True
        
        # If we arrive here, nothing matched after exhausting the list of reference signatures
        else:
            return False
    
    def _find_conversion_pairs(self, rule, reference_signatures, matching_subject_states):
        # Function that returns a list of 4-tuples containing the subject state, indices of the changed components, and lists of converted components and conformations.
        
        # Get the type of signatures required
        count_type = reference_signatures[0].count_type
        convert_rule_components, convert_rule_conformations = rule.generate_component_list('object')
        
        # Loop through all the found object states and store any valid splitting of components
        valid_pairs = []
        for current_sub in matching_subject_states:
            sub_comp_list, sub_conf_list = current_sub.generate_component_list()

            # Look for matches with each conversion signature
            for current_sig in reference_signatures:
                
                # Test all possible sets of indices for the conversion
                number_indices_needed = len([*current_sig.subject_count.elements()])
                for current_indices in itertools.combinations(range(0, len(sub_comp_list)), number_indices_needed):
                    
                    # Get lists of anything that's not tagged as being converted with the indices
                    nonconvert_comp_list = [x for i, x in enumerate(sub_comp_list) if i not in current_indices]
                    nonconvert_conf_list = [x for i, x in enumerate(sub_conf_list) if i not in current_indices]
                    
                    # Make lists of new components/conformations and see if they create a matching signature
                    ignored_output, convert_conf_lists = rule._get_all_conformation_combinations([], convert_rule_conformations, [], convert_rule_components)   
                    for current_convert_conf in convert_conf_lists:
                        
                        # Make a new signature with the converted components
                        test_signature = bkcc.CountingSignature(count_type, subject_state = current_sub) 
                        test_signature.count_for_object(nonconvert_comp_list + convert_rule_components, nonconvert_conf_list + current_convert_conf) # Add the converted components to the signature
                        
                        # See if the test signature matches the reference
                        if self._signature_match_conversion(test_signature, current_sig):
                            valid_pairs.append((current_sub, current_indices, convert_rule_components, current_convert_conf))
        
        # Give the list of tuples back
        return valid_pairs
        
    def _signature_match_conversion(self, query_signature, reference_signature):
        # Helper function that asks if a given query signature matches the reference signature in a conversion context
        # Returns True or False
       
        # Get the difference between reference states
        ref_fwd_diff = reference_signature.subject_count - reference_signature.object_count
        ref_rev_diff = reference_signature.object_count - reference_signature.subject_count
        
        # Get the difference between query states states
        query_fwd_diff = query_signature.subject_count - query_signature.object_count
        query_rev_diff = query_signature.object_count - query_signature.subject_count
        
        if ref_fwd_diff == query_fwd_diff and ref_rev_diff == query_rev_diff:
            return True
        else:
            return False
    
    def _create_association(self, graph, subject_state, object_state, new_link, reversible = False):
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
        if reversible:
            graph.add_edge(associated_state, subject_state)
            graph.add_edge(associated_state, object_state)
    
    def _create_dissociation(self, graph, object_state, split_indices, subject_link_list, third_state_link_list, reversible = False):
        # Function to split an object state into the subject state and a remaining third state on a given graph
        # Creates new states if the generated ones cannot be found in existing graph
        # Link lists must be given for components in the split states
                
        # Get indices of the third state 
        remaining_indices = [x for x in range(0, len(object_state.required_drug_list + object_state.required_protein_list)) if x not in split_indices] 
        
        # Create new lists of components for the split-off states 
        object_comp, object_conf = object_state.generate_component_list()
        subject_comp = [object_comp[i] for i in split_indices]
        subject_conf = [object_conf[i] for i in split_indices]
        third_state_comp = [object_comp[i] for i in remaining_indices]
        third_state_conf = [object_conf[i] for i in remaining_indices]
                
        # See if the subject state already exists in the graph
        for current_state in graph.__iter__():
            if self._state_match_to_component_lists(current_state, subject_comp, subject_conf, subject_link_list, match = 'exact'):
                # If the state already exists, use it
                subject_state = current_state
                break
        
        # If no valid subject state alreay exists, create a new one
        else: # no break
            subject_state = bkcc.State()
            subject_state.add_component_list(subject_comp, subject_conf)
            subject_state.internal_links = subject_link_list
        
        # See if the third state already exists in the graph
        for current_state in graph.__iter__():
            if self._state_match_to_component_lists(current_state, third_state_comp, third_state_conf, third_state_link_list, match = 'exact'):
                # If the state already exists, use it
                third_state = current_state
                break
        
        # If no valid subject state alreay exists, create a new one
        else: # no break
            third_state = bkcc.State()
            third_state.add_component_list(third_state_comp, third_state_conf)
            third_state.internal_links = third_state_link_list
            
        # Add edges to the associated state from the object and subject states (NetworkX will add any states that don't alreay exist in the graph)
        graph.add_edge(object_state, subject_state)
        graph.add_edge(object_state, third_state)
        if reversible:
             graph.add_edge(subject_state, object_state)
             graph.add_edge(third_state, object_state)    

    def _create_conversion(self, graph, subject_state, new_component_list, new_conformation_list, new_link_tuples, reversible = False):
        # Function to convert the given components in the state to those specified by the rule
        # Creates new states if the generated ones cannot be found in existing graph
        
        # See if the object state already exists in the graph
        for current_state in graph.__iter__():
            if self._state_match_to_component_lists(current_state, new_component_list, new_conformation_list, new_link_tuples, match = 'exact'):
                # If the state already exists, use it
                object_state = current_state
                break
        
        # If no valid subject state alreay exists, create a new one
        else: # no break
            object_state = bkcc.State()
            object_state.add_component_list(new_component_list, new_conformation_list)
            object_state.internal_links = new_link_tuples
            
        # Add edges to convert states (NetworkX will add any states that don't alreay exist in the graph)
        graph.add_edge(subject_state, object_state)
        if reversible:
             graph.add_edge(object_state, subject_state)          

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
    
    def _compare_components_dissociation_rule_and_link(self, rule, link, reference_state):
        # Tests if the broken link in the given state matches the rule's splitting pattern
        
        # Get lists of components/conformations
        reference_component_list, reference_conformation_list = reference_state.generate_component_list()
        link_element1_component_list, link_element1_conformation_list = self._collect_link_components(link[0], reference_component_list, reference_conformation_list)
        link_element2_component_list, link_element2_conformation_list = self._collect_link_components(link[1], reference_component_list, reference_conformation_list)
        rule_subject_comp, rule_subject_conf = rule.generate_component_list('subject')
        rule_third_state_comp, rule_third_state_conf = rule.generate_component_list('difference')
        
        # See if the link matches the rule in either forward or reverse direction
        if self._compare_component_lists(link_element1_component_list, link_element1_conformation_list, rule_subject_comp, rule_subject_conf, match = 'minimal'):
            if self._compare_component_lists(link_element2_component_list, link_element2_conformation_list, rule_third_state_comp, rule_third_state_conf, match = 'minimal'):
                return True
        elif self._compare_component_lists(link_element1_component_list, link_element1_conformation_list, rule_third_state_comp, rule_third_state_conf, match = 'minimal'):
            if self._compare_component_lists(link_element2_component_list, link_element2_conformation_list, rule_subject_comp, rule_subject_conf, match = 'minimal'):
                return True
        else:
            return False
    
    def _collect_link_components(self, query_tuple_element, reference_component_list, reference_conformation_list, previous_component_list = None, previous_conformation_list = None):
        # Generator expression to collect component/conformation lists from complex tuple trees
        
        # Need to make new list if none passed in call
        if previous_component_list == None:
            previous_component_list = []
        if previous_conformation_list == None:
            previous_conformation_list = []
        
        # Recursive testing if we have a tuple of values single index value
        if isinstance(query_tuple_element, tuple):
            for element in query_tuple_element:
                previous_component_list, previous_conformation_list = self._collect_link_components(element, reference_component_list, reference_conformation_list, previous_component_list, previous_conformation_list) 
    
        # Otherwise, get the referenced component/conformation by the index value 
        elif isinstance(query_tuple_element, int):
            previous_component_list.append(reference_component_list[query_tuple_element])
            previous_conformation_list.append(reference_conformation_list[query_tuple_element])
        
        # Don't want to silently accecpt bad argurments
        else:
            raise ValueError("Function received unexpected values in the query_tuple_element parameter")
        
        return previous_component_list, previous_conformation_list
    
    def _collect_link_values(self, tuple_element, previous_element_list = None):
        # Generator expression to collect index values from complex tuple trees
        
        # Need to make new list if none passed in call
        if previous_element_list == None:
            previous_element_list = []            
        
        # Recursive testing if we have a tuple of values single index value
        if isinstance(tuple_element, tuple) or isinstance(tuple_element, list):
            for element in tuple_element:
                previous_element_list = self._collect_link_values(element, previous_element_list) 
    
        # Otherwise, get the referenced component/conformation by the index value 
        elif isinstance(tuple_element, int):
            previous_element_list.append(tuple_element)
            
        # Don't want to silently accecpt bad argurments
        else:
            raise ValueError("Function received unexpected values of type {} in the tuple_element parameter".format(str(type(tuple_element))))
        
        return previous_element_list
            
    def _combine_internal_link_lists(self, state_1, state_2, return_translation_dicts = False):
        # Function to translate the link lists from states 1 and 2 into that of the associated state
        
        # How many of each type of component are in each state
        num_1_drug = len(state_1.required_drug_list)
        num_1_protein = len(state_1.required_protein_list)
        num_2_drug = len(state_2.required_drug_list)
        num_2_protein = len(state_2.required_protein_list)
        
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
        
    def _split_internal_link_list(self, state_12, split_indices, return_translation_dicts = False):
        # Function to translate the link list from an associated state into two seperate states.
        
        # Calculate which links are going to be broken
        broken_state12_links = []
        state1_links = []
        state2_links = []
        for current_link in state_12.internal_links:
            
            # Run tests to see if link contains elements from the states that are split off
            all_test = self._index_tuple_element_all_split_test(current_link, split_indices)
            any_test = self._index_tuple_element_any_split_test(current_link, split_indices)
            
            # Classify link
            if all_test: # Only a link completely contained in the split-off state will pass
                state1_links.append(current_link)
            elif not any_test: # Only a link completely contained in the remaining state will fail
                state2_links.append(current_link)
            elif not all_test and any_test: # Some, but not all, the components of the link are split, so this is broken
                broken_state12_links.append(current_link)
       
        # Translate the remaining links for the new states
        # How many of each type of component are in the associated state
        num_12_drug = len(state_12.required_drug_list)
        num_12_protein = len(state_12.required_protein_list)
        
        # Get indices of the third state 
        remaining_indices = [x for x in range(0, num_12_drug + num_12_protein) if x not in split_indices] 
        
        # Make 12 to 1/2 translation dictionary, 1 = subject, 2 = third state
        translate_12_to_1 = {}
        for new_index, old_index in enumerate(split_indices):
            translate_12_to_1[old_index] = new_index
        translate_12_to_2 = {}
        for new_index, old_index in enumerate(remaining_indices):
            translate_12_to_2[old_index] = new_index
            
        # Translate the link lists
        translated_state1_links = [self._translate_link_tuple(link, translate_12_to_1) for link in state1_links]
        translated_state2_links = [self._translate_link_tuple(link, translate_12_to_2) for link in state2_links]
        
        # Return the requested information
        if return_translation_dicts:
            # Make 12 to 1 and 12 to 2 dictionaries
            translate_1_to_12 = {value: key for key, value in translate_12_to_1.items()}
            translate_2_to_12 = {value: key for key, value in translate_12_to_2.items()}
            
            return broken_state12_links, translated_state1_links, translated_state2_links, [translate_1_to_12, translate_2_to_12, translate_12_to_1, translate_12_to_2]
        else:
            return broken_state12_links, translated_state1_links, translated_state2_links
    
    def _index_tuple_element_all_split_test(self, index_tuple_element, split_indices):
        # Generator expression to compare complex tuple trees, asks if all of the index values of the tuple element are included in the split list
                
        # Recursive testing if we have tuples instead of a single index value
        if  isinstance(index_tuple_element, tuple):
            return all(self._index_tuple_element_any_split_test(element, split_indices) for element in index_tuple_element)
        else:
            return index_tuple_element in split_indices
    
    def _index_tuple_element_any_split_test(self, index_tuple_element, split_indices):
        # Generator expression to compare complex tuple trees, asks if any of the index values of the tuple element are included in the split list
                
        # Recursive testing if we have tuples instead of a single index value
        if  isinstance(index_tuple_element, tuple):
            return any(self._index_tuple_element_any_split_test(element, split_indices) for element in index_tuple_element)
        else:
            return index_tuple_element in split_indices
        
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
                
    def _find_dissociation_internal_link(self, rule, state_split_pairs):
        # Function to check if the state and split pair given can be dissociated to leave a new internal link structure implied by the given rule
        
        # Check each pair of states and splits
        valid_state_split_tuples = []
        valid_link_tuples = []
        for state_split_pair in state_split_pairs:
            
            # Try making the split and get back the broken and changed link lists
            broken_links, new_subject_links, new_third_state_links = self._split_internal_link_list(state_split_pair[0], state_split_pair[1])
            
            # Any split that breaks more than one link is not appropreate (Is this true? I think so...)
            if len(broken_links) != 1:
                continue # This is not correct, try next splitting position           
            
            # See if the broken link (in the forward direction) matches the rule's splitting pattern
            if self._compare_components_dissociation_rule_and_link(rule, broken_links[0], state_split_pair[0]):
                
                # Add split to the list of valid reactions            
                valid_state_split_tuples.append(state_split_pair)
                valid_link_tuples.append((new_subject_links, new_third_state_links))

        # Send back the list of things that should work
        return valid_state_split_tuples, valid_link_tuples
    
    def _find_conversion_internal_link(self, rule, possible_conversion_tuples):
    # Function to check if the internal link structure of indicated conversion could be described by rule.

        # Work through each tuple of possible solutions and figure out if the internal links could work.
        valid_conversion_tuples = []
        for current_possible_tuple in possible_conversion_tuples:
            subject_state, conversion_indices, converted_components, converted_conformations = current_possible_tuple
            
            # Can only deal with 1 to 1 component conversion for the moment
            if not len(conversion_indices) == len(converted_components):
                raise NotImplementedError("Conversion rules not yet able to handle reactions with changing numbers of components")

            # The components being converted must all be somehow linked directly to match the rule
            # Split the state to get the converted components by themselves
            broken_links, conversion_links, nonconvert_links = self._split_internal_link_list(subject_state, conversion_indices)
            
            # See if the conversion links have at least one mention of each component
            if len(conversion_indices) > 1: # Single component conversions don't need to pass this test
                if not all([x in self._collect_link_values(conversion_links) for x in conversion_indices]):
                    continue # Failed, short-circult test to try next tuple of possible conversions
            
            # Test all possible permutations of the index order, we don't know which one is right yet
            index_translations = itertools.permutations(conversion_indices)
            for current_translation_indices in index_translations:
                 
                # Get lists of components/conformations to compare
                new_comp_list, new_conf_list = subject_state.generate_component_list()
                rule_comp_list, rule_conf_list = rule.generate_component_list('subject')
                
                # See if all the ordered components in the rule match with those from the state
                match_list = []
                for rule_index, state_index in enumerate(current_translation_indices):
                    match_list.append(new_comp_list[state_index] == rule_comp_list[rule_index])
                    if not rule_conf_list[rule_index] == []: # Just ignore the state's conformation if there is an "any" conformation in the rule
                        match_list.append(new_conf_list[state_index] == rule_conf_list[rule_index])
                
                # If they all match, make the converstion and add to the converstion tuple list
                if all(match_list):
                    converted_rule_comp_list, converted_rule_conf_list = rule.generate_component_list('object')
                    for rule_index, state_index in enumerate(current_translation_indices):
                        new_comp_list[state_index] = converted_rule_comp_list[rule_index]
                        if not converted_rule_conf_list[rule_index] == []: # Just ignore the state's conformation if there is an "any" conformation in the rule
                            new_conf_list[state_index] = converted_rule_conf_list[rule_index]
                    
                    # Add to the working converstion tuples, links not changed yet (but might be if we implement non-1:1 conversions)
                    valid_conversion_tuples.append((subject_state, new_comp_list, new_conf_list, subject_state.internal_links))
                
        # Give back the list of stuff that should be working
        return valid_conversion_tuples
    
    def _find_competition_internal_link(self, rule, matching_subject_states):
    # Function to check if the internal link structure of indicated competitive states is valid.

        # Work through each tuple of possible solutions and figure out if the internal links could match.
        valid_competitive_states = []
        for current_state in matching_subject_states:
            
            # Which components are referenced by the rule
            rule_comp_sub, rule_conf_sub = rule.generate_component_list('subject')
            rule_comp_obj, rule_conf_obj = rule.generate_component_list('object')
            
            # Test all combinations of subject and object states
            = 
            
            # The competitive comonents must all be linked directly to the same thing
            # Split the state to get the converted components by themselves
            broken_links, conversion_links, nonconvert_links = self._split_internal_link_list(subject_state, conversion_indices)
            
            # See if the conversion links have at least one mention of each component
            if len(conversion_indices) > 1: # Single component conversions don't need to pass this test
                if not all([x in self._collect_link_values(conversion_links) for x in conversion_indices]):
                    continue # Failed, short-circult test to try next tuple of possible conversions
            
            # Test all possible permutations of the index order, we don't know which one is right yet
            index_translations = itertools.permutations(conversion_indices)
            for current_translation_indices in index_translations:
                 
                # Get lists of components/conformations to compare
                new_comp_list, new_conf_list = subject_state.generate_component_list()
                rule_comp_list, rule_conf_list = rule.generate_component_list('subject')
                
                # See if all the ordered components in the rule match with those from the state
                match_list = []
                for rule_index, state_index in enumerate(current_translation_indices):
                    match_list.append(new_comp_list[state_index] == rule_comp_list[rule_index])
                    if not rule_conf_list[rule_index] == []: # Just ignore the state's conformation if there is an "any" conformation in the rule
                        match_list.append(new_conf_list[state_index] == rule_conf_list[rule_index])
                
                # If they all match, make the converstion and add to the converstion tuple list
                if all(match_list):
                    converted_rule_comp_list, converted_rule_conf_list = rule.generate_component_list('object')
                    for rule_index, state_index in enumerate(current_translation_indices):
                        new_comp_list[state_index] = converted_rule_comp_list[rule_index]
                        if not converted_rule_conf_list[rule_index] == []: # Just ignore the state's conformation if there is an "any" conformation in the rule
                            new_conf_list[state_index] = converted_rule_conf_list[rule_index]
                    
                    # Add to the working converstion tuples, links not changed yet (but might be if we implement non-1:1 conversions)
                    valid_conversion_tuples.append((subject_state, new_comp_list, new_conf_list, subject_state.internal_links))
                
        # Give back the list of stuff that should be working
        return valid_conversion_tuples
    
    # Graphing utility function
    def _main_graph_dump(self, label='default'):
        plt.cla()
        G = self.network.main_graph
        nx.draw(G, with_labels=True, font_weight='bold')
        plt.savefig(label+'.png')
        
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

