"""Classes representing physical or conceptual objects needed to describe the
biochemical system, the model, and the experiments performed. 

"""

import uuid
import itertools
import networkx as nx
import sympy as sp
from sympy.core.basic import Basic as spBaseClass 
from collections import Counter
from traits.api import HasTraits, Str, List, Tuple, Int, Instance, Enum, Either, Bool
from bikipy.bikicore.exceptions import ComponentNotValidError, RuleNotValidError

# Define classes for the different components of the biochemical system
class Drug(HasTraits):
    
    # Traits initialization
    name = Str()
    symbol = Str()
    ID = Instance(uuid.UUID) 
    
    # Want to give a new state an ID right away
    def __init__(self, *args, **kwargs):
        self.ID = uuid.uuid4()
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery 
    
class Protein(HasTraits):
    
    # Traits initialization
    name = Str()
    symbol = Str()
    conformation_names = List(Str()) 
    conformation_symbols = List(Str())
    ID = Instance(uuid.UUID) 
    
    # Want to give a new state an ID right away
    def __init__(self, *args, **kwargs):
        self.ID = uuid.uuid4()
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery 
    
    # Check if the trait values are valid. Use after user editing of the object, for example.
    def check_protein_traits(self):
        
        # Check if the number of conformation names and symbols are the same
        if len(self.conformation_names) != len(self.conformation_symbols):
            raise ComponentNotValidError('The number of protein conformation names and symbols must be the same')

# Define class for the different states - nodes on network graph
class State(HasTraits):
    
    # Traits initialization
    name = Str()
    symbol = Str()
    number = Int()
    variable = Instance(spBaseClass) # Must be a sympy object
    ID = Instance(uuid.UUID)
    required_drug_list = List(Instance(Drug))
    required_protein_list = List(Instance(Protein))
    req_protein_conf_lists = List(List(Int()))
    internal_links = List(Tuple())
    
    # Want to give a new state an ID right away
    def __init__(self, *args, **kwargs):
        self.ID = uuid.uuid4()
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery 
    
    def autosymbol(self, method = 0):
        # Create the generated symbol name that will usually be the one used in display. 
        # Symbols are created from the components of the state
        # Parameters:
            # method - integer, default = 0. Choose one of the different methods of sorting the components

        # Get the component lists of the state - NOTE: This is stupid, i unzip then rezip. Future cleanup.
        (indices, components, conformations) = zip(*self.enumerate_components(True))

        # Order list using method 0
        if method == 0:
            # We will sort the proteins by conformation, then move any linked proteins together,
            # then move drugs next to their linked protein. 
            
            # Make some sorting tuples
            tuple_list = [*zip(components, indices)]

            # Step 1 - Sort proteins by conformation
            # Put in standard timsort order
            tuple_list.sort(key = lambda x: x[0].symbol)
            
            # Divide list into sublists with the same components
            group_tuples = []
            for group_key, group_iterator in itertools.groupby(tuple_list, key = lambda x: x[0]): # Group by the first element of the tuple
                group_tuples.append(list(group_iterator))
            
            # Work on one group at a time
            sorted_groups = []
            for current_group in group_tuples:

                # Insertion sort the group of components
                group_tuple_list = []
                for current_tuple in current_group:

                    # Get current conformation
                    current_conf = conformations[current_tuple[1]]
                    
                    # Compare tuple against elements until we find one that belongs behind it
                    for test_position in range(len(group_tuple_list)):
                        
                        # Get test conformation
                        test_conf = conformations[group_tuple_list[test_position[1]]]
                        
                        # Insert shorter conformations lists before the long one
                        if len(current_conf) < len(test_conf):
                            group_tuple_list.insert(test_position, current_tuple)
                            break
                        
                        # Test values where the length of the conformation lists are equal
                        elif len(current_conf) == len(test_conf):
                            
                                # Inset this conformation if any of it's values are smaller than the test position's
                                if any([current_conf[i] < test_conf[i] for i in range(len(current_conf))]):
                                    group_tuple_list.insert(test_position, current_tuple)
                                    break
                        
                    # Add to the end if no insertion conditions are found
                    else:
                         group_tuple_list.append(current_tuple)   
                                
                # Add each sorted group to the growing list
                sorted_groups.extend(group_tuple_list)
            
            # Step 2 - Place linked proteins next to each other, then drugs linked to proteins
            
            # Find simple links in the internal link list - any links to multiple components must (should?) also have a corrosponding simple link tuple
            simple_links = [link for link in self.internal_links if isinstance(link[0], int) and isinstance(link[1], int)]
            
            # Find links that are from proteins to proteins
            protein_links = [link for link in simple_links if isinstance(components[link[0]], Protein) and isinstance(components[link[1]], Protein)]
            
            # Find links that are from drugs to proteins. If needed, flip link to make sure the drug index is listed first.
            drug_links_p1 = [(x, y) for (x, y) in simple_links if isinstance(components[x], Drug) and isinstance(components[y], Protein)]
            drug_links_p2 = [(y, x) for (x, y) in simple_links if isinstance(components[x], Protein) and isinstance(components[y], Drug)]    
            
            # Make a list of links to use for moving, in the order desired
            links_for_moving = protein_links + drug_links_p1 + drug_links_p2
            
            # Move stuff around
            for current_link in links_for_moving:# Local function to move sorting tuple
           
                # Find first element and remove it
                [pop_index] = [i for i, x in enumerate(sorted_groups) if x[1] == current_link[0]]
                moving_tuple = sorted_groups.pop(pop_index)
                
                # Find the second element and insert at that position
                [insert_index] = [i for i, x in enumerate(sorted_groups) if x[1] == current_link[1]]
                sorted_groups.insert(insert_index, moving_tuple)
                
            # Step 3 - Convert to name string
            
            # Process each sorting tuple in order
            symbol_list = []
            for current_tuple in sorted_groups:
                current_component = current_tuple[0]
                
                # Get the symbol for the component
                comp_str = current_component.symbol
                
                # Get the symbol(s) for the conformation(s)
                if isinstance(current_component, Protein):
                    conf_str_list = [current_component.conformation_symbols[i] for i in conformations[current_tuple[1]]]
                    conf_str = ','.join(conf_str_list)
                else:
                    conf_str = ''
                
                # Stick the parts together and save in list
                symbol_list.append(comp_str + conf_str)
                
            # Assemble and assign to the state
            self.symbol = ''.join(symbol_list)
        # End Method 0
    
    def autoname(self):
        # Creates an automated name for state
        # Run after numbering
               
        self.name = 'State {}'.format(self.number)
    
    def autovariable(self):
        # Creates a sympy symbol object that represents the state's concentration in rate equations. 
        # Call after numbering by the network function.
        # NOTE: we call the sympy symbol a "variable", and the single-letter notation for chemical components a
        #   "symbol", in line with usage in basic algebra and chemistry language. Sorry, symbolic mathmatics.
        
        self.variable = sp.symbols('S_{}'.format(self.number))
            
    def generate_component_list(self, return_components_only = False):
        # Returns a list of components and a list of conformations
        # The list is ordered all drugs first, then proteins
        
        # Concatenate lists of components and mark the drug conformations as None 
        component_sorted = [*self.required_drug_list, *self.required_protein_list]
        if not return_components_only:
            conformation_sorted = [*itertools.repeat(None, len(self.required_drug_list)), *self.req_protein_conf_lists]
            return component_sorted, conformation_sorted
        elif return_components_only:
            return component_sorted
        
    def add_component_list(self, incoming_components, incoming_conformations):
        # Adds component and conformation lists to the state
         
        # Seperate into lists of drugs and proteins
        self.required_drug_list.extend([x for x in incoming_components if isinstance(x, Drug)])
        self.required_protein_list.extend([x for x in incoming_components if isinstance(x, Protein)])
        self.req_protein_conf_lists.extend([incoming_conformations[i] for i, x in enumerate(incoming_components) if isinstance(x, Protein)])
    
    def enumerate_components(self, return_conformations = False):
        # A generator that returns an integer index and the corrsoponding component (and conformation if requested) of all the state's components
        # The components are ordered all drugs first, then proteins
        
        # Loop through the drugs first
        comp_num = 0
        for drug in self.required_drug_list:
            if return_conformations:
                yield (comp_num, self.required_drug_list[comp_num], None)
            else:
                yield (comp_num, self.required_drug_list[comp_num])
            comp_num += 1
        
        # Now we loop through proteins after saving the number of drugs
        offset_num = comp_num # Note: integers are immutable
        for protein in self.required_protein_list:
            if return_conformations:
                yield (comp_num, self.required_protein_list[comp_num - offset_num], self.req_protein_conf_lists[comp_num - offset_num])
            else:
                yield (comp_num, self.required_protein_list[comp_num - offset_num])
            comp_num += 1
    
    def get_component_by_number(self, request_number, return_conformations = False):
        # Return a state component (and conformations if requested) by number
        
        # Need to deal with two lists of components
        offset_num = len(self.required_drug_list)
        max_num = offset_num + len(self.required_protein_list)
        
        # Not going to implement negative indexing here
        if request_number < 0:
            raise IndexError('index out of bounds, negative indexing not supported')
        
        # If number is less than the offset, return a drug
        elif request_number < offset_num:
            if return_conformations:
                return self.required_drug_list[request_number], None
            else:
                return self.required_drug_list[request_number]
        
        # If number is equal or more than offset but still valid, return a protein
        elif request_number < max_num:
            if return_conformations:
                return self.required_protein_list[request_number - offset_num], self.req_protein_conf_lists[request_number - offset_num]
            else:
                return self.required_protein_list[request_number - offset_num]
        
        # Can't ask for more components than are listed in the state
        elif request_number >= max_num:
            raise IndexError('index out of bounds, value too high')
        
# Define classes for the different types of state transitions - edges on network graph
class StateTransition(HasTraits):    
    #Superclass for all transition objects
    
    # Traits initialization
    name = Str()
    number = Int(None)
    variable = Instance(spBaseClass) # Must be a sympy object
    ID = Instance(uuid.UUID)
    
     # Want to give a new state an ID right away, determined by the network generation code
    def __init__(self, *args, **kwargs):
        self.ID = uuid.uuid4()
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery
        
    def autovariable(self):
        # Creates a sympy symbol object that represents the edges's rate constant in rate equations. 
        # Call after numbering by the network function.
        # NOTE: we call the sympy symbol a "variable", and the single-letter notation for chemical components a
        #   "symbol", in line with usage in basic algebra and chemistry language. Sorry, symbolic mathmatics.
        
        self.variable = sp.symbols('k_{}'.format(self.number))
    
class Conversion(StateTransition):
    
    # Traits initialization
    reference_direction = Bool()
    
    # Keep track of the original conversion direction for better numbering
    def __init__(self, reference_direction = True, *args, **kwargs):
        self.reference_direction = reference_direction
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery 
    
class Association(StateTransition):
    pass
    
class Dissociation(StateTransition):
    pass

class RE_Conversion(StateTransition):
    pass
   
class RE_Dissociation(StateTransition):
    pass

class RE_Association(StateTransition):
    pass
    
# Define class for the rules used to build the network graph
class Rule(HasTraits):
    # Rules must be created with a reference to the model they belong to, and 
    # and must be re-created if the components in the model change. 
    
    # To make sense of the variable names, think of rules as simple sentences 
    # with a grammatical subject and and object. It is unfortunate that 
    # object-oriented programming also shares terminology with grammar.
    
    # Traits initialization
    subject_conf = List(Either(List(Int), None))
    object_conf = List(Either(List(Int), None))
    # Note: rule_subject and rule_object defined dynamically in __init__ as Enum(all listed Drug and Protein objects in the Model)
    
    # Possible rules
    _rule_choices = [' associates with ',
                     ' dissociates from ',
                     ' reversibly associates with ',
                     ' reversibly dissociates from ',
                     ' associates and dissociates in rapid equlibrium with ',
                     ' dissociates and reassociates in rapid equlibrium from ',
                     ' converts to ',
                     ' reversibly converts to ',
                     ' converts in rapid equlibrium to ',
                     ' is competitive with ']
                     #' can only bind along with ']
                     #' is constrained to be the same value as ',
                     #' does not exist.'] #Not sure how to implement these rules right now, maybe need a refactor into different types of rules? 
    rule = Enum(*_rule_choices)
    
    # Create a list of possible component choices before completeing Traits initalization
    def __init__(self, model, *args, **kwargs):
        self._possible_components = [*model.drug_list, *model.protein_list]
        self.add_trait('rule_subject', List(Enum(*self._possible_components)))
        self.add_trait('rule_object', List(Enum(*self._possible_components))) 
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery
    
    # Method to check if the rule is valid
    def check_rule_traits(self):
        # Checks only for generally valid configurations of components and associated configurations. 
        # Specific requirements of each individual rule type are checked in the network generation code on the Model class 
        
        # Check each subject's traits
        for current_subject, current_subject_conf in zip(self.rule_subject, self.subject_conf):

            # Drugs do not have conformations
            if isinstance(current_subject, Drug):
                if current_subject_conf != None:
                    raise RuleNotValidError('Drugs cannot have listed conformations')
                        
            # Proteins cannot choose more conformations than are available and may not have a None value
            if isinstance(current_subject, Protein):
                if current_subject_conf == None:
                    raise RuleNotValidError('Proteins must have at least one conformation participating in rule')
                if current_subject_conf != None and any(index >= len(current_subject.conformation_names) for index in current_subject_conf):
                    raise RuleNotValidError('Rule cannot be given a conformation index value corrosponding to more than the number of conformations available to the protein')
                    
        # Check each object's traits. It occurs to me just now that this could be refactored to avoid a bit of code repeat. Low priorty.
        for current_object, current_object_conf in zip(self.rule_object, self.object_conf):
            
            # Drugs do not have conformations
            if isinstance(current_object, Drug):
                if current_object_conf != None:
                    raise RuleNotValidError('Drugs cannot have listed conformations')
                    
            # Proteins cannot choose more conformations than are available and may not have a None value          
            if isinstance(current_object, Protein):
                if current_object_conf == None:
                    raise RuleNotValidError('Proteins must have at least one conformation participating in rule')
                if current_object_conf != None and any(index >= len(current_object.conformation_names) for index in current_object_conf):
                    raise RuleNotValidError('Rule cannot be given a conformation index value corrosponding to more than the number of conformations available to the protein')
                    
        if self.rule == ' dissociates from ' or self.rule == ' reversibly dissociates from ' \
                or self.rule == ' dissociates and reassociates in rapid equlibrium from ':
            
            check_signature = CountingSignature('conformations included')
            check_signature.count_for_subject(self.rule_subject, self.subject_conf)
            check_signature.count_for_object(self.rule_object, self.object_conf)
            check_signature.third_state_count = check_signature.object_count - check_signature.subject_count
            if not check_signature.object_count == check_signature.third_state_count + check_signature.subject_count:
                raise RuleNotValidError('Dissociation rules require that the components in the subject state be fully included in the object state')
    
    # Method returns the components involved in the rule as a standardized list
    def generate_component_list(self, what_to_include = 'both'):
        # Returns a list of components and a list of conformations
        # The list is ordered all drugs first, then proteins with specific conformations, then proteins with any conformation allowed
        # what_to_include = 'subject', 'object', or 'both'
        
        # Concatenate lists of components
        if what_to_include == 'both':
            component_list = [*self.rule_subject, *self.rule_object]
            conformation_list = [*self.subject_conf, *self.object_conf]
        elif what_to_include == 'subject':
            component_list = self.rule_subject
            conformation_list = self.subject_conf
        elif what_to_include == 'object':
            component_list = self.rule_object
            conformation_list = self.object_conf
        elif what_to_include == 'difference':
            all_tuples = zip(self.rule_object, self.object_conf)
            sub_tuples = zip(self.rule_subject, self.subject_conf)
            component_list, conformation_list  = zip(*[x for x in all_tuples if x not in sub_tuples])
        else:
            raise ValueError("Function supplied with incorrect option for parameter 'what_to_include'")
        
        # Sort the list
        # Seperate into partial lists of drugs and proteins
        drug_comp_part = [x for x in component_list if isinstance(x, Drug)]
        drug_conf_part = [conformation_list[i] for i, x in enumerate(component_list) if isinstance(x, Drug)]
        protein_comp_part = [x for x in component_list if isinstance(x, Protein)]
        protein_conf_part = [conformation_list[i] for i, x in enumerate(component_list) if isinstance(x, Protein)]
        
        # Put the proteins with a [] conformation at the back of the list
        put_in_back = [i for i, x in enumerate(protein_conf_part) if x == []]
        protein_comp_backpart = []
        protein_conf_backpart = []
        put_in_back.reverse() # If we loop through the list starting from the back, we avoid re-indexing problems with pop()
        for i in put_in_back:
            protein_comp_backpart.append(protein_comp_part.pop(i))
            protein_conf_backpart.append(protein_conf_part.pop(i))
        
        # Add all parts back together and return
        component_sorted = [*drug_comp_part, *protein_comp_part, *protein_comp_backpart]
        conformation_sorted = [*drug_conf_part, *protein_conf_part, *protein_conf_backpart]
        return component_sorted, conformation_sorted
    
    # Method returns tuples of signatures with all possibilites of the appropreate components/conformations
    def generate_signature_list(self): 
        # Returns a list of all possible valid signatures that can be created with the rule
        # Automatically returns component only signatures when no part of the rule askes for it
       
        # Make new list
        signature_list = []
        
        # Association rules
        if self.rule == ' associates with ' or self.rule == ' reversibly associates with ' \
                    or self.rule == ' associates and dissociates in rapid equlibrium with ':
            
            # Different rules need different output tuples - Perhaps there is a chance to combine this code with some in the model, but I currently don't know how to refactor it better.
            # If we have all None or [] conformations, this is pretty simple and we return a component only signature
            if all([x == None or x == [] for x in self.subject_conf]) and all([x == None or x == [] for x in self.object_conf]):
                new_sig = CountingSignature('components only')
                new_sig.count_for_subject(*self.generate_component_list('subject'))
                new_sig.count_for_object(*self.generate_component_list('object'))
                new_sig.count_for_third_state(*self.generate_component_list('both'))
                signature_list.append(new_sig)
            
            # If we have to worry about about conformations, but there is no "any" conformations ([]), this is also pretty simple
            elif all([x != [] for x in self.subject_conf]) and all([x != [] for x in self.object_conf]):
                new_sig = CountingSignature('conformations included')
                new_sig.count_for_subject(*self.generate_component_list('subject'))
                new_sig.count_for_object(*self.generate_component_list('object'))
                new_sig.count_for_third_state(*self.generate_component_list('both'))
                signature_list.append(new_sig)
           
            # We must have at least one "any" conformation but not all, so we make a list of all possible conformations that could fit the rule
            else:
                # Save the lists of components/conformations
                sub_comp, sub_conf = self.generate_component_list('subject')
                obj_comp, obj_conf = self.generate_component_list('object')
                third_state_comp = sub_comp + obj_comp
                
                # Get all the possible conformation combinations that can be made
                new_sub_conf_list, new_obj_conf_list = self._get_all_conformation_combinations(sub_conf, obj_conf, sub_comp, obj_comp)
                
                # Now create signatures from the lists                   
                for new_sub_conf, new_obj_conf in zip(new_sub_conf_list, new_obj_conf_list):
                    new_sig = CountingSignature('conformations included')
                    new_sig.count_for_subject(sub_comp, new_sub_conf)
                    new_sig.count_for_object(obj_comp, new_obj_conf)
                    combinded_third_state_conf = new_sub_conf + new_obj_conf
                    new_sig.count_for_third_state(third_state_comp, combinded_third_state_conf)
                    signature_list.append(new_sig)
        
        # Dissociation rules            
        elif self.rule == ' dissociates from ' or self.rule == ' reversibly dissociates from ' \
                or self.rule == ' dissociates and reassociates in rapid equlibrium from ':
                
            # If we have all None or [] conformations, this is pretty simple and we return a component only signature
            if all([x == None or x == [] for x in self.object_conf]):
                new_sig = CountingSignature('components only')
                new_sig.count_for_subject(*self.generate_component_list('subject'))
                new_sig.count_for_object(*self.generate_component_list('object'))
                new_sig.third_state_count = new_sig.object_count - new_sig.subject_count
                signature_list.append(new_sig)
            
            # If we have to worry about about conformations, but there is no "any" conformations ([]), this is also pretty simple
            elif all([x != [] for x in self.subject_conf]) and all([x != [] for x in self.object_conf]):
                new_sig = CountingSignature('conformations included')
                new_sig.count_for_subject(*self.generate_component_list('subject'))
                new_sig.count_for_object(*self.generate_component_list('object'))
                new_sig.third_state_count = new_sig.object_count - new_sig.subject_count
                signature_list.append(new_sig)
           
            # We must have at least one "any" conformation but not all, so we make a list of all possible conformations that could fit the rule
            else:
                # Save the lists of components/conformations
                sub_comp, sub_conf = self.generate_component_list('subject')
                obj_comp, obj_conf = self.generate_component_list('object')
                # Get all the possible conformation combinations that can be made
                new_sub_conf_list, new_obj_conf_list = self._get_all_conformation_combinations(sub_conf, obj_conf, sub_comp, obj_comp)
                
                # Now create signatures from the lists                   
                for new_sub_conf, new_obj_conf in zip(new_sub_conf_list, new_obj_conf_list):
                    new_sig = CountingSignature('conformations included')
                    new_sig.count_for_subject(sub_comp, new_sub_conf)
                    new_sig.count_for_object(obj_comp, new_obj_conf)
                    new_sig.third_state_count = new_sig.object_count - new_sig.subject_count # On a counter object, subtracting something that doesn't exist is silently allowed and does not create "negative" objects
                    if new_sig.object_count == new_sig.subject_count + new_sig.third_state_count: # We filter out mismatched subject and object conformations here
                        signature_list.append(new_sig)
        
        # Conversion rules
        elif self.rule == ' converts to ' or self.rule == ' reversibly converts to ' \
                or self.rule == ' converts in rapid equlibrium to ':
            
            # If we have all None or [] conformations, this is pretty simple and we return a component only signature
            if all([x == None or x == [] for x in self.subject_conf]) and all([x == None or x == [] for x in self.object_conf]):
                new_sig = CountingSignature('components only')
                new_sig.count_for_subject(*self.generate_component_list('subject'))
                new_sig.count_for_object(*self.generate_component_list('object'))
                signature_list.append(new_sig)
            
            # If we have to worry about about conformations, but there is no "any" conformations ([]), this is also pretty simple
            elif all([x != [] for x in self.subject_conf]) and all([x != [] for x in self.object_conf]):
                new_sig = CountingSignature('conformations included')
                new_sig.count_for_subject(*self.generate_component_list('subject'))
                new_sig.count_for_object(*self.generate_component_list('object'))
                signature_list.append(new_sig)
                
            # We must have at least one "any" conformation but not all, so we make a list of all possible conformations that could fit the rule
            else:
                # Save the lists of components/conformations
                sub_comp, sub_conf = self.generate_component_list('subject')
                obj_comp, obj_conf = self.generate_component_list('object')
                # Get all the possible conformation combinations that can be made
                new_sub_conf_list, new_obj_conf_list = self._get_all_conformation_combinations(sub_conf, obj_conf, sub_comp, obj_comp)
                
                # Now create signatures from the lists                   
                for new_sub_conf, new_obj_conf in zip(new_sub_conf_list, new_obj_conf_list):
                    new_sig = CountingSignature('conformations included')
                    new_sig.count_for_subject(sub_comp, new_sub_conf)
                    new_sig.count_for_object(obj_comp, new_obj_conf)
                    
                    # We don't want to add a signature that doesn't change anything between the states
                    if new_sig.subject_count - new_sig.object_count != Counter() or new_sig.object_count - new_sig.subject_count != Counter(): # Could have dissapearing states, check subtraction in both directions
                        signature_list.append(new_sig)
        
        # Competition rule
        elif self.rule == ' is competitive with ':   
            
            # If we have all None or [] conformations, this is pretty simple and we return a component only signature
            if all([x == None or x == [] for x in self.subject_conf]) and all([x == None or x == [] for x in self.object_conf]):
                new_sig = CountingSignature('components only')
                new_sig.count_for_subject(*self.generate_component_list('subject'))
                new_sig.count_for_object(*self.generate_component_list('object'))
                signature_list.append(new_sig)
            
            # If we have to worry about about conformations, but there is no "any" conformations ([]), this is also pretty simple
            elif all([x != [] for x in self.subject_conf]) and all([x != [] for x in self.object_conf]):
                new_sig = CountingSignature('conformations included')
                new_sig.count_for_subject(*self.generate_component_list('subject'))
                new_sig.count_for_object(*self.generate_component_list('object'))
                signature_list.append(new_sig)
           
            # We must have at least one "any" conformation but not all, so we make a list of all possible conformations that could fit the rule
            else:
                # Save the lists of components/conformations
                sub_comp, sub_conf = self.generate_component_list('subject')
                obj_comp, obj_conf = self.generate_component_list('object')
                
                # Get all the possible conformation combinations that can be made
                new_sub_conf_list, new_obj_conf_list = self._get_all_conformation_combinations(sub_conf, obj_conf, sub_comp, obj_comp)
                
                # Now create signatures from the lists                   
                for new_sub_conf, new_obj_conf in zip(new_sub_conf_list, new_obj_conf_list):
                    new_sig = CountingSignature('conformations included')
                    new_sig.count_for_subject(sub_comp, new_sub_conf)
                    new_sig.count_for_object(obj_comp, new_obj_conf)
                    signature_list.append(new_sig)
        
        else:
            raise ValueError("Rule type not recognized")

        return signature_list
    
    def _get_all_conformation_combinations(self, sub_conf, obj_conf, sub_comp, obj_comp):
    # Helper method to find all possible combinations of conformations when there is at least one "any" conformation given 
    
    # NOTE: Changed to only create single conformations, but not sure if this is what we want in the end. Left multiple conformation code as comments.    
    
    # Find the "any" conformations
        sub_any_index = []
        for index, current_conf in enumerate(sub_conf):
            if current_conf == []:
                sub_any_index.append(index)
        obj_any_index = []
        for index, current_conf in enumerate(obj_conf):
            if current_conf == []:
                obj_any_index.append(index)                
        
        # Get the number of allowed conformations for each "any"
        allowed_conformations = []
        for sub_index in sub_any_index:
            allowed_conformations.append(range(len(sub_comp[sub_index].conformation_names)))
        for obj_index in obj_any_index:
            allowed_conformations.append(range(len(obj_comp[obj_index].conformation_names)))
        
        # Get all combinations of each allowed list
        allowed_conf_combos = []
        for allowed_conf in allowed_conformations: # Do this for each list of "any" conformations
            new_combos = []
            new_combos.extend(itertools.combinations(allowed_conf, 1))
#            for num_conf in range(1, len(allowed_conf) + 1): # Need to repeat for single, double, etc. allowed number of conformations
#                new_combos.extend(itertools.combinations(allowed_conf, num_conf))
            allowed_conf_combos.append(new_combos) # Add the list of combos to the main list
        
        # Now we construct a list of conformations lists for each combination
        new_sub_conf_list = []
        new_obj_conf_list = []
        for comboset in itertools.product(*allowed_conf_combos):
            
            # We need to replace the "any" conformations in the conformation lists
            new_sub_conf = sub_conf.copy()
            for combo_index, sub_index in enumerate(sub_any_index):
                new_sub_conf[sub_index] = comboset[combo_index]
            new_obj_conf = obj_conf.copy()
            for combo_index, obj_index in enumerate(obj_any_index):
                new_obj_conf[obj_index] = comboset[combo_index + len(sub_any_index)]
                
            # Add to lists
            new_sub_conf_list.append(new_sub_conf)
            new_obj_conf_list.append(new_obj_conf)
                
        return new_sub_conf_list, new_obj_conf_list
    
# Define class as a container for the network graphs of states
# Do we really want a seperate container that's just added onto a Model class? Don't know yet
class Network(HasTraits):
    # Biochemical networks are stored and manipulated using NetworkX graphs.
    
    # Traits initialization
    main_graph = Instance(nx.DiGraph)
    main_graph_blacklist = List(Instance(State))
    display_graphs = List(Instance(nx.DiGraph))
    solving_graphs = List(Instance(nx.DiGraph))
    
    # Initial network is an empty main_graph and empty lists of derivitive graphs
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery
        self.main_graph = nx.DiGraph()
        self.display_graphs = []
        self.solving_graphs = []
        
    def autosymbol(self):
        # Give each state in the main graph a symbol
        
        # Call autosymbol on each state
        for current_state in self.main_graph:
            current_state.autosymbol()
         
    def autonumber(self):
        # Give each state and edge state-transition object in the main graph a number
        # Sort states/edges by length for numbering, but no more organization than that
    
        # Get the states as a list sorted by length of component lists
        state_list = sorted(self.main_graph.__iter__(), key = lambda x: len(x.required_drug_list) + len(x.required_protein_list))
       
        # Number states
        for current_index, current_state in enumerate(state_list):
            current_state.number = current_index + 1 # Avoid using number 0, since a null state often has a specific meaning
       
        # Number edge ST objects via their connection to the states (so in general, low number edges on low number states)
        for current_state in state_list:
            
            # Find edges with current state as tail and visit each one
            edge_tuples = [x for x in self.main_graph.out_edges(current_state, 'reaction_type')]
            for current_edge_tuple in edge_tuples:
                STobj = current_edge_tuple[2]
                Reverse_STobj = None # Reset reverse
                
                # If edge is already numbered, go to the next edge
                if STobj.number != None:
                    continue
                
                # Otherwise, we need to give a new number out
                new_number = self._get_next_edge_number()

                # Look for opposite rule
                try:
                    Reverse_STobj = self.main_graph.edges[current_edge_tuple[1], current_edge_tuple[0]]['reaction_type']
                except KeyError: # This is OK, there's just not an edge there
                    assign_reverse = False 
                else:
                    if Reverse_STobj.number == None: # Only assign a new value if it's not already numbered
                        assign_reverse = True
                    else:
                        assign_reverse = False
  
                # Handle association rules
                if isinstance(STobj, Association):
                    
                    # Give new number as positive value
                    STobj.number = new_number
                    
                    # Give the opposite reaction the negative number, if needed
                    if assign_reverse:
                        Reverse_STobj.number = -new_number 
                
                # Handle dissociation rules
                elif isinstance(STobj, Dissociation):
                    
                    # Give new number as negative value
                    STobj.number = -new_number
                    
                    # Give the opposite reaction the positive number, if needed
                    if assign_reverse:
                        Reverse_STobj.number = new_number # Give the opposite reaction the negative number, if needed
                
                # Handle Conversion rules
                elif isinstance(STobj, Conversion):
                    
                    # If this is the reference direction, positive numbers go forward
                    if STobj.reference_direction:
                        # Give new number as positive value
                        STobj.number = new_number
                    
                        # Give the opposite reaction the negative number, if needed
                        if assign_reverse:
                            Reverse_STobj.number = -new_number
                    
                    # If this is not the reference direction, negative numbers go forward
                    else: 
                        # Give new number as negative value
                        STobj.number = -new_number
                    
                        # Give the opposite reaction the positive number, if needed
                        if assign_reverse:
                            Reverse_STobj.number = new_number
    
    def autovariable(self):
        # Give each state and edge state-transition object in the main graph a variable
        # Call after numbering
        
        # Call autovariable on each state
        for current_state in self.main_graph:
            current_state.autovariable()
            
        # Call autovairable on each edge 
        for u, v, STobj in self.main_graph.edges.data('reaction_type'):
            STobj.autovariable() 
    
    def autoname(self):
        # Give each state in the main graph a variable
        # Call after numbering
        
        # Call autoname on each state
        for current_state in self.main_graph:
            current_state.autoname()
    
    def _get_next_edge_number(self, graph = None):
        # Helper function to see what the next availiable number in the graph is
        
        # Work on main graph by default, but any graph can be given
        if graph == None: 
            graph = self.main_graph 
            
        # Make a list of the numbers already used - negative numbers for opposite reactions are converted to their positive counterpart
        used_numbers = [abs(STobj.number) for (u, v, STobj) in graph.edges.data('reaction_type') if STobj.number != None]
        
        # Test for presence in the number list until an open value is found
        test_number = 1 # Don't use 0, as null states have a specific meaning in chemical kinetics
        while test_number in used_numbers:
            test_number += 1
        
        # Give the number back, positive values only
        return test_number


class CountingSignature(HasTraits):
    # Used for storing information about the number of elements involved in state to state transistion reactions
    # Two use modes with 1st parameter "count_type":
        # 'components only' - Don't look at conformations of proteins at all
        # 'conformations included' - Count using (component, conformations) tuples
    
    # Traits initialization
    count_type = Either('components only', 'conformations included')
    subject_count = Instance(Counter)
    object_count = Instance(Counter)
    third_state_count = Instance(Counter)
    
    # Can be created with any number of missing state or rule information
    def __init__(self, count_type, subject_state = None, object_state = None, third_state = None, *args, **kwargs):
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery
        self.count_type = count_type
        if subject_state != None:
            self.subject_count = self._count_state(subject_state)
        if object_state != None:
            self.object_count = self._count_state(object_state)
        if third_state != None:
            self.third_state_count = self._count_state(third_state)
      
    def count_for_subject(self, components, conformations = None):
    # Directly count a list of components/conformations into the subject part of the signature
        self.subject_count = self._count_list(components, conformations)
    
    def count_for_object(self, components, conformations = None):
    # Directly count a list of components/conformations into the object part of the signature
        self.object_count = self._count_list(components, conformations)
    
    def count_for_third_state(self, components, conformations = None):
    # Directly count a list of components/conformations into the third state part of the signature   
        self.third_state_count = self._count_list(components, conformations)
    
    # Count the component/conformations given by a particular state
    def _count_state(self, state):
        components, conformations = state.generate_component_list()
        return self._count_list(components, conformations)
    
    # Send the hashable lists to collections.Counter and send the counting dictionary back to the signature
    def _count_list(self, components, conformations):
        # Use the counter directly on the list of components 
        if self.count_type == 'components only':
            return Counter(components)
       
        # Use the counter on a zipped list of component, conformations) tuples 
        elif self.count_type == 'conformations included':
            # Need an actual list of conformations for this type of signature
            if conformations == None: #Default value passed from a calling function
                raise ValueError("Function requries a list of conformations")
                
            # Convert the conformation lists to hashable tuples
            hashable_conformations = [(*x,) if x != None else x for x in conformations]
            return Counter([*zip(components, hashable_conformations)])
        
    