"""Classes representing physical or conceptual objects needed to describe the
biochemical system, the model, and the experiments performed. 

"""

import uuid
import itertools
import networkx as nx
from collections import Counter
from traits.api import HasTraits, Str, List, Tuple, Int, Instance, Enum, Either
from bikipy.bikicore.exceptions import ComponentNotValidError, RuleNotValidError

# Define classes for the different components of the biochemical system
class Drug(HasTraits):
    
    # Traits initialization
    name = Str()
    symbol = Str()
    ID = Instance(uuid.UUID) 
    
class Protein(HasTraits):
    
    # Traits initialization
    name = Str()
    symbol = Str()
    conformation_names = List(Str()) 
    conformation_symbols = List(Str())
    ID = Instance(uuid.UUID) 
    
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
    ID = Instance(uuid.UUID)
    required_drug_list = List(Instance(Drug))
    required_protein_list = List(Instance(Protein))
    req_protein_conf_lists = List(List(Int()))
    internal_links = List(Tuple())
    
    # Want to give a new state an ID right away, determined by the network generation code
    def __init__(self, *args, **kwargs):
        self.ID = uuid.uuid4()
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery 
    
    # Names and symbols are created from the properties of the state, but the state's number is generated at a model level or carried over from parent models
    def autosymbol(self):
        pass
    def autoname(self):
        pass
    
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
    number = Int()
    ID = Instance(uuid.UUID)
    
class ConformationalChange(StateTransition):
    
    # Traits initialization
    name = Str('activation')
    
class Association(StateTransition):
    
    # Traits initialization
    name = Str('association')

class Dissociation(StateTransition):
    
    # Traits initialization
    name = Str('dissociation')

### Leave these for now, but they shouold be simple enough to make when the basic structure of the library is finished
class RE_ConformationalChange(StateTransition):
    
    # Traits initialization
    name = Str('activation')    

class RE_Dissociation(StateTransition):
    
    # Traits initialization
    name = Str('dissociation')
    
class RE_Association(StateTransition):
    
    # Traits initialization
    name = Str('association')
    
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
                     ' does not exist in the same complex as ',
                     ' can only be in complex along with ']
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
        elif self.rule == ' converts to ':
            
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
    display_graphs = List(Instance(nx.DiGraph))
    solving_graphs = List(Instance(nx.DiGraph))
    
    # Initial network is an empty main_graph and empty lists of derivitive graphs
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery
        self.main_graph = nx.DiGraph()
        self.display_graphs = []
        self.solving_graphs = []
        
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
        
        
    