"""Classes representing physical or conceptual objects needed to describe the
biochemical system, the model, and the experiments performed. 

"""

import uuid
import itertools
import networkx as nx
from traits.api import HasTraits, Str, List, Int, Instance, Enum, Either
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
    
    # Want to give a new state an ID right away, determined by the network generation code
    def __init__(self, *args, **kwargs):
        self.ID = uuid.uuid4()
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery 
    
    # Names and symbols are created from the properties of the state, but the state's number is generated at a model level or carried over from parent models
    def autosymbol(self):
        pass
    def autoname(self):
        pass
    
    def generate_component_list(self):
        # Returns a list of components and a list of conformations
        # The list is ordered all drugs first, then proteins
        
        # Concatenate lists of components and mark the drug conformations as None 
        component_sorted = [*self.required_drug_list, *self.required_protein_list]
        conformation_sorted = [*itertools.repeat(None, len(self.required_drug_list)), *self.req_protein_conf_lists]
        return component_sorted, conformation_sorted
     
    def add_component_list(self, incoming_components, incoming_conformations):
        # Adds component and conformation lists to the state
         
        # Seperate into lists of drugs and proteins
        self.required_drug_list.extend([x for x in incoming_components if isinstance(x, Drug)])
        self.required_protein_list.extend([x for x in incoming_components if isinstance(x, Protein)])
        self.req_protein_conf_lists.extend([incoming_conformations[i] for i, x in enumerate(incoming_components) if isinstance(x, Protein)])
    
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
                     ' associates and dissociates in rapid equlibrium with ',
                     ' converts to state ',
                     ' reversibly converts to state ',
                     ' converts in rapid equlibrium to state ',
                     ' does not exist in the same complex as ',
                     ' can only be in complex along with ']
                     #' is constrained to be the same value as ',
                     #' does not exist.'] #Not sure how to implement these rules right now, maybe need a refactor into different types of rules? 
    rule = Enum(*_rule_choices)
    
    # Create a list of possible component choices before completeing Traits initalization
    def __init__(self, model, *args, **kwargs):
        self._components = [*model.drug_list, *model.protein_list]
        self.add_trait('rule_subject', List(Enum(*self._components)))
        self.add_trait('rule_object', List(Enum(*self._components))) 
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
        
    