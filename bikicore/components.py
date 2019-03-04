"""Classes representing physical or conceptual objects needed to describe the
biochemical system, the model, and the experiments performed. 

"""

import uuid
import networkx as nx
import numpy as np
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
    subject_conf = Either(List(Int), None)
    object_conf = Either(List(Int), None)
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
        self.add_trait('rule_subject', Enum(*self._components))
        self.add_trait('rule_object', Enum(*self._components)) 
        super().__init__(*args, **kwargs) # Make sure to call the HasTraits initialization machinery
    
    # Method to check if the rule is valid
    def check_rule_traits(self):
        # Checks only for generally valid configurations of components and associated configurations. 
        # Specific requirements of each individual rule type are checked in the network generation code on the Model class 
        
        # Drugs do not have conformations
        if isinstance(self.rule_subject, Drug):
            if self.subject_conf != None:
                raise RuleNotValidError('Drugs cannot have listed conformations')
                # Drugs do not have conformations
        if isinstance(self.rule_object, Drug):
            if self.object_conf != None:
                raise RuleNotValidError('Drugs cannot have listed conformations')
                
        # Proteins must have at least one conformation chosen and cannot choose more conformations than are available
        if isinstance(self.rule_subject, Protein):
            if self.subject_conf == None:
                raise RuleNotValidError('Proteins must have at least one conformation participating in rule')
            if self.subject_conf == []:
                raise RuleNotValidError('Proteins must have at least one conformation participating in rule')
            if self.subject_conf != None and any(index >= len(self.rule_subject.conformation_names) for index in self.subject_conf):
                raise RuleNotValidError('Rule cannot be given a conformation index value corrosponding to more than the number of conformations available to the protein')
        if isinstance(self.rule_object, Protein):
            if self.object_conf == None:
                raise RuleNotValidError('Proteins must have at least one conformation participating in rule')
            if self.object_conf == []:
                raise RuleNotValidError('Proteins must have at least one conformation participating in rule')
            if self.object_conf != None and any(index >= len(self.rule_object.conformation_names) for index in self.object_conf):
                raise RuleNotValidError('Rule cannot be given a conformation index value corrosponding to more than the number of conformations available to the protein')

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
        
    