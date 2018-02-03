"""Classes representing physical or conceptual objects needed to describe the
biochemical system, the model, and the experiments performed. 

"""

import uuid
from traits.api import HasTraits, Str, List, Int, Instance, Enum
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
    name = Str('1')
    number = Int(1)
    ID = Instance(uuid.UUID)

# Define classes for the different types of state transitions - edges on network graph
class StateTransition(HasTraits):    
    #Superclass for all transition objects
    
    # Traits initialization
    number = Int(1)
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

# Define classes for the components used to build the model
class Rule(HasTraits):
    # To make sense of the variable names, think of rules as simple sentences 
    #    with a grammatical subject and and object. It is unfortunate that 
    #    object-oriented programming also shares terminology with grammar.
    
    # Traits initialization
    subject_conf = Enum(None, 'all', 'select')
    subject_conf_list = List(Int)
    object_conf = Enum(None, 'all', 'select')
    object_conf_list = List(Int)
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
        
        # Drugs do not have conformations
        if isinstance(self.rule_subject, Drug):
            if self.subject_conf != None and not (self.subject_conf == 'select' and self.subject_conf_list == []):
                raise RuleNotValidError('Drugs cannot have selected conformations')
        if isinstance(self.rule_object, Drug):
            if self.object_conf != None and not (self.object_conf == 'select' and self.object_conf_list == []):
                raise RuleNotValidError('Drugs cannot have selected conformations')
                
        # Proteins must have at least one conformation chosen and cannot choose more conformations than are available
        if isinstance(self.rule_subject, Protein):
            if self.subject_conf == None or (self.subject_conf == 'select' and self.subject_conf_list == []):
                raise RuleNotValidError('Proteins must have at least one conformation participating in rule')
            if self.subject_conf == 'select' and any(index >= len(self.rule_subject.conformation_names) for index in self.subject_conf_list):
                raise RuleNotValidError('Rule cannot be given a conformation index value corrosponding to more than the number of conformations available to the protein')
        if isinstance(self.rule_object, Protein):
            if self.object_conf == None or (self.object_conf == 'select' and self.object_conf_list == []):
                raise RuleNotValidError('Proteins must have at least one conformation participating in rule')
            if self.object_conf == 'select' and any(index >= len(self.rule_object.conformation_names) for index in self.object_conf_list):
                raise RuleNotValidError('Rule cannot be given a conformation index value corrosponding to more than the number of conformations available to the protein')
        