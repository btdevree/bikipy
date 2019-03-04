"""Test suite for classes in components.py

"""
import pytest
import bikipy.bikicore.components as bkcc
import bikipy.bikicore.model as bkcm
from bikipy.bikicore.exceptions import ComponentNotValidError, RuleNotValidError

#---- Testing fixtures ----

# Create a default Drug object for reuse in tests
@pytest.fixture()
def default_Drug_instance():
    ddi = bkcc.Drug()
    ddi.name = 'adrenaline'
    ddi.symbol = 'A'
    return ddi
    
# Create a default Protein object for reuse in tests
@pytest.fixture()
def default_Protein_instance():
    dpi = bkcc.Protein()
    dpi.name = 'beta adrenergic receptor'
    dpi.symbol = 'R'
    dpi.conformation_names = ['inactive', 'active']
    dpi.conformation_symbols = ['', '*']
    return dpi

# Create a default ConformationalChange object for reuse in tests
@pytest.fixture()
def default_ConfChange_instance():
    return bkcc.ConformationalChange()
    
# Create a default Association object for reuse in tests
@pytest.fixture()
def default_Association_instance():
    return bkcc.Association()

# Create a default Dissociation object for reuse in tests
@pytest.fixture()
def default_Dissociation_instance():
    return bkcc.Dissociation()
    
# Create a default RE_ConformationalChange object for reuse in tests
@pytest.fixture()
def default_RE_ConfChange_instance():
    return bkcc.RE_ConformationalChange()
    
# Create a default RE_Association object for reuse in tests
@pytest.fixture()
def default_RE_Association_instance():
    return bkcc.RE_Association()

# Create a default RE_Dissociation object for reuse in tests
@pytest.fixture()
def default_RE_Dissociation_instance():
    return bkcc.RE_Dissociation()

# Create a default Model object from bkcm for reuse in tests
@pytest.fixture()
def default_Model_instance(default_Drug_instance, default_Protein_instance):
    newmodel = bkcm.Model(1, 'default model', None)
    newmodel.drug_list.append(default_Drug_instance)
    newmodel.protein_list.append(default_Protein_instance)
    return newmodel

# Create a default Rule object for reuse in tests
@pytest.fixture()
def default_Rule_instance(default_Model_instance):
    return bkcc.Rule(default_Model_instance)
    
# Create a default Network object for reuse in tests
@pytest.fixture()
def default_Network_instance():
    return bkcc.Network()
    
# Create a default State object for reuse in tests
@pytest.fixture()
def default_State_instance():
    return bkcc.State()


# ---- Unit tests ----

# --Tests for Drug objects--

# Test if Drug objects have the required properties
def test_Drug_has_name(default_Drug_instance):
    assert hasattr(default_Drug_instance, 'name')
def test_Drug_has_symbol(default_Drug_instance):
    assert hasattr(default_Drug_instance, 'symbol')
    
# --Tests for Protein objects--

# Test if Protein objects have the required properties
def test_Protein_has_name(default_Protein_instance):
    assert hasattr(default_Protein_instance, 'name')
def test_Protein_has_symbol(default_Protein_instance):
    assert hasattr(default_Protein_instance, 'symbol')
def test_Protein_has_conformation_names(default_Protein_instance):
    assert hasattr(default_Protein_instance, 'conformation_names')
def test_Protein_has_conformation_symbols(default_Protein_instance):
    assert hasattr(default_Protein_instance, 'conformation_symbols')
    
# Test if argument configurations that Traits can't handle throw exceptions correctly
def test_Protein_right_number_conformation(default_Protein_instance):
    default_Protein_instance.check_protein_traits() # No error expected

def test_Protein_wrong_number_conformation_names(default_Protein_instance):
    too_many_names = ['active', 'inactive', 'lazy']
    default_Protein_instance.conformation_names = too_many_names
    with pytest.raises(ComponentNotValidError):
        default_Protein_instance.check_protein_traits() # Error expected

# --Tests for State objects--

# Test if State objects have the required properties
def test_State_has_name(default_State_instance):
    assert hasattr(default_State_instance, 'name')
def test_State_has_number(default_State_instance):
    assert hasattr(default_State_instance, 'number')
def test_State_has_ID(default_State_instance):
    assert hasattr(default_State_instance, 'ID')
def test_State_has_req_drug_list(default_State_instance):
    assert hasattr(default_State_instance, 'required_drug_list')
def test_State_has_req_protein_list(default_State_instance):
    assert hasattr(default_State_instance, 'required_protein_list')
def test_State_has_req_protein_conf_list(default_State_instance):
    assert hasattr(default_State_instance, 'req_protein_conf_lists')


# --Tests for ConformationalChange objects--

# Test if ConformationalChange objects have the required properties
def test_ConfChange_has_name(default_ConfChange_instance):
    assert hasattr(default_ConfChange_instance, 'name')
def test_ConfChange_has_number(default_ConfChange_instance):
    assert hasattr(default_ConfChange_instance, 'number')
def test_ConfChange_has_ID(default_ConfChange_instance):
    assert hasattr(default_ConfChange_instance, 'ID')
    
#Test if Association objects have the required properties
def test_Association_has_name(default_Association_instance):
    assert hasattr(default_Association_instance, 'name')
def test_Association_has_number(default_Association_instance):
    assert hasattr(default_Association_instance, 'number')
def test_Association_has_ID(default_Association_instance):
    assert hasattr(default_Association_instance, 'ID')
    
#Test if Dissociation objects have the required properties
def test_Dissociation_has_name(default_Dissociation_instance):
    assert hasattr(default_Dissociation_instance, 'name')
def test_Dissociation_has_number(default_Dissociation_instance):
    assert hasattr(default_Dissociation_instance, 'number')
def test_Dissociation_has_ID(default_Dissociation_instance):
    assert hasattr(default_Dissociation_instance, 'ID')
    
# Test if RE_ConformationalChange objects have the required properties
def test_RE_ConfChange_has_name(default_RE_ConfChange_instance):
    assert hasattr(default_RE_ConfChange_instance, 'name')
def test_RE_ConfChange_has_number(default_RE_ConfChange_instance):
    assert hasattr(default_RE_ConfChange_instance, 'number')
def test_RE_ConfChange_has_ID(default_RE_ConfChange_instance):
    assert hasattr(default_RE_ConfChange_instance, 'ID')
    
#Test if RE_Association objects have the required properties
def test_RE_Association_has_name(default_RE_Association_instance):
    assert hasattr(default_RE_Association_instance, 'name')
def test_RE_Association_has_number(default_RE_Association_instance):
    assert hasattr(default_RE_Association_instance, 'number')
def test_RE_Association_has_ID(default_RE_Association_instance):
    assert hasattr(default_RE_Association_instance, 'ID')
    
#Test if RE_Dissociation objects have the required properties
def test_RE_Dissociation_has_name(default_RE_Dissociation_instance):
    assert hasattr(default_RE_Dissociation_instance, 'name')
def test_RE_Dissociation_has_number(default_RE_Dissociation_instance):
    assert hasattr(default_RE_Dissociation_instance, 'number')
def test_RE_Dissociation_has_ID(default_RE_Dissociation_instance):
    assert hasattr(default_RE_Dissociation_instance, 'ID')       
    
# --Tests for Rule objects--

#Test if Rule objects have the required properties
def test_Rule_has_subject(default_Rule_instance):
    assert hasattr(default_Rule_instance, 'rule_subject')
def test_Rule_has_subject_conf(default_Rule_instance):
    assert hasattr(default_Rule_instance, 'subject_conf')
def test_Rule_has_object(default_Rule_instance):
    assert hasattr(default_Rule_instance, 'rule_object')
def test_Rule_has_object_conf(default_Rule_instance):
    assert hasattr(default_Rule_instance, 'object_conf')
def test_Rule_has_rule(default_Rule_instance):
    assert hasattr(default_Rule_instance, 'rule')

#Test if the rule raises RuleNotValidError correctly
def test_Rule_check_drug_subject1(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience
    
    # Set up the rule properly and call check_rule_traits
    dri.rule_subject = default_Drug_instance
    dri.subject_conf = None
    dri.check_rule_traits() #No error expected
    
def test_Rule_check_drug_object2(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience
    
    # Set up the rule properly and call check_rule_traits
    dri.rule_object = default_Drug_instance
    dri.object_conf = None
    dri.check_rule_traits() #No error expected

def test_Rule_check_drug_subject3(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect usage, no conformations are allowed for drugs 
    dri.rule_subject = default_Drug_instance
    dri.subject_conf = [0]    
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_drug_subject4(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect usage, None, not an empty list is reqired for drugs
    dri.rule_subject = default_Drug_instance
    dri.subject_conf = []
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_drug_object5(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect usage, no conformations are allowed for drugs 
    dri.rule_object = default_Drug_instance
    dri.object_conf = [0]    
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_drug_object6(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect usage, None, not an empty list is reqired for drugs
    dri.rule_object = default_Drug_instance
    dri.object_conf = []
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected
    
def test_Rule_check_protein_subject7(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect conf list given, None is not allowed for a protein
    dri.rule_subject = default_Protein_instance
    dri.subject_conf = None
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_protein_subject8(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect conf list given, empty list is not allowed
    dri.rule_subject = default_Protein_instance
    dri.subject_conf = []
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_protein_subject9(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect conf list given that contains invalid index numbers
    dri.rule_subject = default_Protein_instance
    dri.subject_conf = [1, 2, 3] # Only 0 and 1 valid for default_Protein_instance
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_protein_object10(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Usual usage of conf list 
    dri.rule_object = default_Protein_instance
    dri.object_conf = [0]
    dri.check_rule_traits() # No error expected

def test_Rule_check_protein_object11(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect conf list given that contains invalid index numbers
    dri.rule_object = default_Protein_instance
    dri.object_conf = [1, 2, 3] # Only 0 and 1 valid for default_Protein_instance
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_protein_object12(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect conf list given, empty list is not allowed
    dri.rule_object = default_Protein_instance
    dri.object_conf = []
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected
        
def test_Rule_check_protein_object13(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect conf list given, None is not allowed for a protein
    dri.rule_object = default_Protein_instance
    dri.object_conf = None
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

# --Tests for Network objects--

#Test if Network objects have the required properties
def test_Network_has_main_graph(default_Network_instance):
    assert hasattr(default_Network_instance, 'main_graph')
     
#continue with writing tests.....