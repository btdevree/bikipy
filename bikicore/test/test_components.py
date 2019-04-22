"""Test suite for classes in components.py

"""
import pytest
import bikipy.bikicore.components as bkcc
import bikipy.bikicore.model as bkcm
from bikipy.bikicore.exceptions import ComponentNotValidError, RuleNotValidError

#---------------------------- Testing fixtures --------------------------------

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

# Create a default model with configured Rule object for reuse in tests
@pytest.fixture()
def default_Model_irreversable_association(default_Model_instance):
    dmi = default_Model_instance
    
    #Setup testing rule for drug association - "A associates with R([])"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [dmi.drug_list[0]]
    r1.subject_conf = [None]
    r1.rule = ' associates with '
    r1.rule_object = [dmi.protein_list[0]]
    r1.object_conf = [[]]
    r1.check_rule_traits()
    dmi.rule_list = [r1]
    return dmi

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

# Create a default CountingSignature object for reuse in tests
@pytest.fixture()
def default_CountingSignature_instance():
    return bkcc.CountingSignature('components only')

# ------------------------------ Unit tests -----------------------------------

# ------Tests for Drug objects------

# Test if Drug objects have the required properties
def test_Drug_has_name(default_Drug_instance):
    assert hasattr(default_Drug_instance, 'name')
def test_Drug_has_symbol(default_Drug_instance):
    assert hasattr(default_Drug_instance, 'symbol')
def test_Drug_has_ID(default_Drug_instance):
    assert hasattr(default_Drug_instance, 'ID')
    
# -------Tests for Protein objects-------

# Test if Protein objects have the required properties
def test_Protein_has_name(default_Protein_instance):
    assert hasattr(default_Protein_instance, 'name')
def test_Protein_has_symbol(default_Protein_instance):
    assert hasattr(default_Protein_instance, 'symbol')
def test_Protein_has_ID(default_Protein_instance):
    assert hasattr(default_Protein_instance, 'ID')
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

# -------Tests for State objects-------

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
def test_State_has_req_protein_conf_lists(default_State_instance):
    assert hasattr(default_State_instance, 'req_protein_conf_lists')
def test_State_has_internal_links(default_State_instance):
    assert hasattr(default_State_instance, 'internal_links')

# See if we return correctly ordered lists from our state 
def test_generate_component_list_state(default_State_instance, default_Protein_instance, default_Drug_instance):
    dsi = default_State_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance
    
    # Create a state
    dsi.required_drug_list = [ddi, ddi]
    dsi.required_protein_list = [dpi, dpi]
    dsi.req_protein_conf_lists = [[0], [0,1]]
    
    # Compare the returned list to the expected one
    expected_components = [ddi, ddi, dpi, dpi]
    expected_conformations = [None, None, [0], [0,1]]
    returned_components, returned_conformations = dsi.generate_component_list()
    assert returned_components == expected_components
    assert returned_conformations == expected_conformations

# Test if the add_component_list adds a list of components/conformations successfully to a state
def test_add_component_list(default_State_instance, default_Protein_instance, default_Drug_instance):
    dsi = default_State_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Make some lists of components/conformations
    test_comp = [ddi, dpi, ddi]
    test_conf = [[None], [0,1], [None]]
    
    # Ask the state to add these
    dsi.add_component_list(test_comp, test_conf)
    
    # The lists should be added to the state
    assert dsi.required_drug_list == [ddi, ddi]
    assert dsi.required_protein_list == [dpi]
    assert dsi.req_protein_conf_lists == [[0,1]]

# Test if the state component generator is giving us the values back in correct order
def test_enumerate_components(default_State_instance, default_Protein_instance, default_Drug_instance):
    dsi = default_State_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a state
    dsi.required_drug_list = [ddi, ddi]
    dsi.required_protein_list = [dpi, dpi]
    dsi.req_protein_conf_lists = [[0], [0,1]]
    
    # Test output
    expected_index = range(4)
    expected_components = [ddi, ddi, dpi, dpi]
    
    for output_tuple, test_tuple in zip(dsi.enumerate_components(), zip(expected_index, expected_components)):
        assert output_tuple == test_tuple
        
# Test if the state component generator is giving us the values back in correct order with conformations
def test_enumerate_components_and_conformations(default_State_instance, default_Protein_instance, default_Drug_instance):
    dsi = default_State_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a state
    dsi.required_drug_list = [ddi, ddi]
    dsi.required_protein_list = [dpi, dpi]
    dsi.req_protein_conf_lists = [[0], [0,1]]
    
    # Test output
    expected_index = range(4)
    expected_components = [ddi, ddi, dpi, dpi]
    expected_conformations = [None, None, [0], [0,1]]
    
    for output_tuple, test_tuple in zip(dsi.enumerate_components(return_conformations = True), zip(expected_index, expected_components, expected_conformations)):
        assert output_tuple == test_tuple
    
# Test for pulling the correct components out of a state by number
def test_get_component_by_number(default_State_instance, default_Protein_instance, default_Drug_instance):
    dsi = default_State_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a state
    dsi.required_drug_list = [ddi, ddi]
    dsi.required_protein_list = [dpi, dpi]
    dsi.req_protein_conf_lists = [[0], [0,1]]
    
    # Test output
    assert ddi == dsi.get_component_by_number(0)
    assert ddi == dsi.get_component_by_number(1)
    assert dpi == dsi.get_component_by_number(2)
    assert dpi == dsi.get_component_by_number(3)

# -------Tests for ConformationalChange objects-------

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
    
# ------Tests for Rule objects------

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
    dri.rule_subject = [default_Drug_instance]
    dri.subject_conf = [None]
    dri.check_rule_traits() #No error expected
    
def test_Rule_check_drug_object2(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience
    
    # Set up the rule properly and call check_rule_traits
    dri.rule_object = [default_Drug_instance]
    dri.object_conf = [None]
    dri.check_rule_traits() #No error expected

def test_Rule_check_drug_subject3(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect usage, no conformations are allowed for drugs 
    dri.rule_subject = [default_Drug_instance]
    dri.subject_conf = [[0]]    
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_drug_subject4(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect usage, None, not an empty list is reqired for drugs
    dri.rule_subject = [default_Drug_instance]
    dri.subject_conf = [[]]
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_drug_object5(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect usage, no conformations are allowed for drugs 
    dri.rule_object = [default_Drug_instance]
    dri.object_conf = [[0]]    
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_drug_object6(default_Rule_instance, default_Drug_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect usage, None, not an empty list is reqired for drugs
    dri.rule_object = [default_Drug_instance]
    dri.object_conf = [[]]
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected
    
def test_Rule_check_protein_subject7(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect conf list given, None is not allowed for a protein
    dri.rule_subject = [default_Protein_instance]
    dri.subject_conf = [None]
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_protein_subject8(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Empty conf list given, empty list is allowed in rules
    dri.rule_subject = [default_Protein_instance]
    dri.subject_conf = [[]]
    dri.check_rule_traits() #No error expected

def test_Rule_check_protein_subject9(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect conf list given that contains invalid index numbers
    dri.rule_subject = [default_Protein_instance]
    dri.subject_conf = [[1, 2, 3]] # Only 0 and 1 valid for default_Protein_instance
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_protein_object10(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Usual usage of conf list 
    dri.rule_object = [default_Protein_instance]
    dri.object_conf = [[0]]
    dri.check_rule_traits() # No error expected

def test_Rule_check_protein_object11(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect conf list given that contains invalid index numbers
    dri.rule_object = [default_Protein_instance]
    dri.object_conf = [[1, 2, 3]] # Only 0 and 1 valid for default_Protein_instance
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

def test_Rule_check_protein_object12(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Empty conf list given, empty list is allowed in rules
    dri.rule_object = [default_Protein_instance]
    dri.object_conf = [[]]
    dri.check_rule_traits() #No error expected
        
def test_Rule_check_protein_object13(default_Rule_instance, default_Protein_instance):
    dri = default_Rule_instance # For typing convenience

    # Incorrect conf list given, None is not allowed for a protein
    dri.rule_object = [default_Protein_instance]
    dri.object_conf = [None]
    with pytest.raises(RuleNotValidError):
        dri.check_rule_traits() # Error expected

# Test if the generate_component_list gives back a correctly ordered list
def test_generate_component_list_rule(default_Model_irreversable_association, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversable_association
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    r1 = dmi.rule_list[0]
    
    # Already have A associates with R([]), add some extra, out of order stuff to the rule that needs to be sorted
    r1.rule_subject.insert(0, dpi)
    r1.subject_conf.insert(0, [])
    r1.rule_object.insert(0, ddi)
    r1.object_conf.insert(0, None)
    r1.rule_object.append(dpi)
    r1.object_conf.append([1])
    r1.check_rule_traits() # Double check that the changes still create a valid rule
    
    # Compare the returned list to the expected one
    # R([])A associates with AR([])R(1) should be returned as a list AAR(1)R([])R([])
    expected_components = [ddi, ddi, dpi, dpi, dpi]
    expected_conformations = [None, None, [1], [], []]
    returned_components, returned_conformations = r1.generate_component_list()
    assert returned_components == expected_components
    assert returned_conformations == expected_conformations

# Test if the generate_signature_list gives back a 'components only' signature
def test_generate_signature_list_components_only(default_Model_irreversable_association, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversable_association
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    r1 = dmi.rule_list[0]

    # We have A associates with R([]), with R able to take '' and '*' conformations
    sig_list = r1.generate_signature_list()
    
    # expect one 'components only' signature, A + R --> AR
    assert len(sig_list) == 1
    assert sig_list[0].count_type == 'components only'
    
    # A + R --> AR  Check if we got the expected result
    assert sig_list[0].subject_count[ddi] == 1
    assert sig_list[0].subject_count[dpi] == 0
    assert sig_list[0].object_count[ddi] == 0
    assert sig_list[0].object_count[dpi] == 1
    assert sig_list[0].third_state_count[ddi] == 1
    assert sig_list[0].third_state_count[dpi] == 1

# Test if the generate_signature_list gives back a 'conformations included' signature
def test_generate_signature_list_conformations_included(default_Model_irreversable_association, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversable_association
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    r1 = dmi.rule_list[0]

    # Edit so that we have A associates with R(['','*'])
    r1.object_conf = [[0,1]]
    sig_list = r1.generate_signature_list()
    
    # expect only one signature, A + R(0,1) --> AR(0,1)
    assert len(sig_list) == 1
    assert sig_list[0].count_type == 'conformations included'
    
    assert sig_list[0].subject_count[(ddi, None)] == 1
    assert sig_list[0].subject_count[(dpi, (0,1))] == 0
    assert sig_list[0].object_count[(ddi, None)] == 0
    assert sig_list[0].object_count[(dpi, (0,1))] == 1
    assert sig_list[0].third_state_count[(ddi, None)] == 1
    assert sig_list[0].third_state_count[(dpi, (0,1))] == 1
    assert sig_list[0].third_state_count[(dpi, (0))] == 0 # Check if a non-existing state is zero (actually returned as False, I think)

# Test if the generate_signature_list gives back a 'conformations included' signature
def test_generate_signature_list_any_conformations_included(default_Model_irreversable_association, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversable_association
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    r1 = dmi.rule_list[0]

    # Edit so that we have A associates with R([])R(['','*'])
    r1.rule_object = [dpi, dpi]
    r1.object_conf = [[], [0,1]]
    sig_list = r1.generate_signature_list()
    
    # Expect three signatures
    assert len(sig_list) == 3
    assert sig_list[0].count_type == 'conformations included'
    # A + R(0)R(0,1) --> AR(0)R(0,1)
    assert sig_list[0].subject_count[(ddi, None)] == 1
    assert sig_list[0].subject_count[(dpi, (0,))] == 0
    assert sig_list[0].subject_count[(dpi, (1,))] == 0
    assert sig_list[0].subject_count[(dpi, (0,1))] == 0
    assert sig_list[0].object_count[(ddi, None)] == 0
    assert sig_list[0].object_count[(dpi, (0,))] == 1
    assert sig_list[0].object_count[(dpi, (1,))] == 0
    assert sig_list[0].object_count[(dpi, (0,1))] == 1
    assert sig_list[0].third_state_count[(ddi, None)] == 1
    assert sig_list[0].third_state_count[(dpi, (0,))] == 1
    assert sig_list[0].third_state_count[(dpi, (1,))] == 0
    assert sig_list[0].third_state_count[(dpi, (0,1))] == 1
    # A + R(1)R(0,1) --> AR(1)R(0,1)
    assert sig_list[1].subject_count[(ddi, None)] == 1
    assert sig_list[1].subject_count[(dpi, (0,))] == 0
    assert sig_list[1].subject_count[(dpi, (1,))] == 0
    assert sig_list[1].subject_count[(dpi, (0,1))] == 0
    assert sig_list[1].object_count[(ddi, None)] == 0
    assert sig_list[1].object_count[(dpi, (0,))] == 0
    assert sig_list[1].object_count[(dpi, (1,))] == 1
    assert sig_list[1].object_count[(dpi, (0,1))] == 1
    assert sig_list[1].third_state_count[(ddi, None)] == 1
    assert sig_list[1].third_state_count[(dpi, (0,))] == 0
    assert sig_list[1].third_state_count[(dpi, (1,))] == 1
    assert sig_list[1].third_state_count[(dpi, (0,1))] == 1
    # A + R(0,1)R(0,1) --> AR(0,1)R(0,1)
    assert sig_list[2].subject_count[(ddi, None)] == 1
    assert sig_list[2].subject_count[(dpi, (0,))] == 0
    assert sig_list[2].subject_count[(dpi, (1,))] == 0
    assert sig_list[2].subject_count[(dpi, (0,1))] == 0
    assert sig_list[2].object_count[(ddi, None)] == 0
    assert sig_list[2].object_count[(dpi, (0,))] == 0
    assert sig_list[2].object_count[(dpi, (1,))] == 0
    assert sig_list[2].object_count[(dpi, (0,1))] == 2
    assert sig_list[2].third_state_count[(ddi, None)] == 1
    assert sig_list[2].third_state_count[(dpi, (0,))] == 0
    assert sig_list[2].third_state_count[(dpi, (1,))] == 0
    assert sig_list[2].third_state_count[(dpi, (0,1))] == 2
    
    #NOTE: Kinda complex function, perhaps more tests would be appropreate? 
    #Also, this is gonna need to be refactored to split out the "any" conformation replacemetn from the rule-specific signature generation 
    
# ------Tests for Network objects------

#Test if Network objects have the required properties
def test_Network_has_main_graph(default_Network_instance):
    assert hasattr(default_Network_instance, 'main_graph')
    
# ------Tests for CountingSignature objects------

#Test if CountingSignature objects have the required properties
def test_CountingSignature_has_count_type(default_CountingSignature_instance):
    assert hasattr(default_CountingSignature_instance, 'count_type')
def test_CountingSignature_has_subject_count(default_CountingSignature_instance):
    assert hasattr(default_CountingSignature_instance, 'subject_count')
def test_CountingSignature_has_object_count(default_CountingSignature_instance):
    assert hasattr(default_CountingSignature_instance, 'object_count')
def test_CountingSignature_has_third_state_count(default_CountingSignature_instance):
    assert hasattr(default_CountingSignature_instance, 'third_state_count')
    
# Test if we get the intended collections.Counter dictionaries when we create CountingSignatures 
def test_CountingSignature_components(default_Protein_instance, default_Drug_instance):
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a few states A + R(0,1) --> AR(0,1)
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    
    s2 = bkcc.State()
    s2.required_protein_list = [dpi]
    s2.req_protein_conf_lists = [[0,1]]
    
    s3 = bkcc.State()
    s3.required_drug_list = [ddi]
    s3.required_protein_list = [dpi]
    s3.req_protein_conf_lists = [[0,1]]
    
    # Create a signature from the states
    test_signature = bkcc.CountingSignature('components only', s1, s2, s3)
    
    # Check if we got the expected result
    assert test_signature.subject_count[ddi] == 1
    assert test_signature.subject_count[dpi] == 0
    assert test_signature.object_count[ddi] == 0
    assert test_signature.object_count[dpi] == 1
    assert test_signature.third_state_count[ddi] == 1
    assert test_signature.third_state_count[dpi] == 1

# Test if we get the intended collections.Counter dictionaries when we create CountingSignatures 
def test_CountingSignature_coformations(default_Protein_instance, default_Drug_instance):
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a few states A + R(0,1) --> AR(0,1)
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    
    s2 = bkcc.State()
    s2.required_protein_list = [dpi]
    s2.req_protein_conf_lists = [[0,1]]
    
    s3 = bkcc.State()
    s3.required_drug_list = [ddi]
    s3.required_protein_list = [dpi]
    s3.req_protein_conf_lists = [[0,1]]
    
    # Create a signature from the states
    test_signature = bkcc.CountingSignature('conformations included', s1, s2, s3)
    
    # Check if we got the expected result
    assert test_signature.subject_count[(ddi, None)] == 1
    assert test_signature.subject_count[(dpi, (0,1))] == 0
    assert test_signature.object_count[(ddi, None)] == 0
    assert test_signature.object_count[(dpi, (0,1))] == 1
    assert test_signature.third_state_count[(ddi, None)] == 1
    assert test_signature.third_state_count[(dpi, (0,1))] == 1
    assert test_signature.third_state_count[(dpi, (0))] == 0 # Check if a non-existing state is zero (actually returned as False, I think)

# Test if we get the intended collections.Counter dictionaries when we create CountingSignatures and add components directly into them
def test_CountingSignature_direct_counts(default_Protein_instance, default_Drug_instance):
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create signature from lists that are equivilant to states A + R(0,1) --> AR(0,1)
    test_signature = bkcc.CountingSignature('conformations included')
    test_signature.count_for_subject([ddi], [None])
    test_signature.count_for_object([dpi], [[0,1]])
    test_signature.count_for_third_state([ddi, dpi], [None, [0,1]])

    # Check if we got the expected result
    assert test_signature.subject_count[(ddi, None)] == 1
    assert test_signature.subject_count[(dpi, (0,1))] == 0
    assert test_signature.object_count[(ddi, None)] == 0
    assert test_signature.object_count[(dpi, (0,1))] == 1
    assert test_signature.third_state_count[(ddi, None)] == 1
    assert test_signature.third_state_count[(dpi, (0,1))] == 1
    assert test_signature.third_state_count[(dpi, (0))] == 0 # Check if a non-existing state is zero (actually returned as False, I think)
    
#continue with writing tests.....