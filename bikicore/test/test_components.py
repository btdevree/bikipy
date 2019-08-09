"""Test suite for classes in components.py

"""
import pytest
import sympy as sp
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

# Create a modified Protein object for reuse in tests
@pytest.fixture()
def modified_Protein_instance():
    mpi = bkcc.Protein()
    mpi.name = 'phosphorylated beta adrenergic receptor'
    mpi.symbol = 'pR'
    mpi.conformation_names = ['inactive', 'active']
    mpi.conformation_symbols = ['', '*']
    return mpi

# Create a default Conversion object for reuse in tests
@pytest.fixture()
def default_Conversion_instance():
    return bkcc.Conversion()
    
# Create a default Association object for reuse in tests
@pytest.fixture()
def default_Association_instance():
    return bkcc.Association()

# Create a default Dissociation object for reuse in tests
@pytest.fixture()
def default_Dissociation_instance():
    return bkcc.Dissociation()
    
# Create a default RE_Conversion object for reuse in tests
@pytest.fixture()
def default_RE_Conversion_instance():
    return bkcc.RE_Conversion()
    
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

# Create a modifed Model object from bkcm for reuse in tests
@pytest.fixture()
def modified_Model_instance(default_Model_instance, modified_Protein_instance):
    newmodel = default_Model_instance
    newmodel.name = 'modified model'
    newmodel.protein_list.append(modified_Protein_instance)
    return newmodel

# Create a default model with configured Rule object for reuse in tests
@pytest.fixture()
def default_Model_irreversible_association(default_Model_instance):
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

# Create a default model with configured Rule object for reuse in tests
@pytest.fixture()
def default_Model_irreversible_dissociation(default_Model_instance):
    dmi = default_Model_instance
    
    #Setup testing rule for a correct drug dissociation - "A dissociates from AR([])"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [dmi.drug_list[0]]
    r1.subject_conf = [None]
    r1.rule = ' dissociates from '
    r1.rule_object = [dmi.drug_list[0], dmi.protein_list[0]]
    r1.object_conf = [None, []]
    r1.check_rule_traits()
    dmi.rule_list = [r1]
    return dmi

# Create a default model with configured Rule object for reuse in tests
@pytest.fixture()
def default_Model_irreversible_conversion(default_Model_instance):
    dmi = default_Model_instance
    
    #Setup testing rule for a conformation conversion - "R(0) converts to R(1)"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [dmi.protein_list[0]]
    r1.subject_conf = [[0]]
    r1.rule = ' converts to '
    r1.rule_object = [dmi.protein_list[0]]
    r1.object_conf = [[1]]
    r1.check_rule_traits()
    dmi.rule_list = [r1]
    return dmi

# Create a modified model with configured Rule object for reuse in tests
@pytest.fixture()
def modified_Model_irreversible_conversion(modified_Model_instance):
    mmi = modified_Model_instance
    
    #Setup testing rule for a protein conversion - "R converts to pR"
    r1 = bkcc.Rule(mmi)
    r1.rule_subject = [mmi.protein_list[0]]
    r1.subject_conf = [[]]
    r1.rule = ' converts to '
    r1.rule_object = [mmi.protein_list[1]]
    r1.object_conf = [[]]
    r1.check_rule_traits()
    mmi.rule_list = [r1]
    return mmi

# Create a default Rule object for reuse in tests
@pytest.fixture()
def default_Rule_instance(default_Model_instance):
    return bkcc.Rule(default_Model_instance)

# Create a modified Rule object for reuse in tests
@pytest.fixture()
def modifed_Rule_instance(modified_Model_instance):
    return bkcc.Rule(modified_Model_instance)

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

# Create a default model with main graph for reuse
@pytest.fixture()
def default_two_state_antagonist_model_with_main_graph(default_Model_irreversible_association):
    dmi = default_Model_irreversible_association
    
    # Have rule A dissociates with R, change to reversible association
    dmi.rule_list[0].rule = ' reversibly associates with '
    
    # Add rule for a conformation conversion - "R(0) converts to R(1)"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [dmi.protein_list[0]]
    r1.subject_conf = [[0]]
    r1.rule = ' reversibly converts to '
    r1.rule_object = [dmi.protein_list[0]]
    r1.object_conf = [[1]]
    r1.check_rule_traits()
    dmi.rule_list.append(r1)
    
    # Make network on main graph return
    dmi.generate_network()

    return dmi


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
    with pytest.raises(IndexError):
        dsi.get_component_by_number(-1)
        dsi.get_component_by_number(4)

# Test for pulling the correct components and conformations out of a state by number
def test_get_component_by_number_with_conformation(default_State_instance, default_Protein_instance, default_Drug_instance):
    dsi = default_State_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a state
    dsi.required_drug_list = [ddi, ddi]
    dsi.required_protein_list = [dpi, dpi]
    dsi.req_protein_conf_lists = [[0], [0,1]]
    
    # Test output
    assert ddi, None == dsi.get_component_by_number(0, return_conformations = True)
    assert ddi, None == dsi.get_component_by_number(1, True)
    assert dpi, [0] == dsi.get_component_by_number(2, return_conformations = True)
    assert dpi, [0,1] == dsi.get_component_by_number(3, True)
    with pytest.raises(IndexError):
        dsi.get_component_by_number(-1, return_conformations = True)
        dsi.get_component_by_number(4, True)
        
# Test for autosymbol function with a drug
def test_State_autosymbol_drug(default_State_instance, default_Drug_instance):
    dsi = default_State_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a state
    dsi.required_drug_list = [ddi]
    
    # Call autosymbol
    dsi.autosymbol()
    
    # Test name
    assert dsi.symbol == 'A'
    
# Test for autosymbol function with a protein
def test_State_autosymbol_protein(default_State_instance, default_Protein_instance):
    dsi = default_State_instance
    dpi = default_Protein_instance # For typing convenience
    
    # Create a state
    dsi.required_protein_list = [dpi]
    dsi.req_protein_conf_lists = [[1]]
    
    # Call autosymbol
    dsi.autosymbol()
    
    # Test name
    assert dsi.symbol == 'R*'
    
# Test for autosymbol function with a drug-protein complex
def test_State_autosymbol_drug_and_protein(default_State_instance, default_Protein_instance, default_Drug_instance):
    dsi = default_State_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a state
    dsi.required_drug_list = [ddi]
    dsi.required_protein_list = [dpi]
    dsi.req_protein_conf_lists = [[1]]
    dsi.internal_links = [(0,1)]
    
    # Call autosymbol
    dsi.autosymbol()
    
    # Test name
    assert dsi.symbol == 'AR*'
    
# Test for autosymbol function with a complex drug-protein complex
def test_State_autosymbol_complex_drug_and_protein(default_State_instance, default_Protein_instance, default_Drug_instance, modified_Protein_instance):
    dsi = default_State_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance
    mpi = modified_Protein_instance # For typing convenience
    
    # Create a state
    dsi.required_drug_list = [ddi]
    dsi.required_protein_list = [dpi, mpi]
    dsi.req_protein_conf_lists = [[1], [0, 1]]
    dsi.internal_links = [(0, 2), (2, 1)]
    
    # Call autosymbol
    dsi.autosymbol()
    
    # Test name
    assert dsi.symbol == 'ApR,*R*'
    
    
# -------Tests for State Transition objects-------



# Test if Conversion objects have the required properties
def test_Conversion_has_name(default_Conversion_instance):
    assert hasattr(default_Conversion_instance, 'name')
def test_Conversion_has_number(default_Conversion_instance):
    assert hasattr(default_Conversion_instance, 'number')
def test_Conversion_has_ID(default_Conversion_instance):
    assert hasattr(default_Conversion_instance, 'ID')
    
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
    
# Test if RE_Conversion objects have the required properties
def test_RE_Conversion_has_name(default_RE_Conversion_instance):
    assert hasattr(default_RE_Conversion_instance, 'name')
def test_RE_Conversion_has_number(default_RE_Conversion_instance):
    assert hasattr(default_RE_Conversion_instance, 'number')
def test_RE_Conversion_has_ID(default_RE_Conversion_instance):
    assert hasattr(default_RE_Conversion_instance, 'ID')
    
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

def test_Rule_check_correct_dissociation(default_Model_instance, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience

    #Setup testing rule for a correct drug dissociation - "A dissociates from AR([])"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [ddi]
    r1.subject_conf = [None]
    r1.rule = ' dissociates from '
    r1.rule_object = [ddi, dpi]
    r1.object_conf = [None, []]
    r1.check_rule_traits()
        
def test_Rule_check_incorrect_dissociation(default_Model_instance, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience

    #Setup testing rule for incorrect drug dissociation - "A dissociates from R([])"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [ddi]
    r1.subject_conf = [None]
    r1.rule = ' dissociates from '
    r1.rule_object = [dpi]
    r1.object_conf = [[]]
    with pytest.raises(RuleNotValidError):
        r1.check_rule_traits()

# Test if the generate_component_list gives back a correctly ordered list
def test_generate_component_list_rule(default_Model_irreversible_association, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversible_association
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
    
# Test if the generate_component_list gives back a difference list
def test_generate_component_list_rule_difference_mode(default_Model_irreversible_dissociation, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversible_dissociation
    dpi = default_Protein_instance # For typing convenience
    r1 = dmi.rule_list[0]

    # Already have A dissociates from AR([])
    
    # Compare the returned list to the expected one
    expected_components = [dpi]
    expected_conformations = [[]]
    returned_components, returned_conformations = r1.generate_component_list('difference')
    assert returned_components == expected_components
    assert returned_conformations == expected_conformations

# Test if the generate_signature_list gives back a 'components only' signature
def test_generate_signature_list_components_only_association(default_Model_irreversible_association, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversible_association
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
def test_generate_signature_list_conformations_included_association(default_Model_irreversible_association, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversible_association
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

# Test if the generate_signature_list gives back multiple 'conformations included' signatures
def test_generate_signature_list_any_conformations_included_association(default_Model_irreversible_association, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversible_association
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    r1 = dmi.rule_list[0]

    # Edit so that we have A associates with R([])R(['','*'])
    r1.rule_object = [dpi, dpi]
    r1.object_conf = [[], [0,1]]
    sig_list = r1.generate_signature_list()
    
    # Expect two signatures
    assert len(sig_list) == 2
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

    #NOTE: Kinda complex function, perhaps more tests would be appropreate? 
    #Also, this is gonna need to be refactored to split out the "any" conformation replacemetn from the rule-specific signature generation 

# Test if the generate_signature_list gives back a 'components only' signature
def test_generate_signature_list_components_only_dissociation(default_Model_irreversible_dissociation, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversible_dissociation
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    r1 = dmi.rule_list[0]

    # We have A dissociates from AR([]), with R able to take '' and '*' conformations
    sig_list = r1.generate_signature_list()
    
    # expect one 'components only' signature, A dissociates from AR --> R
    assert len(sig_list) == 1
    assert sig_list[0].count_type == 'components only'
    
    # A dissociates from AR([]) --> R([])  Check if we got the expected result
    assert sig_list[0].subject_count[ddi] == 1
    assert sig_list[0].subject_count[dpi] == 0
    assert sig_list[0].object_count[ddi] == 1
    assert sig_list[0].object_count[dpi] == 1
    assert sig_list[0].third_state_count[ddi] == 0
    assert sig_list[0].third_state_count[dpi] == 1

# Test if the generate_signature_list gives back a 'conformations included' signature
def test_generate_signature_list_conformations_included_dissociation(default_Model_irreversible_dissociation, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversible_dissociation
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    r1 = dmi.rule_list[0]

    # Edit so that we have A dissociates from AR(['','*'])
    r1.object_conf = [None, [0,1]]
    sig_list = r1.generate_signature_list()
    
    # expect one 'conformations included' signature, A dissociates from AR(0,1) --> R(0,1)
    assert len(sig_list) == 1
    assert sig_list[0].count_type == 'conformations included'
    
    # A dissociates from AR(['','*']) --> R(['','*'])  Check if we got the expected result
    assert sig_list[0].subject_count[(ddi, None)] == 1
    assert sig_list[0].subject_count[(dpi, (0,1))] == 0
    assert sig_list[0].object_count[(ddi, None)] == 1
    assert sig_list[0].object_count[(dpi, (0,1))] == 1
    assert sig_list[0].third_state_count[(ddi, None)] == 0
    assert sig_list[0].third_state_count[(dpi, (0,1))] == 1
    assert sig_list[0].third_state_count[(dpi, (0))] == 0 # Check if a non-existing state is zero (actually returned as False, I think)

# Test if the generate_signature_list gives back multiple 'conformations included' signatures
def test_generate_signature_list_any_conformations_included_dissociation(default_Model_irreversible_dissociation, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_irreversible_dissociation
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    r1 = dmi.rule_list[0]

    # Edit so that we have A dissociates from AR([])R(['','*'])
    r1.rule_object = [ddi, dpi, dpi]
    r1.object_conf = [None, [], [0,1]]
    sig_list = r1.generate_signature_list()
    
    # Expect three signatures
    assert len(sig_list) == 2
    assert sig_list[0].count_type == 'conformations included'
    # A leaves AR(0)R(0,1) --> R(0)R(0,1)
    assert sig_list[0].subject_count[(ddi, None)] == 1
    assert sig_list[0].subject_count[(dpi, (0,))] == 0
    assert sig_list[0].subject_count[(dpi, (1,))] == 0
    assert sig_list[0].subject_count[(dpi, (0,1))] == 0
    assert sig_list[0].object_count[(ddi, None)] == 1
    assert sig_list[0].object_count[(dpi, (0,))] == 1
    assert sig_list[0].object_count[(dpi, (1,))] == 0
    assert sig_list[0].object_count[(dpi, (0,1))] == 1
    assert sig_list[0].third_state_count[(ddi, None)] == 0
    assert sig_list[0].third_state_count[(dpi, (0,))] == 1
    assert sig_list[0].third_state_count[(dpi, (1,))] == 0
    assert sig_list[0].third_state_count[(dpi, (0,1))] == 1
    # A leaves AR(1)R(0,1) --> R(1)R(0,1)
    assert sig_list[1].subject_count[(ddi, None)] == 1
    assert sig_list[1].subject_count[(dpi, (0,))] == 0
    assert sig_list[1].subject_count[(dpi, (1,))] == 0
    assert sig_list[1].subject_count[(dpi, (0,1))] == 0
    assert sig_list[1].object_count[(ddi, None)] == 1
    assert sig_list[1].object_count[(dpi, (0,))] == 0
    assert sig_list[1].object_count[(dpi, (1,))] == 1
    assert sig_list[1].object_count[(dpi, (0,1))] == 1
    assert sig_list[1].third_state_count[(ddi, None)] == 0
    assert sig_list[1].third_state_count[(dpi, (0,))] == 0
    assert sig_list[1].third_state_count[(dpi, (1,))] == 1
    assert sig_list[1].third_state_count[(dpi, (0,1))] == 1

# Test if the generate_signature_list gives back a 'components only' signature
def test_generate_signature_list_components_only_conversion(modified_Model_irreversible_conversion, default_Protein_instance, modified_Protein_instance):
    mmi = modified_Model_irreversible_conversion
    dpi = default_Protein_instance # For typing convenience
    mpi = modified_Protein_instance
    r1 = mmi.rule_list[0]

    # We have rule R([]) converts to pR([])"
    sig_list = r1.generate_signature_list()
    
    # Expect one 'components only' signature, R converts to pR
    assert len(sig_list) == 1
    assert sig_list[0].count_type == 'components only'
    
    # A dissociates from AR([]) --> R([])  Check if we got the expected result
    assert sig_list[0].subject_count[dpi] == 1
    assert sig_list[0].subject_count[mpi] == 0
    assert sig_list[0].object_count[dpi] == 0
    assert sig_list[0].object_count[mpi] == 1
    with pytest.raises(TypeError): # counter will still be set to None
        sig_list[0].third_state_count[dpi]
        sig_list[0].third_state_count[mpi]

# Test if the generate_signature_list gives back a 'conformations included' signature
def test_generate_signature_list_conformations_included_conversion(default_Model_irreversible_conversion, default_Protein_instance):
    dmi = default_Model_irreversible_conversion
    dpi = default_Protein_instance # For typing convenience
    r1 = dmi.rule_list[0]

    # We have rule R([0]) converts to R([1])"
    sig_list = r1.generate_signature_list()
    
    # Expect one 'components only' signature
    assert len(sig_list) == 1
    assert sig_list[0].count_type == 'conformations included'
    
    # A dissociates from AR([]) --> R([])  Check if we got the expected result
    assert sig_list[0].subject_count[(dpi, (0,))] == 1
    assert sig_list[0].subject_count[(dpi, (1,))] == 0
    assert sig_list[0].object_count[(dpi, (0,))] == 0
    assert sig_list[0].object_count[(dpi, (1,))] == 1
    with pytest.raises(TypeError): # counter will still be set to None
        sig_list[0].third_state_count[(dpi, (0,))]
        sig_list[0].third_state_count[(dpi, (1,))]

# Test if the generate_signature_list gives back many 'conformations included' signatures
def test_generate_signature_list_any_conformations_included_conversion(modified_Model_irreversible_conversion, default_Protein_instance, modified_Protein_instance):
    mmi = modified_Model_irreversible_conversion
    dpi = default_Protein_instance # For typing convenience
    mpi = modified_Protein_instance
    r1 = mmi.rule_list[0]

    # We have rule R([]) converts to pR([]), modify to R([1])R([]) converts to R([1])pR([])"
    r1.rule_subject = [dpi, dpi]
    r1.subject_conf = [[1], []]
    r1.rule_object = [dpi, mpi]
    r1.object_conf = [[1], []]
    r1.check_rule_traits()
    sig_list = r1.generate_signature_list()
    
    # Expect four signatures
    assert len(sig_list) == 4
    assert sig_list[0].count_type == 'conformations included'
    # R([1])R([0]) converts to R([1])pR([0])
    assert sig_list[0].subject_count[(dpi, (0,))] == 1
    assert sig_list[0].subject_count[(dpi, (1,))] == 1
    assert sig_list[0].subject_count[(mpi, (0,))] == 0
    assert sig_list[0].subject_count[(mpi, (1,))] == 0
    assert sig_list[0].object_count[(dpi, (0,))] == 0
    assert sig_list[0].object_count[(dpi, (1,))] == 1
    assert sig_list[0].object_count[(mpi, (0,))] == 1
    assert sig_list[0].object_count[(mpi, (1,))] == 0
    with pytest.raises(TypeError): # counter will still be set to None
        sig_list[0].third_state_count[(dpi, (0,))]
        sig_list[0].third_state_count[(dpi, (1,))]
        sig_list[0].third_state_count[(mpi, (0,))]
        sig_list[0].third_state_count[(mpi, (1,))]
    # R([1])R([0]) converts to R([1])pR([1])
    assert sig_list[1].subject_count[(dpi, (0,))] == 1
    assert sig_list[1].subject_count[(dpi, (1,))] == 1
    assert sig_list[1].subject_count[(mpi, (0,))] == 0
    assert sig_list[1].subject_count[(mpi, (1,))] == 0
    assert sig_list[1].object_count[(dpi, (0,))] == 0
    assert sig_list[1].object_count[(dpi, (1,))] == 1
    assert sig_list[1].object_count[(mpi, (0,))] == 0
    assert sig_list[1].object_count[(mpi, (1,))] == 1
    with pytest.raises(TypeError): # counter will still be set to None
        sig_list[1].third_state_count[(dpi, (0,))]
        sig_list[1].third_state_count[(dpi, (1,))]
        sig_list[1].third_state_count[(mpi, (0,))]
        sig_list[1].third_state_count[(mpi, (1,))]
    # R([1])R([1]) converts to R([1])pR([0])
    assert sig_list[2].subject_count[(dpi, (0,))] == 0
    assert sig_list[2].subject_count[(dpi, (1,))] == 2
    assert sig_list[2].subject_count[(mpi, (0,))] == 0
    assert sig_list[2].subject_count[(mpi, (1,))] == 0
    assert sig_list[2].object_count[(dpi, (0,))] == 0
    assert sig_list[2].object_count[(dpi, (1,))] == 1
    assert sig_list[2].object_count[(mpi, (0,))] == 1
    assert sig_list[2].object_count[(mpi, (1,))] == 0
    with pytest.raises(TypeError): # counter will still be set to None
        sig_list[2].third_state_count[(dpi, (0,))]
        sig_list[2].third_state_count[(dpi, (1,))]
        sig_list[2].third_state_count[(mpi, (0,))]
        sig_list[2].third_state_count[(mpi, (1,))]
    #R([1])R([1]) converts to R([1])pR([1])
    assert sig_list[3].subject_count[(dpi, (0,))] == 0
    assert sig_list[3].subject_count[(dpi, (1,))] == 2
    assert sig_list[3].subject_count[(mpi, (0,))] == 0
    assert sig_list[3].subject_count[(mpi, (1,))] == 0
    assert sig_list[3].object_count[(dpi, (0,))] == 0
    assert sig_list[3].object_count[(dpi, (1,))] == 1
    assert sig_list[3].object_count[(mpi, (0,))] == 0
    assert sig_list[3].object_count[(mpi, (1,))] == 1
    with pytest.raises(TypeError): # counter will still be set to None
        sig_list[3].third_state_count[(dpi, (0,))]
        sig_list[3].third_state_count[(dpi, (1,))]
        sig_list[3].third_state_count[(mpi, (0,))]
        sig_list[3].third_state_count[(mpi, (1,))]

# Test if the generate_signature_list gives back correct signatures
def test_generate_signature_list_any_to_conf_conversion(default_Model_irreversible_conversion, default_Protein_instance):
    dmi = default_Model_irreversible_conversion
    dpi = default_Protein_instance # For typing convenience
    r1 = dmi.rule_list[0]

    # We have rule R([0]) converts to R([1]), modify to R([]) converts to R([1])"
    r1.rule_subject = [dpi]
    r1.subject_conf = [[]]
    r1.rule_object = [dpi]
    r1.object_conf = [[1]]
    r1.check_rule_traits()
    sig_list = r1.generate_signature_list()
    
    # Expect one signature
    assert len(sig_list) == 1
    assert sig_list[0].count_type == 'conformations included'
    assert sig_list[0].subject_count[(dpi, (0,))] == 1
    assert sig_list[0].subject_count[(dpi, (1,))] == 0
    assert sig_list[0].object_count[(dpi, (0,))] == 0
    assert sig_list[0].object_count[(dpi, (1,))] == 1
    with pytest.raises(TypeError): # counter will still be set to None
        sig_list[0].third_state_count[(dpi, (0,))]
        sig_list[0].third_state_count[(dpi, (1,))]


# ------Tests for Network objects------


#Test if Network objects have the required properties
def test_Network_has_main_graph(default_Network_instance):
    assert hasattr(default_Network_instance, 'main_graph')
def test_Network_has_main_blacklist(default_Network_instance):
    assert hasattr(default_Network_instance, 'main_graph_blacklist')

# Test for autosymbol function on nodes
def test_Network_autosymbol_node(default_two_state_antagonist_model_with_main_graph):
    dam = default_two_state_antagonist_model_with_main_graph
    
    # Call autosymbol using the main graph on the model
    dam.network.autosymbol()
    
    # Sort results
    singletons = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 1]
    two_components = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 2]
    
    # The singleton states should be 1, 2, and 3. Should add comparision operators to states so they can be sorted better, but for now we sort by component length only
    acceptable_symbols = ['A', 'R', 'R*'] # we should be able to consume this without error
    for teststate in singletons:
        acceptable_symbols.remove(teststate.symbol)
    assert len(acceptable_symbols) == 0
    
    # The two component states should be 4 and 5. 
    acceptable_symbols = ['AR', 'AR*'] # we should be able to consume this without error
    for teststate in two_components:
        acceptable_symbols.remove(teststate.symbol)
    assert len(acceptable_symbols) == 0

# Test for autonumber function on nodes
def test_Network_autonumber_node(default_two_state_antagonist_model_with_main_graph):
    dam = default_two_state_antagonist_model_with_main_graph
    
    # Call autonumber using the main graph on the model
    dam.network.autonumber()
    
    # Sort results
    singletons = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 1]
    two_components = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 2]
    
    # The singleton states should be 1, 2, and 3. Should add comparision operators to states so they can be sorted better, but for now we sort by component length only
    acceptable_numbers = [1, 2, 3] # we should be able to consume this without error
    for teststate in singletons:
        acceptable_numbers.remove(teststate.number)
    assert len(acceptable_numbers) == 0
    
    # The two component states should be 4 and 5. 
    acceptable_numbers = [4, 5] # we should be able to consume this without error
    for teststate in two_components:
        acceptable_numbers.remove(teststate.number)
    assert len(acceptable_numbers) == 0

# Test for autonumber function on edges
def test_Network_autonumber_edge(default_two_state_antagonist_model_with_main_graph):
    dam = default_two_state_antagonist_model_with_main_graph
    
    # Call autonumber using the main graph on the model
    dam.network.autonumber()
    
    # Sort results
    singletons = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 1]
    two_components = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 2]
    
    singleton_tail_edges = [x for u, v, x in dam.network.main_graph.edges.data('reaction_type') if u in singletons]
    two_comp_tail_edges = [x for u, v, x in dam.network.main_graph.edges.data('reaction_type') if u in two_components]
    
    # The singleton states should have edges numbered 1, 2, 3, and one of: <-1, -2, -3>, and up to two copies of each number
    acceptable_numbers = [1, 1, 2, 2, 3, 3, -1, -1, -2, -2, -3, -3]
    for testedge in singleton_tail_edges:
        acceptable_numbers.remove(testedge.number)
    assert len(acceptable_numbers) == 6 # Six total leftovers
    assert len([x for x in acceptable_numbers if x < 0]) == 5 # 5 of which are negative
    assert len([x for x in acceptable_numbers if x > 0]) == 1 # 1 positive leftover
    
    # The one positive and negative pair left over will not be used, but a new pair of 4, -4 will be for two components edges
    leftover_value = [x for x in acceptable_numbers if x > 0][0]
    acceptable_numbers.remove(leftover_value)
    acceptable_numbers.remove(-leftover_value)
    acceptable_numbers.extend([4, -4])

    # The two_comp states should have edges that use up the remaining negative numbers and the additional pair of 4, -4
    for testedge in two_comp_tail_edges:
        acceptable_numbers.remove(testedge.number)
    assert len(acceptable_numbers) == 0
    
# Test for autovariable function on nodes
def test_Network_autovariable_node(default_two_state_antagonist_model_with_main_graph):
    dam = default_two_state_antagonist_model_with_main_graph
    
    # Call autovariable using the main graph on the model, must have numbered graphs
    dam.network.autonumber()
    dam.network.autovariable()
    
    # Sort results
    singletons = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 1]
    two_components = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 2]
    
    # The singleton states should be S1, S2, and S3. Should add comparision operators to states so they can be sorted better, but for now we sort by component length only
    acceptable_vars = [*sp.symbols('S_1, S_2, S_3')] # we should be able to consume this without error
    for teststate in singletons:
        acceptable_vars.remove(teststate.variable)
    assert len(acceptable_vars) == 0
    
    # The two-component states should be S4 and S5. 
    acceptable_vars = [*sp.symbols('S_4, S_5')] # we should be able to consume this without error
    for teststate in two_components:
        acceptable_vars.remove(teststate.variable)
    assert len(acceptable_vars) == 0

# Test for autovariable function on edges
def test_Network_autovariable_edge(default_two_state_antagonist_model_with_main_graph):
    dam = default_two_state_antagonist_model_with_main_graph
    
     # Call autovariable using the main graph on the model, must have numbered graphs
    dam.network.autonumber()
    dam.network.autovariable()
    
    # Sort results
    singletons = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 1]
    two_components = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 2]
    
    singleton_tail_edges = [x for u, v, x in dam.network.main_graph.edges.data('reaction_type') if u in singletons]
    two_comp_tail_edges = [x for u, v, x in dam.network.main_graph.edges.data('reaction_type') if u in two_components]
    
    # The singleton states should have edges with variables k1, k2, k3, and k-3 - I am pretty sure the example network must always be in this order, but check order if it fails
    acceptable_vars = [*sp.symbols('k_1, k_1, k_2, k_2, k_3, k_-3')]
    for testedge in singleton_tail_edges:
        print(testedge, testedge.variable)
        acceptable_vars.remove(testedge.variable)
    assert len(acceptable_vars) == 0
    
    # Edges coming from two-component states should be k-1, k-2, k3, k-3, k4, and k-4
    acceptable_vars = [*sp.symbols('k_-1, k_-1, k_-2, k_-2, k_4, k_-4')]
    for testedge in two_comp_tail_edges:
        acceptable_vars.remove(testedge.variable)
    assert len(acceptable_vars) == 0
    
# Test for autoname function on nodes
def test_Network_autoname_node(default_two_state_antagonist_model_with_main_graph):
    dam = default_two_state_antagonist_model_with_main_graph
    
    # Call autoname using the main graph on the model, must have numbered graphs
    dam.network.autonumber()
    dam.network.autoname()
    
    # Sort results
    singletons = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 1]
    two_components = [x for x in dam.network.main_graph.__iter__() if len(x.required_drug_list) + len(x.required_protein_list) == 2]
    
    # The singleton states should be S1, S2, and S3. Should add comparision operators to states so they can be sorted better, but for now we sort by component length only
    acceptable_names = ['State 1', 'State 2', 'State 3'] # we should be able to consume this without error
    for teststate in singletons:
        acceptable_names.remove(teststate.name)
    assert len(acceptable_names) == 0
    
    # The two-component states should be S4 and S5. 
    acceptable_names = ['State 4', 'State 5'] # we should be able to consume this without error
    for teststate in two_components:
        acceptable_names.remove(teststate.name)
    assert len(acceptable_names) == 0
    
    
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
def test_CountingSignature_components_association(default_Protein_instance, default_Drug_instance):
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a few states A + R(0,1) --> AR(0,1), rule not actually tested here, just signature function
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
def test_CountingSignature_conformations(default_Protein_instance, default_Drug_instance):
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a few states A + R(0,1) --> AR(0,1), rule not actually tested here, just signature function
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


# ------Tests for  objects------


