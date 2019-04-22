"""Test suite for classes in model.py

"""
import pytest
import collections
import networkx as nx
import bikipy.bikicore.model as bkcm
import bikipy.bikicore.components as bkcc

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

# Create a default Model object with empty graph from bkcm for reuse in tests
@pytest.fixture()
def default_Model_instance(default_Drug_instance, default_Protein_instance):
    newmodel = bkcm.Model(1, 'Default Null Model', None)
    newmodel.drug_list.append(default_Drug_instance)
    newmodel.protein_list.append(default_Protein_instance)
    return newmodel

@pytest.fixture()
def default_Model_list():
    model_list = [bkcm.Model(1, 'Default-1', None),
                bkcm.Model(2, 'Default-2', None),
                bkcm.Model(4, 'Default-4', None)]
    return model_list

# Run a lot of tests about state matching with a simple association rule
@pytest.fixture()
def model_for_matching_tests(default_Model_instance):
    dmi = default_Model_instance
    
    # Create a new Network object without calling dmi.generate_network()
    dmi.network = bkcc.Network()
        
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
    
# ---- Unit tests ----

# --Tests for Model objects--

# Test if Model object has the required properties
def test_Model_has_number(default_Model_instance):
    assert hasattr(default_Model_instance, 'number')
    
def test_Model_has_name(default_Model_instance):
    assert hasattr(default_Model_instance, 'name')

def test_Model_has_parent(default_Model_instance):
    assert hasattr(default_Model_instance, 'parent_model')

def test_Model_has_ID(default_Model_instance):
    assert hasattr(default_Model_instance, 'ID')

def test_Model_has_drug_list(default_Model_instance):
    assert hasattr(default_Model_instance, 'drug_list')
    
def test_Model_has_protein_list(default_Model_instance):
    assert hasattr(default_Model_instance, 'protein_list')

def test_Model_has_rule_list(default_Model_instance):
    assert hasattr(default_Model_instance, 'rule_list')

def test_Model_has_compartment_list(default_Model_instance):
    assert not hasattr(default_Model_instance, 'compartment_list') #not implemented yet

# -- Tests for model class methods --

# A graph of non-connected singletons should be made with an empty rule list
def test_Model_generate_network_null(default_Model_instance):
    dmi = default_Model_instance
    dmi.generate_network()
    
    # Create comparision graph shape
    testgraph = nx.DiGraph([])
    testgraph.add_nodes_from([1,2])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(dmi.network.main_graph, testgraph)
    
# Test that the 'associates with' rule creates a valid shaped graph 
def test_Model_generate_network_Irr_association(default_Model_instance):
    dmi = default_Model_instance
    
    #Setup rule for simple drug association - "A associates with R"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [dmi.drug_list[0]]
    r1.subject_conf = [None]
    r1.rule = ' associates with '
    r1.rule_object = [dmi.protein_list[0]]
    r1.object_conf = [[0,1]]
    r1.check_rule_traits()
    
    # Attach rule to model and generate network
    dmi.rule_list = [r1]
    dmi.generate_network()
    
    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1,2,3])
    testgraph.add_edges_from([(1,3), (2,3)])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(dmi.network.main_graph, testgraph)
    
# Test that the 'associates with' rule creates a valid shaped graph without unwanted dimerization
def test_Model_generate_network_Irr_association_dimers(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mmt.drug_list[0]]
    s1.required_protein_list = [mmt.protein_list[0]]
    s1.req_protein_conf_lists = [[0]]
    mmt.network.main_graph.add_node(s1)
    
    s2 = bkcc.State()
    s2.required_drug_list = []
    s2.required_protein_list = [mmt.protein_list[0]]
    s2.req_protein_conf_lists = [[0]]
    mmt.network.main_graph.add_node(s2)
    
    # Run method
    subject_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'subject')
    object_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'object')
    
    # Check output to double-check that the states match as intended
    assert subject_list == [s1]
    assert object_list == [s1, s2]
    
    # Apply the existing association rule - "A associates with R"
    mmt.apply_rules_to_network()
    
    # Create comparision graph shape
    # We should not allow AR to associate with R to make ARR
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1,2])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(mmt.network.main_graph, testgraph)
   
# -- Helper method tests in model.py --    
            
# Test for creation methods
def test_create_new(default_Model_list):
    # Add in when copy logic is created
    pass
    
# Test logic for private helper methods
def test_modelnum_finder(default_Model_list):
    # Choose first positive integer not already in the list
    next_number = bkcm._find_next_model_number(default_Model_list)
    assert next_number == 3
    
    # Add a model with number 3 to the list, next number should be 5
    default_Model_list.append(bkcm.Model(3, 'Default-3', None))
    next_number = bkcm._find_next_model_number(default_Model_list)
    assert next_number == 5
    
    # Delete first model from list, next number should be 1
    del default_Model_list[0]
    next_number = bkcm._find_next_model_number(default_Model_list)
    assert next_number == 1

# If there are no states in the graph, empty lists of states should be returned
def test1_find_states_that_match_rule(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    # Empty network, no code here
    
    # Run method
    subject_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'subject')
    object_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'object')
    
    # Check output
    assert subject_list == []
    assert object_list == []
    
# A should be returned as subject
def test2_find_states_that_match_rule(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mmt.drug_list[0]]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = []
    mmt.network.main_graph.add_node(s1)
    
    # Run method
    subject_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'subject')
    object_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'object')
    
    # Check output
    assert subject_list == [s1]
    assert object_list == []
    
# R should be returned as object
def test3_find_states_that_match_rule(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = []
    s1.required_protein_list = [mmt.protein_list[0]]
    s1.req_protein_conf_lists = [[]]
    mmt.network.main_graph.add_node(s1)
    
    # Run method
    subject_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'subject')
    object_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'object')
    
    # Check output
    assert subject_list == []
    assert object_list == [s1]
    
# Both and R and A should be returned
def test4_find_states_that_match_rule(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mmt.drug_list[0]]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = []
    mmt.network.main_graph.add_node(s1)

    s2 = bkcc.State()
    s2.required_drug_list = []
    s2.required_protein_list = [mmt.protein_list[0]]
    s2.req_protein_conf_lists = [[]]
    mmt.network.main_graph.add_node(s2)
    
   # Run method
    subject_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'subject')
    object_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'object')
    
    # Check output
    assert subject_list == [s1]
    assert object_list == [s2]

# If the graph is already associated, R and A should be returned as object and subject, respectivly, plus the complex returned for both
def test5_find_states_that_match_rule(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mmt.drug_list[0]]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = []
    mmt.network.main_graph.add_node(s1)
    
    s2 = bkcc.State()
    s2.required_drug_list = []
    s2.required_protein_list = [mmt.protein_list[0]]
    s2.req_protein_conf_lists = [[]]
    mmt.network.main_graph.add_node(s2)
    
    s3 = bkcc.State()
    s3.required_drug_list = [mmt.drug_list[0]]
    s3.required_protein_list = [mmt.protein_list[0]]
    s3.req_protein_conf_lists = [[]]
    mmt.network.main_graph.add_node(s3)
    
    # Run method
    subject_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'subject')
    object_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'object')
    
    # Check output - Note that we do not care about the order of the lists
    assert collections.Counter(subject_list) == collections.Counter([s1, s3])
    assert collections.Counter(object_list) == collections.Counter([s2, s3])

# R(0, 1) should be returned as object
def test6_find_states_that_match_rule(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = []
    s1.required_protein_list = [mmt.protein_list[0]]
    s1.req_protein_conf_lists = [[0]]
    mmt.network.main_graph.add_node(s1)
    
    # Setup test network
    s2 = bkcc.State()
    s2.required_drug_list = []
    s2.required_protein_list = [mmt.protein_list[0]]
    s2.req_protein_conf_lists = [[1]]
    mmt.network.main_graph.add_node(s2)
    
    # Setup test network
    s3 = bkcc.State()
    s3.required_drug_list = []
    s3.required_protein_list = [mmt.protein_list[0]]
    s3.req_protein_conf_lists = [[0, 1]]
    mmt.network.main_graph.add_node(s3)
    
    # Edit rule
    mmt.rule_list[0].object_conf = [[0,1]]
    
    # Run method
    subject_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'subject')
    object_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'object')
    
    # Check output
    assert subject_list == []
    assert object_list == [s3]
    
# R should not be returned as matching a dimer in the rule
def test7_find_states_that_match_rule(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = []
    s1.required_protein_list = [mmt.protein_list[0]]
    s1.req_protein_conf_lists = [[0]]
    mmt.network.main_graph.add_node(s1)
    
    # Edit rule to associate A with RR
    mmt.rule_list[0].rule_object = [mmt.protein_list[0], mmt.protein_list[0]]
    mmt.rule_list[0].object_conf = [[], []]
    
    # Run method
    subject_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'subject')
    object_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'object')
    
    # Check output, no states should be found
    assert subject_list == []
    assert object_list == []
    
# Do we return valid state pairs with a simple signature
def test_find_association_pairs(model_for_matching_tests, default_Protein_instance, default_Drug_instance):
    mmt = model_for_matching_tests    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # We already have the rule A + R([]) --> AR([]), get the signature
    ref_sig = mmt.rule_list[0].generate_signature_list()
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    
    s2 = bkcc.State()
    s2.required_protein_list = [dpi]
    s2.req_protein_conf_lists = [[0,1]]
    
    s3 = bkcc.State()
    s3.required_drug_list = [ddi]
    s3.required_protein_list = [dpi]
    s3.req_protein_conf_lists = [[0,1]]
    
    # Make some state lists
    test_sub_list = [s1, s3]
    test_obj_list = [s2, s3]
    
    # Run function
    valid_tuples = mmt._find_association_pairs(ref_sig, test_sub_list, test_obj_list)
    
    # Check if we got the expected result
    assert valid_tuples == [(s1, s2)]
    
    # Do we return valid state pairs with a conformation signature
def test_find_association_pairs_with_conformation(model_for_matching_tests, default_Protein_instance, default_Drug_instance):
    mmt = model_for_matching_tests    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Modify rule to give A + R([1]) --> AR([1]), get the signature
    mmt.rule_list[0].object_conf = [[1]]
    ref_sig = mmt.rule_list[0].generate_signature_list()
    print(ref_sig[0].count_type)
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    
    s2 = bkcc.State()
    s2.required_protein_list = [dpi]
    s2.req_protein_conf_lists = [[1]]
    
    s3 = bkcc.State()
    s3.required_drug_list = [ddi]
    s3.required_protein_list = [dpi]
    s3.req_protein_conf_lists = [[1]]
    
    s4 = bkcc.State()
    s4.required_protein_list = [dpi]
    s4.req_protein_conf_lists = [[0]]
    
    s5 = bkcc.State()
    s5.required_drug_list = [ddi]
    s5.required_protein_list = [dpi]
    s5.req_protein_conf_lists = [[0]]
    
    # Make some state lists (not all these should actually be identified by the matching code, but should be rejected nontheless)
    test_sub_list = [s1, s3, s4, s5]
    test_obj_list = [s2, s3, s4, s5]
    
    # Run function
    valid_tuples = mmt._find_association_pairs(ref_sig, test_sub_list, test_obj_list)
    print(valid_tuples)
    print('s5', s5)
    # Check if we got the expected result
    assert valid_tuples == [(s1, s2)]

    