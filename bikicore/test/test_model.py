"""Test suite for classes in model.py

"""
import pytest
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
    
    # Test that the 'associates with' rule creates a valid shaped graph 

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
def test1_find_states_with_matching_components(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    # Empty network, no code here
    
    # Run method
    subject_list, object_list = mmt._find_states_with_matching_components(mmt.rule_list[0])
    
    # Check output
    assert subject_list == []
    assert object_list == []
    
# A should be returned as subject
def test2_find_states_with_matching_components(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mmt.drug_list[0]]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = [[]]
    mmt.network.main_graph.add_node(s1)
    
    # Run method
    subject_list, object_list = mmt._find_states_with_matching_components(mmt.rule_list[0])
    
    # Check output
    assert subject_list == [s1]
    assert object_list == []
    
# R should be returned as object
def test3_find_states_with_matching_components(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = []
    s1.required_protein_list = [mmt.protein_list[0]]
    s1.req_protein_conf_lists = [[]]
    mmt.network.main_graph.add_node(s1)
    
    # Run method
    subject_list, object_list = mmt._find_states_with_matching_components(mmt.rule_list[0])
    
    # Check output
    assert subject_list == []
    assert object_list == [s1]
    
# Both and R and A should be returned
def test4_find_states_with_matching_components(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mmt.drug_list[0]]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = [[]]
    mmt.network.main_graph.add_node(s1)

    s2 = bkcc.State()
    s2.required_drug_list = []
    s2.required_protein_list = [mmt.protein_list[0]]
    s2.req_protein_conf_lists = [[]]
    mmt.network.main_graph.add_node(s2)
    
    # Run method
    subject_list, object_list = mmt._find_states_with_matching_components(mmt.rule_list[0])
    
    # Check output
    assert subject_list == [s1]
    assert object_list == [s2]

# If the graph is already associated, R and A should be returned as object and subject, respectivly, plus the complex returned for both
def test5_find_states_with_matching_components(model_for_matching_tests):
    mmt = model_for_matching_tests
        
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mmt.drug_list[0]]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = [[]]
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
    subject_list, object_list = mmt._find_states_with_matching_components(mmt.rule_list[0])
    
    # Check output - Note that we do not care about the order of the lists
    assert set(subject_list) == set([s1, s3])
    assert set(object_list) == set([s2, s3])

# R(0, 1) should be returned as object
def test6_find_states_with_matching_components(model_for_matching_tests):
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
    print(mmt.rule_list[0].object_conf)
    mmt.rule_list[0].object_conf = [[0,1]]
    print(mmt.rule_list[0].object_conf)
    
    # Run method
    subject_list, object_list = mmt._find_states_with_matching_components(mmt.rule_list[0])
    
    # Check output
    assert subject_list == []
    assert object_list == [s3]

# Test out the counting helper method
def test_count_components(default_Model_instance):
    dmi = default_Model_instance
    
    # Create a list of components
    test_drug = dmi.drug_list[0]
    test_protein = dmi.protein_list[0]
    
    test_list = [test_drug, test_protein, test_drug]
    
    # Ask the function to count the components
    count_dict = dmi._count_components(test_list)
    
    # See if it returns what we expect
    assert count_dict[test_drug] == 2
    assert count_dict[test_protein] == 1
    
    