"""Test suite for classes in model.py

"""
import pytest
import networkx as nx
import bikipy.bikicore.model as bkcm
import bikipy.bikicore.components as bkcc

#---- Testing fixtures ----

# Create a default Model object for reuse in tests
@pytest.fixture()
def default_Model_instance():
    return bkcm.Model(1, 'Default-1', None)
    
@pytest.fixture()
def default_Model_list():
    model_list = [bkcm.Model(1, 'Default-1', None),
                bkcm.Model(2, 'Default-2', None),
                bkcm.Model(4, 'Default-4', None)]
    return model_list
    
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

# Create a default Model object from bkcm for reuse in tests
@pytest.fixture()
def default_Model_instance(default_Drug_instance, default_Protein_instance):
    newmodel = bkcm.Model(1, 'default model', None)
    newmodel.drug_list.append(default_Drug_instance)
    newmodel.protein_list.append(default_Protein_instance)
    return newmodel
    
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

# An empty graph should be made with an empty rule list
def test_Model_generate_network_null(default_Model_instance):
    dmi = default_Model_instance
    dmi.generate_network()
    
    # Create comparision graph shape
    testgraph = nx.DiGraph([]) #empty graph
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(dmi.network.main_graph, testgraph)
    
# Test that the 'associates with' rule creates a valid graph 
def test_Model_generate_network_Irr_association(default_Model_instance):
    dmi = default_Model_instance
    
    #Setup rule for simple drug association
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = dmi.drug_list[0]
    r1.subject_conf = None
    r1.rule = ' associates with '
    r1.rule_object = dmi.protein_list[0]
    r1.object_conf = 'all'
    r1.check_rule_traits()
    
    # Attach rule to model and generate network
    dmi.rule_list = [r1]
    dmi.generate_network()
    
    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1,2])
    testgraph.add_edge(1,2)
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(dmi.network.main_graph, testgraph)


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

    
    