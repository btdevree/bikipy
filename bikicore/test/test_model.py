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

# Create a second Drug object for reuse in tests
@pytest.fixture()
def second_Drug_instance():
    ddi = bkcc.Drug()
    ddi.name = 'alprenolol'
    ddi.symbol = 'B'
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

# Create a secondary Model object with empty graph from bkcm for reuse in tests
@pytest.fixture()
def second_Model_instance(default_Drug_instance, second_Drug_instance, default_Protein_instance):
    newmodel = bkcm.Model(1, 'Second Null Model', None)
    newmodel.drug_list.append(default_Drug_instance)
    newmodel.drug_list.append(second_Drug_instance)
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

# Setup for many dissociation tests
@pytest.fixture()
def model_for_dissociation_matching(default_Model_instance):
    dmi = default_Model_instance
    
    # Create a new Network object without calling dmi.generate_network()
    dmi.network = bkcc.Network()
    
    #Setup rule for simple drug disassociation - "A disassociates from AR"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [dmi.drug_list[0]]
    r1.subject_conf = [None]
    r1.rule = ' dissociates from '
    r1.rule_object = [dmi.drug_list[0], dmi.protein_list[0]]
    r1.object_conf = [None, []]
    r1.check_rule_traits()
    dmi.rule_list = [r1]
    return dmi

# Setup for conversion tests
@pytest.fixture()
def default_Model_conversion(default_Model_instance):
    dmi = default_Model_instance
    
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [dmi.protein_list[0]]
    r1.subject_conf = [[0]]
    r1.rule = ' converts to '
    r1.rule_object = [dmi.protein_list[0]]
    r1.object_conf = [[1]]
    r1.check_rule_traits()
    dmi.rule_list = [r1]
    return dmi

@pytest.fixture()
def default_Model_two_rule_antagonist(default_Model_instance):
    dmi = default_Model_instance
    
    #Setup rule for an antagonist association - "A associates reversibly with R([])"
    r0 = bkcc.Rule(dmi)
    r0.rule_subject = [dmi.drug_list[0]]
    r0.subject_conf = [None]
    r0.rule = ' reversibly associates with '
    r0.rule_object = [dmi.protein_list[0]]
    r0.object_conf = [[]]
    r0.check_rule_traits()

    #Setup testing rule for a state conversion - "R(0) reversibly converts to R(1)"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [dmi.protein_list[0]]
    r1.subject_conf = [[0]]
    r1.rule = ' reversibly converts to '
    r1.rule_object = [dmi.protein_list[0]]
    r1.object_conf = [[1]]
    r1.check_rule_traits()
    
    # Add to model
    dmi.rule_list = [r0, r1]
    return dmi

@pytest.fixture()
def default_Model_three_rule_antagonists(second_Model_instance):
    smi = second_Model_instance
    
    #Setup rule for an antagonist association - "A associates reversibly with R([])"
    r0 = bkcc.Rule(smi)
    r0.rule_subject = [smi.drug_list[0]]
    r0.subject_conf = [None]
    r0.rule = ' reversibly associates with '
    r0.rule_object = [smi.protein_list[0]]
    r0.object_conf = [[]]
    r0.check_rule_traits()

    #Setup testing rule for a state conversion - "R(0) reversibly converts to R(1)"
    r1 = bkcc.Rule(smi)
    r1.rule_subject = [smi.protein_list[0]]
    r1.subject_conf = [[0]]
    r1.rule = ' reversibly converts to '
    r1.rule_object = [smi.protein_list[0]]
    r1.object_conf = [[1]]
    r1.check_rule_traits()
    
    #Setup rule for an antagonist association - "B associates reversibly with R([])"
    r2 = bkcc.Rule(smi)
    r2.rule_subject = [smi.drug_list[1]]
    r2.subject_conf = [None]
    r2.rule = ' reversibly associates with '
    r2.rule_object = [smi.protein_list[0]]
    r2.object_conf = [[]]
    r2.check_rule_traits()
    
    # Add to model
    smi.rule_list = [r0, r1, r2]
    return smi

@pytest.fixture()
def default_Model_four_rule_competitive_antagonists(second_Model_instance):
    smi = second_Model_instance
    
    #Setup rule for an antagonist association - "A associates reversibly with R([])"
    r0 = bkcc.Rule(smi)
    r0.rule_subject = [smi.drug_list[0]]
    r0.subject_conf = [None]
    r0.rule = ' reversibly associates with '
    r0.rule_object = [smi.protein_list[0]]
    r0.object_conf = [[]]
    r0.check_rule_traits()

    #Setup testing rule for a state conversion - "R(0) reversibly converts to R(1)"
    r1 = bkcc.Rule(smi)
    r1.rule_subject = [smi.protein_list[0]]
    r1.subject_conf = [[0]]
    r1.rule = ' reversibly converts to '
    r1.rule_object = [smi.protein_list[0]]
    r1.object_conf = [[1]]
    r1.check_rule_traits()
    
    #Setup rule for an antagonist association - "B associates reversibly with R([])"
    r2 = bkcc.Rule(smi)
    r2.rule_subject = [smi.drug_list[1]]
    r2.subject_conf = [None]
    r2.rule = ' reversibly associates with '
    r2.rule_object = [smi.protein_list[0]]
    r2.object_conf = [[]]
    r2.check_rule_traits()
    
    #Setup testing rule for a competition - "A is competitive with B"
    r3 = bkcc.Rule(smi)
    r3.rule_subject = [smi.drug_list[0]]
    r3.subject_conf = [None]
    r3.rule = ' is competitive with '
    r3.rule_object = [smi.drug_list[1]]
    r3.object_conf = [None]
    r3.check_rule_traits()
    
    # Add to model
    smi.rule_list = [r0, r1, r2, r3]
    return smi

@pytest.fixture()
def default_Model_competition(second_Model_instance):
    smi = second_Model_instance
    
    #Setup testing rule for a competition - "A is competitive with B"
    r0 = bkcc.Rule(smi)
    r0.rule_subject = [smi.drug_list[0]]
    r0.subject_conf = [None]
    r0.rule = ' is competitive with '
    r0.rule_object = [smi.drug_list[1]]
    r0.object_conf = [None]
    r0.check_rule_traits()
    
    # Add to model
    smi.rule_list = [r0]
    return smi


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
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3]) # A, R(0), R(1)
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(dmi.network.main_graph, testgraph)
    
# Test that the 'associates with' rule creates a valid shaped graph 
def test_Model_apply_rules_irr_association(model_for_matching_tests):
    mmt = model_for_matching_tests
    
    # Has rule for simple drug association - "A associates with R"
    
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mmt.drug_list[0]]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = []
    mmt.network.main_graph.add_node(s1)
    
    s2 = bkcc.State()
    s2.required_drug_list = []
    s2.required_protein_list = [mmt.protein_list[0]]
    s2.req_protein_conf_lists = [[0]]
    mmt.network.main_graph.add_node(s2)
    
    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3])
    testgraph.add_edges_from([(1, 3), (2, 3)])
    
    # Apply the existing association rule, bypass singleton generation with .generate_network() 
    mmt.apply_rules_to_network()
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(mmt.network.main_graph, testgraph)
    
# Test that the 'associates with' rule creates a valid shaped graph without unwanted dimerization
def test_Model_apply_rules_irr_association_dimers(model_for_matching_tests):
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
    
# Test that the 'dissociates from' rule creates a valid shaped graph 
def test_Model_apply_rules_irr_disassociation(model_for_dissociation_matching):
    mdm = model_for_dissociation_matching
    
    # Have rule for simple drug dissociation - "A dissociates from AR([])"
    
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mdm.drug_list[0]]
    s1.required_protein_list = [mdm.protein_list[0]]
    s1.req_protein_conf_lists = [[0, 1]]
    s1.internal_links = [(0, 1)]
    mdm.network.main_graph.add_node(s1)
    
    # Apply the existing association rule - "A dissociates from AR"
    mdm.apply_rules_to_network()

    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1,2,3])
    testgraph.add_edges_from([(3,1), (3,2)])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(mdm.network.main_graph, testgraph)

# Test that the 'reversibly dissociates from' rule creates a valid shaped graph 
def test_Model_apply_rules_rev_disassociation(model_for_dissociation_matching):
    mdm = model_for_dissociation_matching
    
    # Have rule for simple drug disassociation - "A dissociates from AR([])", edit to be reversible
    mdm.rule_list[0].rule = ' reversibly dissociates from '
    
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mdm.drug_list[0]]
    s1.required_protein_list = [mdm.protein_list[0]]
    s1.req_protein_conf_lists = [[0, 1]]
    s1.internal_links = [(0, 1)]
    mdm.network.main_graph.add_node(s1)
    
    # Apply the existing association rule - "A dissociates from AR"
    mdm.apply_rules_to_network()

    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3])
    testgraph.add_edges_from([(3, 1), (3, 2), (1, 3), (2, 3)])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(mdm.network.main_graph, testgraph)
    
# Test that the 'reversibly associates with' rule creates a valid shaped graph 
def test_Model_apply_rules_rev_association(model_for_matching_tests):
    mmt = model_for_matching_tests
    
    # Has rule for simple drug association, modify to reversible association - "A reversibly associates with R"
    mmt.rule_list[0].rule = ' reversibly associates with '
    
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mmt.drug_list[0]]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = []
    mmt.network.main_graph.add_node(s1)
    
    s2 = bkcc.State()
    s2.required_drug_list = []
    s2.required_protein_list = [mmt.protein_list[0]]
    s2.req_protein_conf_lists = [[0]]
    mmt.network.main_graph.add_node(s2)
        
    # Apply the association rule
    mmt.apply_rules_to_network()
    
    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3])
    testgraph.add_edges_from([(1, 3), (2, 3), (3, 1), (3, 2)])

    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(mmt.network.main_graph, testgraph)

# Test that the ' associates and dissociates in rapid equlibrium with ' rule creates a valid shaped graph 
def test_Model_apply_rules_rev_RE_association(model_for_matching_tests):
    mmt = model_for_matching_tests
    
    # Has rule for simple drug association, modify to RE association - "A associates and dissociates in rapid equlibrium with R"
    mmt.rule_list[0].rule = ' associates and dissociates in rapid equlibrium with '
    
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mmt.drug_list[0]]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = []
    mmt.network.main_graph.add_node(s1)
    
    s2 = bkcc.State()
    s2.required_drug_list = []
    s2.required_protein_list = [mmt.protein_list[0]]
    s2.req_protein_conf_lists = [[0]]
    mmt.network.main_graph.add_node(s2)
        
    # Apply the association rule
    mmt.apply_rules_to_network()
    
    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3])
    testgraph.add_edges_from([(1, 3), (2, 3), (3, 1), (3, 2)])

    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(mmt.network.main_graph, testgraph)

# Test that the ' dissociates and reassociates in rapid equlibrium with ' rule creates a valid shaped graph 
def test_Model_apply_rules_rev_RE_disassociation(model_for_dissociation_matching):
    mdm = model_for_dissociation_matching
    
    # Have rule for simple drug disassociation - "A disassociates from AR([])", edit to be reversible in rapid equlibrium
    mdm.rule_list[0].rule = ' dissociates and reassociates in rapid equlibrium from '
    
    # Setup test network
    s1 = bkcc.State()
    s1.required_drug_list = [mdm.drug_list[0]]
    s1.required_protein_list = [mdm.protein_list[0]]
    s1.req_protein_conf_lists = [[0, 1]]
    s1.internal_links = [(0, 1)]
    mdm.network.main_graph.add_node(s1)
    
    # Apply the existing association rule - "A disassociates from AR([])"
    mdm.apply_rules_to_network()

    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3])
    testgraph.add_edges_from([(3, 1), (3, 2), (1, 3), (2, 3)])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(mdm.network.main_graph, testgraph)

# Test that the ' is competitive with ' rule creates a valid shaped graph 
def test_Model_apply_rules_competition(second_Model_instance, default_Drug_instance, second_Drug_instance, default_Protein_instance):
    smi = second_Model_instance
    ddi = default_Drug_instance 
    sdi = second_Drug_instance
    dpi = default_Protein_instance# For typing convenience
    
    #Setup testing rule for a competition - "A is competitive with B"
    r0 = bkcc.Rule(smi)
    r0.rule_subject = [smi.drug_list[0]]
    r0.subject_conf = [None]
    r0.rule = ' is competitive with '
    r0.rule_object = [smi.drug_list[1]]
    r0.object_conf = [None]
    r0.check_rule_traits()
    smi.rule_list = [r0]

    # Create a new Network object without calling dmi.generate_network()
    smi.network = bkcc.Network()
        
    # Setup test network
    s1 = bkcc.State() # Allowed
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi]
    s1.req_protein_conf_lists = [[0]]
    s1.internal_links = [(0, 1)]
    smi.network.main_graph.add_node(s1)
    
    s2 = bkcc.State() # Allowed
    s2.required_drug_list = [sdi]
    s2.required_protein_list = [dpi]
    s2.req_protein_conf_lists = [[0]]
    s2.internal_links = [(0, 1)]
    smi.network.main_graph.add_node(s2)

    s3 = bkcc.State() # Should be deleted
    s3.required_drug_list = [ddi, sdi]
    s3.required_protein_list = [dpi]
    s3.req_protein_conf_lists = [[0]]
    s3.internal_links = [(0, 2), (1, 2)]
    smi.network.main_graph.add_node(s3)
    
    # Apply the existing rule - "A is competitive with B"
    smi.apply_rules_to_network()

    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2]) # State 3 should be deleted
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(smi.network.main_graph, testgraph)
    assert smi.network.main_graph_blacklist == [s3] # Make sure this illegal state is on the blacklist

# Test that generate_network can process multiple rules
def test_Model_generate_network_two_rules(default_Model_two_rule_antagonist):
    m2r = default_Model_two_rule_antagonist
    m2r.generate_network()

    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3, 4, 5]) #1=A, 2=R, 3=R*, 4=AR, 5=AR*
    testgraph.add_edges_from([(1, 4), (4, 1), (2, 4), (4, 2), (1, 5), (5, 1), (3, 5), (5, 3), (2, 3), (3, 2), (4, 5), (5, 4)])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(m2r.network.main_graph, testgraph)

# Test that generate_network can process multiple rules
def test_Model_generate_network_three_rules(default_Model_three_rule_antagonists):
    m3r = default_Model_three_rule_antagonists
    m3r.generate_network()

    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) #1=A, 2=B, 3=R, 4=R*, 5=AR, 6=AR*, 7=BR, 8=BR*, 9=ABR, 10=ABR*
    testgraph.add_edges_from([(1, 5), (5, 1), (3, 5), (5, 3), (2, 7), (7, 2), (3, 7), (7, 3),\
                              (1, 6), (6, 1), (4, 6), (6, 4), (2, 8), (8, 2), (4, 8), (8, 4),\
                              (1, 9), (9, 1), (7, 9), (9, 7), (2, 9), (9, 2), (5, 9), (9, 5),\
                              (1, 10), (10, 1), (8, 10), (10, 8), (2, 10), (10, 2), (6, 10), (10, 6),\
                              (3, 4), (4, 3), (5, 6), (6, 5), (7, 8), (8, 7), (9, 10), (10, 9)])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(m3r.network.main_graph, testgraph)

# Test that generate_network can process multiple rules
def test_Model_generate_network_four_rules_competitive(default_Model_four_rule_competitive_antagonists):
    m4r = default_Model_four_rule_competitive_antagonists
    m4r.generate_network()

    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8]) #1=A, 2=B, 3=R, 4=R*, 5=AR, 6=AR*, 7=BR, 8=BR*
    testgraph.add_edges_from([(1, 5), (5, 1), (3, 5), (5, 3), (2, 7), (7, 2), (3, 7), (7, 3),\
                              (1, 6), (6, 1), (4, 6), (6, 4), (2, 8), (8, 2), (4, 8), (8, 4),\
                              (3, 4), (4, 3), (5, 6), (6, 5), (7, 8), (8, 7)])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(m4r.network.main_graph, testgraph)

# Test that the ' converts to ' rule creates a valid shaped graph 
def test_Model_conversion(default_Model_conversion):
    dmc = default_Model_conversion
    
    # Have rule for simple conformational change - R(0) converts to R(1)
    dmc.generate_network()

    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3])
    testgraph.add_edges_from([(2, 3)])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(dmc.network.main_graph, testgraph)
    
# Test that the '' reversibly converts to ' rule creates a valid shaped graph 
def test_Model_rev_conversion(default_Model_conversion):
    dmc = default_Model_conversion
    
    # Have rule for simple conformational change - R(0) converts to R(1), edit to be reversible
    dmc.rule_list[0].rule = ' reversibly converts to '
    dmc.generate_network()

    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3])
    testgraph.add_edges_from([(2, 3), (3, 2)])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(dmc.network.main_graph, testgraph)
   
    
# Test that the '' converts in rapid equlibrium to ' rule creates a valid shaped graph 
def test_Model_RE_conversion(default_Model_conversion):
    dmc = default_Model_conversion
    
    # Have rule for simple conformational change - R(0) converts to R(1), edit to be reversible in rapid equlibrium
    dmc.rule_list[0].rule = ' converts in rapid equlibrium to '
    dmc.generate_network()

    # Create comparision graph shape
    testgraph = nx.DiGraph()
    testgraph.add_nodes_from([1, 2, 3])
    testgraph.add_edges_from([(2, 3), (3, 2)])
    
    # Compare shape of graph
    assert nx.algorithms.isomorphism.is_isomorphic(dmc.network.main_graph, testgraph)
   
   
# --------------------- Helper method tests in model.py ---------------------------    
            
    
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

# Both R(0) and R(1) should be returned as objects from rule matching "[]"
def test8_find_states_that_match_rule(model_for_matching_tests):
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
    
    mmt.rule_list[0].object_conf = [[0], [1]]
    
    # Run method
    subject_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'subject')
    object_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'object')
    
    # Check output
    assert subject_list == [] 
    with pytest.raises(AssertionError): # does not work this way, but would like to get it doing this in the future
        assert object_list == [s1, s2]
    
    
# Both R(0) and R(1) should be returned as objects from rule matching ([0], [1])
def test9_find_states_that_match_rule(model_for_matching_tests):
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
    
    # Run method
    subject_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'subject')
    object_list = mmt._find_states_that_match_rule(mmt.rule_list[0], 'object')
    
    # Check output
    assert subject_list == []
    assert object_list == [s1, s2]

# Test that we return valid state pairs with a simple signature
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
    
# Test that we return valid state pairs with a conformation signature. NOTE Signature matching is not quite sufficient
def test_find_association_pairs_with_conformation(model_for_matching_tests, default_Protein_instance, default_Drug_instance):
    mmt = model_for_matching_tests    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Modify rule to give A + R([1]) --> AR([1]), get the signature
    mmt.rule_list[0].object_conf = [[1]]
    ref_sig = mmt.rule_list[0].generate_signature_list()
    
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
    s3.internal_links = [(0,1)]
    
    s4 = bkcc.State()
    s4.required_protein_list = [dpi]
    s4.req_protein_conf_lists = [[0]]
    
    s5 = bkcc.State()
    s5.required_drug_list = [ddi]
    s5.required_protein_list = [dpi]
    s5.req_protein_conf_lists = [[0]]
    s5.internal_links = [(0,1)]
    
    # Make some state lists (not all these should actually be identified by the matching code, but should be rejected nontheless)
    test_sub_list = [s1, s2, s3, s4, s5]
    test_obj_list = [s1, s2, s3, s4, s5]
    
    # Run functions - need an additional internal link check to get the desired behavior
    signature_valid_tuples = mmt._find_association_pairs(ref_sig, test_sub_list, test_obj_list)
    valid_tuples, link_list = mmt._find_association_internal_link(mmt.rule_list[0], signature_valid_tuples)
    
    # Check if we got the expected result
    assert valid_tuples == [(s1, s2)]

# Test that we return valid state and split index pairs with a simple signature
def test_find_dissociation_pairs(default_Model_instance, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    #Setup rule for simple drug disassociation - "A disassociates from AR"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [ddi]
    r1.subject_conf = [None]
    r1.rule = ' dissociates from '
    r1.rule_object = [ddi, dpi]
    r1.object_conf = [None, []]
    r1.check_rule_traits()
    dmi.rule_list = [r1]
    
    # Get the signature
    ref_sig = r1.generate_signature_list()
    
    # Make some states
    s1 = bkcc.State() # Doesn't work - no protein
    s1.required_drug_list = [ddi]
    
    s2 = bkcc.State()  # Doesn't work - no drug
    s2.required_protein_list = [dpi]
    s2.req_protein_conf_lists = [[0,1]]
    
    s3 = bkcc.State() # Should work, only one split
    s3.required_drug_list = [ddi]
    s3.required_protein_list = [dpi]
    s3.req_protein_conf_lists = [[0,1]]
    
    s4 = bkcc.State() # Should work, two different possible leaving drugs (need link analysis to find the correct one)
    s4.required_drug_list = [ddi, ddi]
    s4.required_protein_list = [dpi]
    s4.req_protein_conf_lists = [[0,1]]
    
    # Make some state lists
    test_obj_list = [s1, s2, s3, s4]
    
    # Run function
    valid_tuples = dmi._find_dissociation_pairs(ref_sig, test_obj_list)

    # Check if we got the expected result
    assert valid_tuples == [(s3, (0,)), (s4, (0,)), (s4, (1,))]
    
# Test that we return valid state and split index pairs with a simple signature
def test_find_dissociation_pairs_with_conformations(default_Model_instance, default_Protein_instance, default_Drug_instance):
    dmi = default_Model_instance
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    #Setup rule for simple drug disassociation - "A disassociates from AR"
    r1 = bkcc.Rule(dmi)
    r1.rule_subject = [ddi]
    r1.subject_conf = [None]
    r1.rule = ' dissociates from '
    r1.rule_object = [ddi, dpi]
    r1.object_conf = [None, [0]]
    r1.check_rule_traits()
    dmi.rule_list = [r1]
    
    # Get the signature
    ref_sig = r1.generate_signature_list()
    
    # Make some states
    s1 = bkcc.State() # Doesn't work - no protein
    s1.required_drug_list = [ddi]
    
    s2 = bkcc.State() # Doesn't work - no drug
    s2.required_protein_list = [dpi]
    s2.req_protein_conf_lists = [[0]]
    
    s3 = bkcc.State() # Should work, only one split
    s3.required_drug_list = [ddi]
    s3.required_protein_list = [dpi]
    s3.req_protein_conf_lists = [[0]]
    
    s4 = bkcc.State() # Should work, two different possible leaving drugs (need link analysis to find the correct one)
    s4.required_drug_list = [ddi, ddi]
    s4.required_protein_list = [dpi]
    s4.req_protein_conf_lists = [[0]]
    
    s5 = bkcc.State() # Doesn't work, wrong conformation
    s5.required_drug_list = [ddi] 
    s5.required_protein_list = [dpi]
    s5.req_protein_conf_lists = [[1]]
    
    s6 = bkcc.State() # Should work, get one good split of just the drug, and one split of drug plus R(1), which we have to clean with link checking
    s6.required_drug_list = [ddi]
    s6.required_protein_list = [dpi, dpi]
    s6.req_protein_conf_lists = [[0], [1]]
    
    s7 = bkcc.State() # Should work, get two good splits of the different drugs (need link analysis to find the correct ones)
    s7.required_drug_list = [ddi, ddi]
    s7.required_protein_list = [dpi, dpi]
    s7.req_protein_conf_lists = [[0], [0]]
    
    # Make some state lists
    test_obj_list = [s1, s2, s3, s4, s5, s6, s7]
    
    # Run function
    valid_tuples = dmi._find_dissociation_pairs(ref_sig, test_obj_list)
    
    # Check if we got the expected result
    assert valid_tuples == [(s3, (0,)), (s4, (0,)), (s4, (1,)), (s6, (0,)), (s6, (0, 2)), (s7, (0,)), (s7, (1,))]

# Test for correct translation of link lists
def test_combine_internal_link_lists(model_for_matching_tests, default_Protein_instance, default_Drug_instance):
    mmt = model_for_matching_tests    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi]
    s1.req_protein_conf_lists = [[1]]
    s1.internal_links = [(0, 1)] # {AR(1)}
    
    s2 = bkcc.State()
    s2.required_drug_list = [ddi]
    s2.required_protein_list = [dpi, dpi]
    s2.req_protein_conf_lists = [[0], [1]]
    s2.internal_links = [(0, 2)] # {AR(1)}R(0)
    
    # Run method
    test_link_list = mmt._combine_internal_link_lists(s1, s2)
    
    # Check if we got the expected result
    assert test_link_list == [(0, 2), (1, 4)]
    
# Test for correct translation of link lists with a nested link
def test_combine_internal_link_lists_nested_link(model_for_matching_tests, default_Protein_instance, default_Drug_instance):
    mmt = model_for_matching_tests    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi]
    s1.req_protein_conf_lists = [[1]]
    s1.internal_links = [(0, 1)] # {AR(1)}
    
    s2 = bkcc.State()
    s2.required_drug_list = [ddi]
    s2.required_protein_list = [dpi, dpi]
    s2.req_protein_conf_lists = [[0], [1]]
    s2.internal_links = [(0, (1,2))] # {A{R(1)R(0)}} A is associated specifically with the dimer of R's.
    
    # Run method
    test_link_list = mmt._combine_internal_link_lists(s1, s2)
    
    # Check if we got the expected result
    assert test_link_list == [(0, 2), (1, (3, 4))]

# Test if the method returns the translation dictionarys
def test_combine_internal_link_lists_with_dicts(model_for_matching_tests, default_Protein_instance, default_Drug_instance):
    mmt = model_for_matching_tests    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi]
    s1.req_protein_conf_lists = [[1]]
    s1.internal_links = [(0, 1)] # {AR(1)}
    
    s2 = bkcc.State()
    s2.required_drug_list = [ddi]
    s2.required_protein_list = [dpi, dpi]
    s2.req_protein_conf_lists = [[0], [1]]
    s2.internal_links = [(0, 2)] # {AR(1)}R(0)
    
    # Run method
    test_link_list, translation_dicts = mmt._combine_internal_link_lists(s1, s2, True)
    
    # Dictionaries are returned in order: 1 to 12, 2 to 12, 12 to 1, 12 to 2
    assert len(translation_dicts) == 4
    assert translation_dicts[0][0] == 0
    assert translation_dicts[0][1] == 2
    assert translation_dicts[1][0] == 1
    assert translation_dicts[1][1] == 3    
    assert translation_dicts[1][2] == 4
    assert translation_dicts[2][0] == 0
    assert translation_dicts[2][2] == 1
    with pytest.raises(KeyError): # These keys should not exist
        translation_dicts[2][1]
        translation_dicts[2][3]
        translation_dicts[2][4]
    assert translation_dicts[3][1] == 0
    assert translation_dicts[3][3] == 1    
    assert translation_dicts[3][4] == 2
    with pytest.raises(KeyError): # These keys should not exist
        translation_dicts[3][0]
        translation_dicts[3][2]

# Test if the internal link finding/testing function returns the expected result for a simple case       
def test_find_association_internal_link(model_for_matching_tests, default_Protein_instance, default_Drug_instance):
    mmt = model_for_matching_tests    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Rule for drug association - "A associates with R([])"
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = []
    s1.internal_links = []
    
    s2 = bkcc.State()
    s2.required_drug_list = []
    s2.required_protein_list = [dpi]
    s2.req_protein_conf_lists = [[1]]
    s2.internal_links = []
    
    # Run function
    test_state_tuples, test_link_list = mmt._find_association_internal_link(mmt.rule_list[0], [(s1, s2)])
    
    # Compare outputs
    assert test_state_tuples == [(s1, s2)]
    assert test_link_list == [(0, 1)]
    
# Test if the internal link finding/testing function returns the expected result for a dimer case       
def test_find_association_internal_link_dimers(model_for_matching_tests, default_Protein_instance, default_Drug_instance):
    mmt = model_for_matching_tests    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Rule for drug association - "A associates with R([])"
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = []
    s1.internal_links = []
    
    s2 = bkcc.State()
    s2.required_drug_list = []
    s2.required_protein_list = [dpi, dpi]
    s2.req_protein_conf_lists = [[1], [0]]
    s2.internal_links = [(0,1)]
        
    # Run function
    test_state_tuples, test_link_list = mmt._find_association_internal_link(mmt.rule_list[0], [(s1, s2)])
    
    # Compare outputs
    assert test_state_tuples == [(s1, s2), (s1, s2)]
    assert test_link_list == [(0, 1), (0, 2)]
    
# Test if the internal link finding/testing function returns the expected result for a complex dimer case       
def test_find_association_internal_link_dimer_already_bound(model_for_matching_tests, default_Protein_instance, default_Drug_instance):
    mmt = model_for_matching_tests    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Edit rule for drug association - "A associates with R([1])"
    mmt.rule_list[0].object_conf = [[1]]
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = []
    s1.internal_links = []
    
    s2 = bkcc.State()
    s2.required_drug_list = [ddi]
    s2.required_protein_list = [dpi, dpi]
    s2.req_protein_conf_lists = [[1], [0]]
    s2.internal_links = [(1, 2), (0, 1)]
        
    # Run function
    test_state_tuples, test_link_list = mmt._find_association_internal_link(mmt.rule_list[0], [(s1, s2)])
    
    # Compare outputs
    assert test_state_tuples == []
    assert test_link_list == []

# Test if the internal link finding/testing function returns the expected result for a complex dimer case       
def test_find_association_internal_link_dimer_open_for_binding(model_for_matching_tests, default_Protein_instance, default_Drug_instance):
    mmt = model_for_matching_tests    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Edit rule for drug association - "A associates with R([1])"
    mmt.rule_list[0].object_conf = [[1]]
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = []
    s1.req_protein_conf_lists = []
    s1.internal_links = []
    
    s2 = bkcc.State()
    s2.required_drug_list = [ddi]
    s2.required_protein_list = [dpi, dpi]
    s2.req_protein_conf_lists = [[1], [1]]
    s2.internal_links = [(1, 2), (0, 1)]
        
    # Run function
    test_state_tuples, test_link_list = mmt._find_association_internal_link(mmt.rule_list[0], [(s1, s2)])
    
    # Compare outputs
    assert test_state_tuples == [(s1, s2)]
    assert test_link_list == [(0,3)]

# Test for correct translation of link lists
def test_split_internal_link_list(model_for_dissociation_matching, default_Protein_instance, default_Drug_instance):
    mdm = model_for_dissociation_matching    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi, dpi]
    s1.req_protein_conf_lists = [[0, 1], [0, 1]]
    s1.internal_links = [(0, 1), (1, 2)]
    
    split_indices = [0]
    
    # Run method
    broken_state12_links, state1_links, state2_links = mdm._split_internal_link_list(s1, split_indices)
    
    # Check if we got the expected result
    assert broken_state12_links == [(0, 1)]
    assert state1_links == []
    assert state2_links == [(0, 1)]

# Test for correct creation of translation dictionaries
def test_split_internal_link_list_translation_dictionaries(model_for_dissociation_matching, default_Protein_instance, default_Drug_instance):
    mdm = model_for_dissociation_matching    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Make some states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi, dpi]
    s1.req_protein_conf_lists = [[0, 1], [0, 1]]
    s1.internal_links = [(0, 1), (1, 2)]
    
    split_indices = [0]
    
    # Run method
    broken_state12_links, state1_links, state2_links, translation_dicts = mdm._split_internal_link_list(s1, split_indices, True)

    # Dictionaries are returned in order: 1 to 12, 2 to 12, 12 to 1, 12 to 2; 1 = subject, 12 = object, 2 = third state
    assert len(translation_dicts) == 4
    assert translation_dicts[0][0] == 0
    assert translation_dicts[1][0] == 1
    assert translation_dicts[1][1] == 2
    assert translation_dicts[2][0] == 0    
    assert translation_dicts[3][1] == 0
    assert translation_dicts[3][2] == 1
    with pytest.raises(KeyError): # These keys should not exist
        translation_dicts[0][1]
        translation_dicts[0][2]
        translation_dicts[1][2]
        translation_dicts[2][1]
        translation_dicts[2][2]
        translation_dicts[3][0]
    
# Test if the internal link finding/testing function returns the expected result for a dissociation case       
def test_find_dissociation_internal_link(model_for_dissociation_matching, default_Protein_instance, default_Drug_instance):
    mdm = model_for_dissociation_matching    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Rule for drug association - "A dissociates from AR([])"
    
    # Make a states
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi, dpi]
    s1.req_protein_conf_lists = [[0, 1], [0, 1]]
    s1.internal_links = [(0, 1), (1, 2)]
    
    # Run function
    test_state_split_list, test_link_lists = mdm._find_dissociation_internal_link(mdm.rule_list[0], [(s1, (0,)), (s1, (0, 1))])
    state1_links, state2_links = test_link_lists[0]
    
    # Compare outputs, should only allow the drug to dissociate, not the dimer
    assert test_state_split_list == [(s1, (0,))]
    assert state1_links == []
    assert state2_links == [(0, 1)] # In state2 index numbers
    
# Test if the internal link finding/testing function returns the expected result for a simple conversion case       
def test_find_conversion_simple_internal_link(default_Model_conversion, default_Protein_instance, default_Drug_instance):
    dmc = default_Model_conversion   
    dpi = default_Protein_instance     
    ddi = default_Drug_instance # For typing convenience
    
    # We have rule R([0]) converts to R([1])
    
    # Make a state
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi, dpi]
    s1.req_protein_conf_lists = [[0], [0]]
    s1.internal_links = [(0, 1)]
    
    # Run function
    #(current_sub, current_indices, convert_rule_components, current_convert_conf)
    possible_tuples = [(s1, (1,), [dpi], [[1]]), (s1, (2,), [dpi], [[1]])]
    test_results = dmc._find_conversion_internal_link(dmc.rule_list[0], possible_tuples)

    # Test output (matching_subject_state, new_component_list, new_conformation_list, new_link_tuples)
    # Both R components should be able to be converted - NOTE: conversion of both required re-running the rule on the new graph
    assert len(test_results) == 2
    assert test_results[0][0] == s1
    assert test_results[0][1] == [ddi, dpi, dpi]
    assert test_results[0][2] == [None, [1], [0]]
    assert test_results[0][3] == [(0, 1)]
    assert test_results[1][0] == s1
    assert test_results[1][1] == [ddi, dpi, dpi]
    assert test_results[1][2] == [None, [0], [1]]
    assert test_results[1][3] == [(0, 1)]
    
# Test if the internal link finding/testing function returns the expected result for a more complex conversion case       
def test_find_conversion_complex_internal_link(default_Model_conversion, default_Protein_instance, default_Drug_instance):
    dmc = default_Model_conversion   
    dpi = default_Protein_instance     
    ddi = default_Drug_instance # For typing convenience
    r1 = dmc.rule_list[0]
    
    # We have rule R([0]) converts to R([1]), change to require binding of A
    r1.rule_subject = [ddi, dpi]
    r1.subject_conf = [None, [0]]
    r1.rule_object = [ddi, dpi]
    r1.object_conf = [None, [1]]
    r1.check_rule_traits()

    # Make a state
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi, dpi]
    s1.req_protein_conf_lists = [[0], [0]]
    s1.internal_links = [(0, 1)]
    
    # Run function
    #(current_sub, current_indices, convert_rule_components, current_convert_conf)
    possible_tuples = [(s1, (0, 1), [ddi, dpi], [None, [1]]), (s1, (0, 2), [ddi, dpi], [None, [1]])]
    test_results = dmc._find_conversion_internal_link(dmc.rule_list[0], possible_tuples)

    # Test output (matching_subject_state, new_component_list, new_conformation_list, new_link_tuples)
    # Only the R bound to A should be allowed to convert
    assert len(test_results) == 1
    assert test_results[0][0] == s1
    assert test_results[0][1] == [ddi, dpi, dpi]
    assert test_results[0][2] == [None, [1], [0]]
    assert test_results[0][3] == [(0, 1)]
    
 # Tests if the broken link is consistant with a dissociation rule
def test_compare_components_dissociation_rule_and_link(model_for_dissociation_matching, default_Protein_instance, default_Drug_instance):
    mdm = model_for_dissociation_matching    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Rule for drug association - "A dissociates from AR([])"
    
    # Make a state
    s1 = bkcc.State()
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi, dpi]
    s1.req_protein_conf_lists = [[0, 1], [0, 1]]
    s1.internal_links = [(0, 1), (1, 2)]
    
    # Try to break the AR link - should work
    assert mdm._compare_components_dissociation_rule_and_link(mdm.rule_list[0], (0, 1), s1)
    
    # Try to break the AR link written backwards - should work
    assert mdm._compare_components_dissociation_rule_and_link(mdm.rule_list[0], (1, 0), s1)
    
    # Try to break the RR link - should not work
    assert not mdm._compare_components_dissociation_rule_and_link(mdm.rule_list[0], (1, 2), s1)
    
def test_collect_link_components(model_for_dissociation_matching, default_Protein_instance, default_Drug_instance):
    mdm = model_for_dissociation_matching    
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    
    # Create a reference list
    ref_comp_list = [ddi, ddi, dpi, dpi, dpi]
    ref_conf_list = [None, None, [0], [1], [0, 1]]
    
    # Run fuction with a single index
    test_comp, test_conf = mdm._collect_link_components(4, ref_comp_list, ref_conf_list)
    assert test_comp == [dpi]
    assert test_conf == [[0, 1]]
    
    # Run fuction with a double index
    test_comp, test_conf = mdm._collect_link_components((1, 3), ref_comp_list, ref_conf_list)
    assert test_comp == [ddi, dpi]
    assert test_conf == [None, [1]]
    
    # Run fuction with a complex index
    test_comp, test_conf = mdm._collect_link_components((1, (0, 2), 3), ref_comp_list, ref_conf_list)
    assert test_comp == [ddi, ddi, dpi, dpi]
    assert test_conf == [None, None, [0], [1]]
    
    # Run fuction with a invalid argument
    with pytest.raises(ValueError):
        test_comp, test_conf = mdm._collect_link_components((1, 1.1), ref_comp_list, ref_conf_list)
        
# Test that we return valid state and conversion index pairs with a simple signature
def test_find_conversion_pairs(default_Model_conversion, default_Protein_instance, default_Drug_instance):
    dmc = default_Model_conversion
    dpi = default_Protein_instance
    ddi = default_Drug_instance # For typing convenience
    r1 = dmc.rule_list[0]
    
    # We have rule R([0]) converts to R([1])
    ref_sig = r1.generate_signature_list()
    
    # Make some states
    s1 = bkcc.State() # Doesn't work - no protein
    s1.required_drug_list = [ddi]
    
    s2 = bkcc.State()  # Doesn't work - wrong conformation
    s2.required_protein_list = [dpi]
    s2.req_protein_conf_lists = [[0,1]]
    
    s3 = bkcc.State() # Doesn't work - no conformation change
    s3.required_protein_list = [dpi]
    s3.req_protein_conf_lists = [[1]]
    
    s4 = bkcc.State() # Should work
    s4.required_protein_list = [dpi]
    s4.req_protein_conf_lists = [[0]]
    
    # Make some state lists
    test_obj_list = [s1, s2, s3, s4]
    
    # Run function
    test_results = dmc._find_conversion_pairs(r1, ref_sig, test_obj_list)

    # Check if we got the expected result
    assert len(test_results) == 1
    assert test_results[0][0] == s4
    assert test_results[0][1] == (0,)
    assert test_results[0][2] == [dpi]
    assert test_results[0][3] == [[1]]
     
# Test that we return states with possible competitive states and the possible indices that match the rule
def test_find_competitive_states(default_Model_competition, default_Drug_instance, second_Drug_instance, default_Protein_instance):
    dmc = default_Model_competition
    ddi = default_Drug_instance 
    sdi = second_Drug_instance
    dpi = default_Protein_instance# For typing convenience
    r1 = dmc.rule_list[0]
    
    # We have rule A is competitive with B
    ref_sig = r1.generate_signature_list()
    
    # Make some states
    s1 = bkcc.State() # OK, only one ligand
    s1.required_drug_list = [ddi]
    s1.required_protein_list = [dpi]
    s1.req_protein_conf_lists = [[0]]
    s1.internal_links = [(0, 1)]
    
    s2 = bkcc.State()  # Not allowed, two ligands bound to one receptor
    s2.required_drug_list = [ddi, sdi]
    s2.required_protein_list = [dpi]
    s2.req_protein_conf_lists = [[0]]
    s2.internal_links = [(0, 2), (1, 2)]
    
    s3 = bkcc.State()  # OK, two ligands bound to different receptors
    s3.required_drug_list = [ddi, sdi]
    s3.required_protein_list = [dpi, dpi]
    s3.req_protein_conf_lists = [[0], [0]]
    s3.internal_links = [(0, 2), (1, 3)]
    
    s4 = bkcc.State()  # Not allowed, three ligands, two recpetors, with both drugs bound to same receptor
    s4.required_drug_list = [ddi, ddi, sdi]
    s4.required_protein_list = [dpi, dpi]
    s4.req_protein_conf_lists = [[0], [0]]
    s4.internal_links = [(0, 3), (1, 4), (2, 3)]
      
    # Test states
    test_state_list = [s1, s2, s3, s4]
    test_results = dmc._find_competitive_states(r1, ref_sig, test_state_list)
    
    # Check if we got the expected results, tuple= (state, (indices of subject), (indices of object))
    assert len(test_results) == 4
    assert test_results[0][0] == s2
    assert test_results[0][1] == (0,)
    assert test_results[0][2] == (1,)
    assert test_results[1][0] == s3 # Will be allowed by link-checking
    assert test_results[1][1] == (0,)
    assert test_results[1][2] == (1,)
    assert test_results[2][0] == s4 
    assert test_results[2][1] == (0,)
    assert test_results[2][2] == (2,)
    assert test_results[3][0] == s4 # Will be allowed by link-checking, but overall state disqualified by previous configuration
    assert test_results[3][1] == (1,)
    assert test_results[3][2] == (2,)

    
# ------Tests for graph edge objects------


# Test that association graphs have the right number of edge objects
def test_Irr_Association_edges(model_for_matching_tests):
    mmt = model_for_matching_tests
    
    # Have rule A associates with R([]), generate corrosponding network
    mmt.generate_network()
    
    # Analysis: We should have two association-ST objects assigned in the network
    
    # Get edge info
    test_tails, test_heads, test_STobjs = zip(*mmt.network.main_graph.edges(data = 'reaction_type'))

    # Count StateTransition objects
    ST_counts = collections.Counter(test_STobjs)
    
    # Find Association objects
    association_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.Association)]
        
    # See if we get the expected number of objects and connections
    assert len(association_STobjs) == 2
    for ST_key in association_STobjs:
        assert ST_counts[ST_key] == 2 # Associations require placement on 2 edges each
        
# Test that association graphs have the right number of edge objects
def test_Rev_Association_edges(model_for_matching_tests):
    mmt = model_for_matching_tests
    
    # Have rule A associates with R([]), change to reversible association and generate corrosponding network
    mmt.rule_list[0].rule = ' reversibly associates with '
    mmt.generate_network()
    
    # Analysis: We should have two association-ST objects and dissociation-ST objects assigned in the network
    
    # Get edge info
    test_tails, test_heads, test_STobjs = zip(*mmt.network.main_graph.edges(data = 'reaction_type'))

    # Count StateTransition objects
    ST_counts = collections.Counter(test_STobjs)
    
    # Find Association objects
    association_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.Association)]
        
    # See if we get the expected number of objects and connections
    assert len(association_STobjs) == 2
    for ST_key in association_STobjs:
        assert ST_counts[ST_key] == 2 # Associations require placement on 2 edges each
        
    # Find Dissociation objects
    dissociation_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.Dissociation)]
        
    # See if we get the expected number of objects and connections
    assert len(dissociation_STobjs) == 2
    for ST_key in dissociation_STobjs:
        assert ST_counts[ST_key] == 2 # Dissociations require placement on 2 edges each
        
# Test that dissociation graphs have the right number of edge objects
def test_Irr_Dissociation_edges(model_for_dissociation_matching):
    mdm = model_for_dissociation_matching
    
    # Setup test network, must have starting associated state
    s1 = bkcc.State()
    s1.required_drug_list = [mdm.drug_list[0]]
    s1.required_protein_list = [mdm.protein_list[0]]
    s1.req_protein_conf_lists = [[0, 1]]
    s1.internal_links = [(0, 1)]
    mdm.network.main_graph.add_node(s1)
    
    # Apply the existing association rule - "A dissociates from AR([0,1])"
    mdm.apply_rules_to_network()

    # Analysis: We should have one dissociation-ST objects assigned in the network
    
    # Get edge info
    print(*zip(*mdm.network.main_graph.edges(data = 'reaction_type')))
    test_tails, test_heads, test_STobjs = zip(*mdm.network.main_graph.edges(data = 'reaction_type'))

    # Count StateTransition objects
    ST_counts = collections.Counter(test_STobjs)
    
    # Find Dissociation objects
    dissociation_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.Dissociation)]
        
    # See if we get the expected number of objects and connections
    assert len(dissociation_STobjs) == 1
    for ST_key in dissociation_STobjs:
        assert ST_counts[ST_key] == 2 # Dissociations require placement on 2 edges each
        
# Test that dissociation graphs have the right number of edge objects
def test_Rev_Dissociation_edges(model_for_dissociation_matching):
    mdm = model_for_dissociation_matching
    
    # Have rule A dissociates from R([0,1])A, change to reversible dissociation and generate corrosponding network
    mdm.rule_list[0].rule = ' reversibly dissociates from '
    
    # Setup test network, must have starting associated state
    s1 = bkcc.State()
    s1.required_drug_list = [mdm.drug_list[0]]
    s1.required_protein_list = [mdm.protein_list[0]]
    s1.req_protein_conf_lists = [[0, 1]]
    s1.internal_links = [(0, 1)]
    mdm.network.main_graph.add_node(s1)
    
    # Apply the existing association rule - "A revesibly dissociates from AR([0,1])"
    mdm.apply_rules_to_network()
    
    # Analysis: We should have one dissociation-ST object and one association-ST object assigned in the network
    
    # Get edge info
    test_tails, test_heads, test_STobjs = zip(*mdm.network.main_graph.edges(data = 'reaction_type'))

    # Count StateTransition objects
    ST_counts = collections.Counter(test_STobjs)
    
    # Find Dissociation objects
    dissociation_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.Dissociation)]
        
    # See if we get the expected number of objects and connections
    assert len(dissociation_STobjs) == 1
    for ST_key in dissociation_STobjs:
        assert ST_counts[ST_key] == 2 # Associations require placement on 2 edges each
    
    # Find Association objects
    association_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.Association)]
        
    # See if we get the expected number of objects and connections
    assert len(association_STobjs) == 1
    for ST_key in association_STobjs:
        assert ST_counts[ST_key] == 2 # Associations require placement on 2 edges each

# Test that conversion graphs have the right number of edge objects
def test_Irr_Conversion(default_Model_conversion):
    dmc = default_Model_conversion
    
    # Have rule R(0) converts to R(1), generate corrosponding network
    dmc.generate_network()
    
    # Analysis: We should have one conversion-ST object assigned in the network
    
    # Get edge info
    test_tails, test_heads, test_STobjs = zip(*dmc.network.main_graph.edges(data = 'reaction_type'))

    # Count StateTransition objects
    ST_counts = collections.Counter(test_STobjs)
    
    # Find Conversion objects
    conversion_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.Conversion)]
        
    # See if we get the expected number of objects and connections
    assert len(conversion_STobjs) == 1
    assert conversion_STobjs[0].reference_direction == True # A single converstion should always be the reference direction
    for ST_key in conversion_STobjs:
        assert ST_counts[ST_key] == 1 # Each Conversion reaction requires only 1 edge
        
# Test that conversion graphs have the right number of edge objects
def test_Rev_Conversion(default_Model_conversion):
    dmc = default_Model_conversion
    
    # Have rule R(0) converts to R(1), change to reversible conversion and generate corrosponding network
    dmc.rule_list[0].rule = ' reversibly converts to '
    dmc.generate_network()
    
    # Analysis: We should have two conversion-ST object assigned in the network
    
    # Get edge info
    test_tails, test_heads, test_STobjs = zip(*dmc.network.main_graph.edges(data = 'reaction_type'))

    # Count StateTransition objects
    ST_counts = collections.Counter(test_STobjs)
    
    # Find Conversion objects
    conversion_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.Conversion)]
        
    # See if we get the expected number of objects and connections
    assert len(conversion_STobjs) == 2
    assert conversion_STobjs[0].reference_direction ^ conversion_STobjs[1].reference_direction # XOR (these values are booleans), only one direction should be reference
    for ST_key in conversion_STobjs:
        assert ST_counts[ST_key] == 1 # Each Conversion reaction requires only 1 edge

# Test that association graphs have the right number of edge objects
def test_RE_Association_edges(model_for_matching_tests):
    mmt = model_for_matching_tests
    
    # Have rule A associates with R([]), change to rapid equlibrium association and generate corrosponding network
    mmt.rule_list[0].rule = ' associates and dissociates in rapid equlibrium with '
    mmt.generate_network()
    
    # Analysis: We should have two RE association-ST objects and RE dissociation-ST objects assigned in the network
    
    # Get edge info
    test_tails, test_heads, test_STobjs = zip(*mmt.network.main_graph.edges(data = 'reaction_type'))

    # Count StateTransition objects
    ST_counts = collections.Counter(test_STobjs)
    
    # Find Association objects
    association_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.RE_Association)]
        
    # See if we get the expected number of objects and connections
    assert len(association_STobjs) == 2
    for ST_key in association_STobjs:
        assert ST_counts[ST_key] == 2 # Associations require placement on 2 edges each
        
    # Find Dissociation objects
    dissociation_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.RE_Dissociation)]
        
    # See if we get the expected number of objects and connections
    assert len(dissociation_STobjs) == 2
    for ST_key in dissociation_STobjs:
        assert ST_counts[ST_key] == 2 # Dissociations require placement on 2 edges each
        
# Test that RE dissociation graphs have the right number of edge objects
def test_RE_Dissociation_edges(model_for_dissociation_matching):
    mdm = model_for_dissociation_matching
    
    # Have rule A dissociates from R([0,1])A, change to rapid equlibrium dissociation and generate corrosponding network
    mdm.rule_list[0].rule = ' dissociates and reassociates in rapid equlibrium from '
    
    # Setup test network, must have starting associated state
    s1 = bkcc.State()
    s1.required_drug_list = [mdm.drug_list[0]]
    s1.required_protein_list = [mdm.protein_list[0]]
    s1.req_protein_conf_lists = [[0, 1]]
    s1.internal_links = [(0, 1)]
    mdm.network.main_graph.add_node(s1)
    
    # Apply the existing association rule - "A revesibly dissociates and reassociates in rapid equlibrium from AR([0,1])"
    mdm.apply_rules_to_network()
    
    # Analysis: We should have one RE dissociation-ST object and one RE association-ST object assigned in the network
    
    # Get edge info
    test_tails, test_heads, test_STobjs = zip(*mdm.network.main_graph.edges(data = 'reaction_type'))

    # Count StateTransition objects
    ST_counts = collections.Counter(test_STobjs)
    
    # Find Dissociation objects
    dissociation_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.RE_Dissociation)]
        
    # See if we get the expected number of objects and connections
    assert len(dissociation_STobjs) == 1
    for ST_key in dissociation_STobjs:
        assert ST_counts[ST_key] == 2 # Associations require placement on 2 edges each
    
    # Find Association objects
    association_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.RE_Association)]
        
    # See if we get the expected number of objects and connections
    assert len(association_STobjs) == 1
    for ST_key in association_STobjs:
        assert ST_counts[ST_key] == 2 # Associations require placement on 2 edges each


# Test that RE conversion graphs have the right number of edge objects
def test_RE_Conversion(default_Model_conversion):
    dmc = default_Model_conversion
    
    # Have rule R(0) converts to R(1), change to Rapid equlibrium conversion and generate corrosponding network
    dmc.rule_list[0].rule = ' converts in rapid equlibrium to '
    dmc.generate_network()
    
    # Analysis: We should have two RE conversion-ST object assigned in the network
    
    # Get edge info
    test_tails, test_heads, test_STobjs = zip(*dmc.network.main_graph.edges(data = 'reaction_type'))

    # Count StateTransition objects
    ST_counts = collections.Counter(test_STobjs)
    
    # Find RE Conversion objects
    conversion_STobjs = [x for x in ST_counts.keys() if isinstance(x, bkcc.RE_Conversion)]
        
    # See if we get the expected number of objects and connections
    assert len(conversion_STobjs) == 2
    assert conversion_STobjs[0].reference_direction ^ conversion_STobjs[1].reference_direction # XOR (these values are booleans), only one direction should be reference
    for ST_key in conversion_STobjs:
        assert ST_counts[ST_key] == 1 # Each Conversion reaction requires only 1 edge