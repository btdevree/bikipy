"""Test suite for classes in model.py

"""
import pytest
import bikipy.bikicore.model as bkcm

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
    
# ---- Unit tests ----

# --Tests for Model objects--

# Test if Model object has the required properties
def test_Drug_has_number(default_Model_instance):
    assert hasattr(default_Model_instance, 'number')
    
def test_Drug_has_name(default_Model_instance):
    assert hasattr(default_Model_instance, 'name')

def test_Drug_has_parent(default_Model_instance):
    assert hasattr(default_Model_instance, 'parent_model')

def test_Drug_has_ID(default_Model_instance):
    assert hasattr(default_Model_instance, 'ID')
    
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

    
    