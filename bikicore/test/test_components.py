"""Test suite for classes in components.py

"""
import pytest
import bikipy.bikicore.components as bkcc
from bikipy.bikicore.exceptions import ComponentNotValidError

#---- Testing fixtures ----

# Create a default Drug object for reuse in tests
@pytest.fixture()
def default_Drug_instance():
    return bkcc.Drug();
    
# Create a default Drug object for reuse in tests
@pytest.fixture()
def default_Protein_instance():
    return bkcc.Protein();

# Create a default State object for reuse in tests
@pytest.fixture()
def default_State_instance():
    return bkcc.State();

# Create a default ConformationalChange object for reuse in tests
@pytest.fixture()
def default_ConfChange_instance():
    return bkcc.ConformationalChange();

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

def test_State_has_IDnumber(default_State_instance):
    assert hasattr(default_State_instance, 'IDnumber')

# --Tests for ConformationalChange objects--

# Test if ConformationalChange objects have the required properties
def test_ConfChange_has_name(default_ConfChange_instance):
    assert hasattr(default_ConfChange_instance, 'name')
    
def test_ConfChange_has_number(default_ConfChange_instance):
    assert hasattr(default_ConfChange_instance, 'number')

def test_ConfChange_has_IDnumber(default_ConfChange_instance):
    assert hasattr(default_ConfChange_instance, 'IDnumber')
    
#Test if Association objects have the required properties
def test_Association_has_name(default_Association_instance):
    assert hasattr(default_Association_instance, 'name')
    
def test_Association_has_number(default_Association_instance):
    assert hasattr(default_Association_instance, 'number')

def test_Association_has_IDnumber(default_ConfChange_instance):
    assert hasattr(default_ConfChange_instance, 'IDnumber')
    
        
