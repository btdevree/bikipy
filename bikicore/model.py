"""Class for the model objects used throughout the program.
"""

from traits.api import HasTraits, Int

#Model class
class Model(HasTraits):
    
    #Initalize traits
    IDnum = Int()
    
    def __init__(self, IDnum):
        self.IDnum = IDnum

#Method to create a new model from scratch
def create_new_model(IDnum):
    new_model = Model(IDnum)
    print(new_model)
    return new_model