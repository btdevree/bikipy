"""Custom exceptions for BiKiPy"""

# Base custom exception
class BikipyException(Exception):
    pass # No special activity needed yet  

class ComponentNotValidError(BikipyException):
    """Exception for when a component object is not in a valid configuration."""
    
    # Hand message to the base exception module
    def __init__(self, message = ''):
        super().__init__(message)

class RuleNotValidError(BikipyException):
    """Exception for when a rule object is not in a valid configuration."""
    
    # Hand message to the base exception module
    def __init__(self, message = ''):
        super().__init__(message)

