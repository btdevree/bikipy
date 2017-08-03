"""Custom exceptions for BiKiPy"""

# Base custom exception
class BikipyException(Exception):
    pass   

class ComponentNotValidError(BikipyException):

    def __init__(self, message = ''):
        self.message = message
    