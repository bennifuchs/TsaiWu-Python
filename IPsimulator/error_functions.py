'''
Created on 05.07.2013
defines various error functions
@author: c8441101
'''


class MaterialParameterError(Exception):
    def __init__(self):
        pass
    
    def __str__(self):
        return "unsupported model/material parameter"


class MaterialInputError(Exception):
    def __init__(self):
        pass
    
    def __str__(self):
        return "unsupported value for input parameter"

    
class DeltaLambdaError(Exception):
    def __init__(self):
        pass
    
    def __str__(self):
        return "plastic multiplier delta_lambda zero or negative"

    
class UmatError(Exception):
    def __init__(self):
        pass
    
    def __str__(self):
        return " error in Umat occured"
    
class OuterNewtonError(Exception):
    def __init__(self):
        pass
    
    def __str__(self):
        return " error in Outer Newton Iteration"

class OuterNewtonConvergenceError(OuterNewtonError):    
    def __str__(self):
        return "Outer Newton Iteration doesn't converge within the given maximum number of iterations"

class OuterNewtonSubsteppingError(OuterNewtonError):    
    def __str__(self):
        return "error in substepping in Outer Newton Iteration"