import numpy as np
from uqcreator.closures.closure import Closure

class StochasticGalerkin(Closure):
    def __init__(self, nStates):
        self.nStates = nStates

    def U(self, out, Lambda):
        out = Lambda
    
    def DU(self, y, Lambda):
        y = np.eye(self.nStates)

    def SolveClosure(self, u, Lambda, level):
        Lambda = u