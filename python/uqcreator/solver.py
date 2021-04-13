from numba import cfunc
from ctypes import CDLL, Structure
from ctypes.util import find_library

class Solver():
    def __init__(self, settings, problem, closure):
        self.libpath = find_library("UQCreator")
        self.settings = settings
        self.problem = problem
        self.closure = closure        
        
    def run(self):
        self.lib = CDLL(self.libpath)
        G = cfunc("array(double, 2d, C)(array(double, 2d, C), array(double, 2d, C), array(double, 1d, C), array(double, 1d, C), uintc)")(self.problem.G)
        F = cfunc("array(double, 2d, C)(array(double, 1d, C))")(self.problem.F)
        computeDt = cfunc("double(uintc)")(self.problem.ComputeDt)
        IC = cfunc("array(double, 1d, C)(array(double, 1d, C), array(double, 1d, C))")(self.problem.IC)
        U = cfunc("void(array(double, 3d, C), array(double, 3d, C))")(self.closure.U)
        DU = cfunc("void(array(double, 2d, C), array(double, 1d, C))")(self.closure.DU)
        SolveClosure = cfunc("void(array(double, 3d, C), array(double, 3d, C), uintc)")(self.closure.DS)
        self.settings.validate()        
        self.lib.run(
            self.settings.convert_to_structure(),            
            G, F, computeDt, IC,
            U, DU, SolveClosure
        )


