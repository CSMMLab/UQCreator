from numba import cfunc
from ctypes import CDLL, Structure
from ctypes.util import find_library
from abc import ABCMeta, abstractmethod

class Solver():
    def __init__(self, settings, problem, closure):
        self.libpath = find_library("UQCreator")
        self.settings = settings
        self.problem = problem
        self.closure = closure        
        
    def run(self):
        self.lib = CDLL(self.libpath)
        G = cfunc("int32(int32)")(self.problem.G)
        F = cfunc("int32(int32)")(self.problem.F)
        computeDt = cfunc("int32(int32)")(self.problem.ComputeDt)
        IC = cfunc("int32(int32)")(self.problem.IC)
        U = cfunc("int32(int32)")(self.closure.U)
        DU = cfunc("int32(int32)")(self.closure.DU)
        DS = cfunc("int32(int32)")(self.closure.DS)
        self.settings.validate()        
        self.lib.run(
            self.settings.convert_to_structure(),            
            G, F, computeDt, IC,
            U, DU, DS
        )


