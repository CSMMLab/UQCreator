from numba import cfunc
from ctypes import CDLL
from ctypes.util import find_library
from abc import ABCMeta, abstractmethod

class Closure(metaclass = ABCMeta):
    def __init__(self):
        pass
    @abstractmethod
    def U(self):
        pass
    @abstractmethod
    def DU(self):
        pass
    @abstractmethod
    def DS(self):
        pass


class Problem(metaclass = ABCMeta):
    def __init__(self):
        pass
    @abstractmethod
    def G(self):
        pass
    @abstractmethod
    def F(self):
        pass
    @abstractmethod
    def ComputeDt(self):
        pass
    @abstractmethod
    def IC(self):
        pass


class Settings():
    def __init__(self):
        self.general = {
            "outputDir" : "",
            "referenceSolution" : "",
            "writeFrequency" : 0,
            "restartFile" : ""
        }
        self.mesh = {
            "dimension" : 2,
            "format" : "SU2",
            "SU2MeshFile" : "", 
            "SU2BC" : [[], []],
            "outputFile" : ""
        }
        self.problem = {
            "timestepping" : "explicitEuler",
            "distribution" : [[], []],
            "CFL" : 0.9,
            "tEnd" : 1.0,
            "residual" : 1e-7
        }
        self.moment_system = {
            "moments" : [[],[]],
            "quadPoints" : [[],[]],
            "maxIterations" : 1000,
            "epsilon" : 1e-7
        }

    def validate(self):
        return True

    def convert_to_structure(self):
        return True

class UQCreator():
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


