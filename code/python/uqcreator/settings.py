from ctypes import Structure, c_int, c_double

class CSettingsStruct(Structure):
    _fields_ = [("x", c_int),
                ("y", c_int)]


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