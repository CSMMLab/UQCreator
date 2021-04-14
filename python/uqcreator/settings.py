from ctypes import Structure, c_uint, c_double, c_char_p, POINTER

class CSettingsStruct(Structure):
    _fields_ = [
        ("cwd", c_char_p),
        ("dim", c_uint),
        ("outputDir", c_char_p),
        ("outputFile", c_char_p),
        ("referenceFile", c_char_p),
        ("writeFrequency", c_uint),
        ("timesteppingType", c_char_p),
        ("nUncertainties", c_uint),
        ("dist", POINTER(c_char_p)),
        ("sigma", POINTER(c_double)),
        ("cfl", c_double),
        ("tEnd", c_double),
        ("residual", c_double),
        ("filter", c_char_p),
        ("nMultiElements", c_uint),
        ("nRefinementLevels", c_uint),
        ("moments", POINTER(c_uint)),
        ("momentDegreeType", POINTER(c_char_p)),
        ("quadOrder", POINTER(c_uint)),
        ("quadType", POINTER(c_char_p)),
        ("nRetardationSteps", c_uint),
        ("retardationSteps", POINTER(c_uint)),
        ("retardationResidual", POINTER(c_double)),
        ("refinementThresholds", POINTER(c_double)),
        ("regularizationStrength", c_double),
        ("filterStrength", c_double),
        ("maxIterations", c_uint),
        ("epsilon", c_double)
    ]


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