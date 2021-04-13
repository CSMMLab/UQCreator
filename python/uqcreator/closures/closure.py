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
    def SolveClosure(self):
        pass