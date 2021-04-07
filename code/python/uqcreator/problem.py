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