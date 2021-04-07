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