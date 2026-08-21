from pyworkflow.object import String
from tomo.objects import SetOfTomograms, Tomogram


class PytomScoreTomogram(Tomogram):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._tomoFile = String() # Path to the corresponding tomogram

    def getTomoFile(self)->str:
        return self._tomoFile.get()

    def setTomoFile(self, tomoFile:str)->None:
        self._tomoFile.set(tomoFile)


class SetOfPytomScoreTomograms(SetOfTomograms):
    ITEM_TYPE = PytomScoreTomogram
    EXPOSE_ITEMS = True