from pyworkflow.object import String
from tomo.objects import SetOfTomograms, Tomogram


class PytomScoreTomogram(Tomogram):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._tomoFile = String() # Path to the corresponding tomogram
        self._jsonFile = String()

    def getTomoFile(self)->str:
        return self._tomoFile.get()

    def setTomoFile(self, tomoFile:str)->None:
        self._tomoFile.set(tomoFile)

    def getJsonFile(self)->str:
        return self._jsonFile.get()

    def setJsonFile(self, jsonFile:str)->None:
        self._jsonFile.set(jsonFile)


class SetOfPytomScoreTomograms(SetOfTomograms):
    ITEM_TYPE = PytomScoreTomogram
    EXPOSE_ITEMS = True
