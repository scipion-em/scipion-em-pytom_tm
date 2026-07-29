from pyworkflow.object import String
from tomo.objects import SetOfTomograms, Tomogram


class PytomScoreTomogram(Tomogram):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._tomoFile = String() # Path to the corresponding tomogram

class SetOfPytomScoreTomograms(SetOfTomograms):
    ITEM_TYPE = PytomScoreTomogram
    EXPOSE_ITEMS = True