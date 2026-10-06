from os.path import join
from typing import Optional, Union

from pwem.objects import Volume
from pwem.protocols import EMProtocol
from pytom_tm.constants import MRC_EXT
from pyworkflow.object import Pointer
from pyworkflow.utils import makePath
from tomo.objects import SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries


class ProtPytomBase(EMProtocol):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    #---------------------Utils functions---------------------

    def _getFormAttrib(self, attribName: str, returnPointer: bool = False) -> Optional[Union[SetOfTiltSeries,
    SetOfTomograms, SetOfCTFTomoSeries, Volume, Pointer]]:
        inTsPointer = getattr(self, attribName, None)
        if not inTsPointer:
            return None
        else:
            return inTsPointer if returnPointer else inTsPointer.get()

    def _getCurrentTomoTmpDir(self, tsId: str) -> str:
        return self._getTmpPath(tsId)

    def _getConvertedOrLinkedName(self, tsId: str, suffix: str = '') -> str:
        return join(self._getCurrentTomoDir(tsId), f'{tsId}{suffix}{MRC_EXT}')

    def _getCurrentTomoDir(self, tsId: str) -> str:
        return self._getExtraPath(tsId)

    def make_dirs(self, tsId: str)-> None:
        tsDir = self._getCurrentTomoDir(tsId)
        tsTmpDir = self._getCurrentTomoTmpDir(tsId)
        makePath(tsDir, tsTmpDir)
