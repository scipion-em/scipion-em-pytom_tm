# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team
# *
# * National Center of Biotechnology, CSIC, Spain
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from os.path import exists
from typing import Union, List, Tuple
from cistem.protocols import CistemProtTsCtffind
from gapstop.objects import SetOfPytomScoreTomograms
from imod.constants import OUTPUT_TILTSERIES_NAME
from imod.protocols import ProtImodExcludeViews
from pwem.objects import VolumeMask
from tomo.objects import SetOfCTFTomoSeries, SetOfTiltSeries, SetOfCoordinates3D, TiltSeries, CTFTomoSeries, \
    SetOfTomograms
from pytom_tm.protocols import ProtPytomTemplateMatching, ProtPytomExtractCoordinates
from pwem.protocols import ProtImportMask, ProtImportVolumes
from pyworkflow.tests import setupTestProject, DataSet
from pyworkflow.utils import magentaStr, cyanStr
from tomo.protocols import ProtImportTs, ProtImportTsCTF, ProtImportTomograms
from tomo.protocols.protocol_import_ctf import ImportChoice
from tomo.protocols.protocol_import_tomograms import OUTPUT_NAME
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer

TS_03 = 'TS_03'
TS_54 = 'TS_54'


class TestPytomTM(TestBaseCentralizedLayer):
    unbinnedSRate = DataSetRe4STATuto.unbinnedPixSize.value
    binFactor4 = 4
    binFactor8 = 8
    nTomos = 2
    expectedTomoDims = [480, 464, 140]
    excludedViewsDict = {
        TS_03: [0, 1, 2, 38, 39],
        TS_54: [0, 1, 38, 39, 40]
    }
    ctfExcludedViewsDict = {
        TS_03: [0, 1, 38, 39],
        TS_54: [0, 39, 40]
    }
    UNMODIFIED = 'unmodified'
    EXC_VIEWS = 'exc. views'
    RE_STACKED = 're-stacked'
    particleDiameter = 10

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls.sRateBin8 = cls.unbinnedSRate * cls.binFactor8

    @classmethod
    def _runPreviousProtocols(cls, eVCtf: bool = False, eVTs: bool = False) \
            -> Tuple[SetOfTiltSeries, SetOfCTFTomoSeries]:
        print(cyanStr('--------------------------------- RUNNING PREVIOUS PROTOCOLS ---------------------------------'))
        importedTs = cls._runImportTs()
        if eVTs:
            cls._excludeSetViews(importedTs)
        importedCtfs = cls._runImportCtf(importedTs)
        if eVCtf:
            cls._excludeSetViews(importedCtfs, excludedViewsDict=cls.ctfExcludedViewsDict)
        cls.tomoNoFidBin8 = cls._runImportTomograms()
        cls.maskBin8 = cls._runImportMaskBin8()
        cls.refBin8 = cls._runImportReferenceBin8()
        print(
            cyanStr('\n-------------------------------- PREVIOUS PROTOCOLS FINISHED ---------------------------------'))
        return importedTs, importedCtfs

    @classmethod
    def _excludeSetViews(cls,
                         inSet: Union[SetOfTiltSeries, SetOfCTFTomoSeries],
                         excludedViewsDict: Union[dict, None] = None) -> None:
        if not excludedViewsDict:
            excludedViewsDict = cls.excludedViewsDict
        objList = [obj.clone(ignoreAttrs=[]) for obj in inSet]
        for obj in objList:
            cls._excIntermediateSetViews(inSet, obj, excludedViewsDict[obj.getTsId()])

    @staticmethod
    def _excIntermediateSetViews(inSet: Union[SetOfTiltSeries, SetOfCTFTomoSeries],
                                 obj: Union[TiltSeries, CTFTomoSeries],
                                 excludedViewsList: List[int]) -> None:
        tiList = [ti.clone() for ti in obj]
        for i, ti in enumerate(tiList):
            if i in excludedViewsList:
                ti._objEnabled = False
                obj.update(ti)
        obj.write()
        inSet.update(obj)
        inSet.write()
        inSet.close()

    @classmethod
    def _runImportTs(cls) -> SetOfTiltSeries:
        print(magentaStr("\n==> Importing the tilt series:"))
        protImportTs = cls.newProtocol(ProtImportTs,
                                       filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                       filesPattern=DataSetRe4STATuto.tsPattern.value,
                                       exclusionWords=DataSetRe4STATuto.exclusionWordsTs03ts54.value,
                                       anglesFrom=2,  # From tlt file
                                       voltage=DataSetRe4STATuto.voltage.value,
                                       magnification=DataSetRe4STATuto.magnification.value,
                                       sphericalAberration=DataSetRe4STATuto.sphericalAb.value,
                                       amplitudeContrast=DataSetRe4STATuto.amplitudeContrast.value,
                                       samplingRate=cls.unbinnedSRate,
                                       doseInitial=DataSetRe4STATuto.initialDose.value,
                                       dosePerFrame=DataSetRe4STATuto.dosePerTiltImg.value,
                                       tiltAxisAngle=DataSetRe4STATuto.tiltAxisAngle.value)

        cls.launchProtocol(protImportTs)
        tsImported = getattr(protImportTs, 'outputTiltSeries', None)
        return tsImported

    @classmethod
    def _runImportCtf(cls, importedTs: SetOfTiltSeries) -> SetOfCTFTomoSeries:
        print(magentaStr("\n==> Importing the CTFs:"))
        protImportCtf = cls.newProtocol(ProtImportTsCTF,
                                        filesPath=cls.ds.getFile(DataSetRe4STATuto.tsPath.value),
                                        filesPattern=DataSetRe4STATuto.ctfPattern.value,
                                        importFrom=ImportChoice.CTFFIND.value,
                                        inputSetOfTiltSeries=importedTs)
        cls.launchProtocol(protImportCtf)
        outputCtfs = getattr(protImportCtf, protImportCtf._possibleOutputs.CTFs.name, None)
        return outputCtfs

    @classmethod
    def _runImportTomograms(cls) -> SetOfTomograms:
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomos = cls.newProtocol(ProtImportTomograms,
                                          filesPath=cls.ds.getFile(DataSetRe4STATuto.tomogramsNoFidPath.value),
                                          filesPattern='*.mrc',
                                          samplingRate=cls.sRateBin8)  # Bin 8
        cls.launchProtocol(protImportTomos)
        outTomos = getattr(protImportTomos, OUTPUT_NAME, None)
        return outTomos

    @classmethod
    def _runImportMaskBin8(cls) -> VolumeMask:
        print(magentaStr("\n==> Resizing the reference to bin 8:"))
        protImportMask3D = cls.newProtocol(ProtImportMask,
                                           maskPath=cls.ds.getFile(DataSetRe4STATuto.maskHivBin8.value),
                                           samplingRate=cls.sRateBin8)  # Bin 8
        cls.launchProtocol(protImportMask3D)
        mask = getattr(protImportMask3D, 'outputMask', None)
        return mask

    @classmethod
    def _runImportReferenceBin8(cls):
        print(magentaStr("\n==> Importing the reference volume:"))
        protImportRef = cls.newProtocol(ProtImportVolumes,
                                        filesPath=cls.ds.getFile(DataSetRe4STATuto.referenceHivBin8.value),
                                        samplingRate=cls.sRateBin8)  # Bin 8
        cls.launchProtocol(protImportRef)
        return getattr(protImportRef, ProtImportVolumes._possibleOutputs.outputVolume.name, None)

    @classmethod
    def _runExcludeViewsProt(cls,
                             inTsSet: SetOfTiltSeries,
                             objLabel: str = None) -> SetOfTiltSeries:
        print(magentaStr("\n==> Running the TS exclusion of views:"))
        protExcViews = cls.newProtocol(ProtImodExcludeViews, inputSetOfTiltSeries=inTsSet)
        if objLabel:
            protExcViews.setObjLabel(objLabel)
        cls.launchProtocol(protExcViews)
        outTsSet = getattr(protExcViews, OUTPUT_TILTSERIES_NAME, None)
        return outTsSet

    @classmethod
    def _runCistemEstimateCtf(cls, inTsSet: SetOfTiltSeries) -> SetOfCTFTomoSeries:
        print(magentaStr("\n==> Estimating the CTF with Cistem:"))
        protEstimateCtf = cls.newProtocol(CistemProtTsCtffind,
                                          inputTiltSeries=inTsSet,
                                          lowRes=50,
                                          highRes=5,
                                          minDefocus=5000,
                                          maxDefocus=50000)
        cls.launchProtocol(protEstimateCtf)
        return getattr(protEstimateCtf, CistemProtTsCtffind._possibleOutputs.CTFs.name, None)

    @classmethod
    def _genReStackedCtf(cls) -> SetOfCTFTomoSeries:
        tsSet = cls._runImportTs()
        # Exclude some views from the TS at metadata level
        cls._excludeSetViews(tsSet, excludedViewsDict=cls.ctfExcludedViewsDict)
        # Re-stack that TS
        reStackedTsSet = cls._runExcludeViewsProt(tsSet)
        # Estimate the CTF using the re-stacked TS
        return cls._runCistemEstimateCtf(reStackedTsSet)

    def _runPytomTM(self,
                    inCtfSet: SetOfCTFTomoSeries,
                    inTsSet: SetOfTiltSeries,
                    ctfSetMsg: str,
                    tsSetMsg: str) \
            -> Union[SetOfPytomScoreTomograms, None]:
        print(magentaStr(f"\n==> Running the Pytom_TM:"
                         f"\n\t- CTFs: {ctfSetMsg}"
                         f"\n\t- Tilt-series = {tsSetMsg}"))
        protPytomTM = self.newProtocol(ProtPytomTemplateMatching,
                                       inTomos=self.tomoNoFidBin8,
                                       inCtfSet=inCtfSet,
                                       inTsSet=inTsSet,
                                       reference=self.refBin8,
                                       mask=self.maskBin8,
                                       doInvertRefContrast=True,
                                       nTiles=8,
                                       coneSampling=15,
                                       rotSymDeg=6)
        objLabel = f'ts {tsSetMsg}, ctf {ctfSetMsg}'
        protPytomTM.setObjLabel(objLabel)
        self.launchProtocol(protPytomTM)
        return getattr(protPytomTM, protPytomTM._possibleOutputs.scoreTomograms.name, None)

    def _checkScoredTomos(self, scoreTomos: SetOfPytomScoreTomograms) -> None:
        # Check the results of the pytom_TM
        self.checkTomograms(inTomoSet=scoreTomos,
                            expectedSetSize=self.nTomos,
                            expectedSRate=self.sRateBin8,
                            expectedDimensions=self.expectedTomoDims)
        # PytomScoreTomogram specific attributes
        for tomo in scoreTomos:
            self.assertTrue(exists(tomo.getTomoFile()))
            self.assertTrue(exists(tomo.getAnglesMap()))
            self.assertTrue(exists(tomo.getAngleList()))
            self.assertGreater(tomo.getTomoNum(), 0)
            self.assertEqual(tomo.getSymmetry(), 'C6')

    def _runPytomExtractCoords(self, scoreTomos: SetOfPytomScoreTomograms) \
            -> Union[SetOfCoordinates3D, None]:
        print(magentaStr("\n==> Pytom_TM == > extracting the coordinates"))
        protPytomExtract = self.newProtocol(ProtPytomExtractCoordinates,
                                              inScoreTomos=scoreTomos,
                                              scoresThreshold=0.09,
                                              percentile=99.9,
                                              partDiameter=self.particleDiameter,
                                              numberOfCoords=-1)
        self.launchProtocol(protPytomExtract)
        return getattr(protPytomExtract, protPytomExtract._possibleOutputs.coordinates.name, None)

    def _checkExtractedCoords(self, coords: SetOfCoordinates3D) -> None:
        self.checkCoordinates(outCoords=coords,
                              expectedBoxSize=self.particleDiameter,
                              expectedSRate=self.sRateBin8,
                              orientedParticles=True,
                              orientedTolPercent=0.01)  # 1%
        self.assertTrue(coords.getSize() > 2000)

    def _runTestPytomTM(self,
                        inCtfSet: SetOfCTFTomoSeries,
                        inTsSet: SetOfTiltSeries,
                        ctfSetMsg: str,
                        tsSetMsg: str) -> None:
        # Run the template matching
        scoreTomos = self._runPytomTM(inCtfSet=inCtfSet,
                                      inTsSet=inTsSet,
                                      ctfSetMsg=ctfSetMsg,
                                      tsSetMsg=tsSetMsg, )
        # Check the scored tomograms
        self._checkScoredTomos(scoreTomos)
        # Run the coordinate extraction
        coords = self._runPytomExtractCoords(scoreTomos)
        # Check the results of gapStop's extraction
        self._checkExtractedCoords(coords)

    def testPytom(self):
        importedTs, importedCtfs = self._runPreviousProtocols()
        self._runTestPytomTM(inCtfSet=importedCtfs,
                             inTsSet=importedTs,
                             ctfSetMsg=self.UNMODIFIED,
                             tsSetMsg=self.UNMODIFIED)

    def testPytom_EV_Ctf(self):
        importedTs, importedCtfs = self._runPreviousProtocols(eVCtf=True)
        self._runTestPytomTM(inCtfSet=importedCtfs,
                             inTsSet=importedTs,
                             ctfSetMsg=self.EXC_VIEWS,
                             tsSetMsg=self.UNMODIFIED)

    def testPytom_EV_Ts(self):
        importedTs, importedCtfs = self._runPreviousProtocols(eVTs=True)
        self._runTestPytomTM(inCtfSet=importedCtfs,
                             inTsSet=importedTs,
                             ctfSetMsg=self.UNMODIFIED,
                             tsSetMsg=self.EXC_VIEWS)

    def testPytom_EV_Restacked_Ctf(self):
        importedTs, _ = self._runPreviousProtocols()
        ctfSetReStacked = self._genReStackedCtf()  # Gen a CTF estimated on a re-stacked TS
        self._runTestPytomTM(inCtfSet=ctfSetReStacked,
                             inTsSet=importedTs,
                             ctfSetMsg=self.EXC_VIEWS,
                             tsSetMsg=self.UNMODIFIED)

    def testPytom_EV_Restacked_Ts(self):
        importedTs, importedCtfs = self._runPreviousProtocols(eVTs=True)
        tsSetReStacked = self._runExcludeViewsProt(importedTs)  # Re-stack the TS
        self._runTestPytomTM(inCtfSet=importedCtfs,
                             inTsSet=tsSetReStacked,
                             ctfSetMsg=self.UNMODIFIED,
                             tsSetMsg=self.EXC_VIEWS)