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
import hashlib
import logging
import os
import traceback
import typing
from enum import Enum
from os.path import join
from typing import List, Optional, Union

import numpy as np
from fidder.protocols.protocol_detect_and_erase_fiducials import MASK_SUFFIX

from pwem.convert.headers import setMRCSamplingRate
from pwem.emlib.image.image_readers import MRCImageReader
from pwem.objects import VolumeMask, Volume
from pwem.protocols import EMProtocol
from pytom_tm import Plugin
from pytom_tm.constants import IN_TOMOS, REF_VOL, IN_MASK, TOMO_MASKS, IN_TS_SET, IN_CTF_SET, MRC_EXT, DEFOCUS_EXT, \
    TILT_ANGLES_EXT, DOSE_EXT, DOSE_SUFFIX, TOMO_SUFFIX, DEFOCUS_HAND_OFF, DEFOCUS_HAND_NEG, DEFOCUS_HAND_POS, \
    BASE_SEED, SCORE_SUFFIX, JSON_EXT, JSON_SUFFIX
from pytom_tm.objects import SetOfPytomScoreTomograms, PytomScoreTomogram
from pyworkflow import BETA
from pyworkflow.object import Pointer, String, Set
from pyworkflow.protocol import PointerParam, BooleanParam, FloatParam, IntParam, StringParam, LEVEL_ADVANCED, \
    EnumParam, GT, GPU_LIST
from pyworkflow.utils import Message, cyanStr, makePath, redStr
from tomo.objects import SetOfTiltSeries, SetOfTomograms, SetOfCTFTomoSeries, SetOfTomoMasks, TiltSeries, TiltImage
from tomo.utils import getObjFromRelation, getCommonTsAndCtfElements, \
    getTsIdsIntersection, getTsIdsDicts, invertContrast, convertOrLink, genDefocusFileFromScipion

logger = logging.getLogger(__name__)


class Pytom_tm_output(Enum):
    scoreTomograms = SetOfPytomScoreTomograms  # output name is scoreTomograms within a specific class SetOfPytomTOmograms


class ProtPytomTemplateMatching(EMProtocol):
    """GPU-accelerated template matching for cryo-electron tomography,
    originally developed in PyTom, as a standalone Python package that is run from the command line.
    """

    _label = 'template matching'
    _devStatus = BETA
    _possibleOutputs = Pytom_tm_output
    _program = 'pytom_match_template.py'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tomoDict = {}
        self.tsDict = {}
        self.ctfDict = {}
        self.tomoMaskDict = {}
        self.refName = ''
        self.failedTsIds = []
        self.failedTsIdsStr = String()
        self.samplingRate = -1

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Tomograms')

        form.addParam(IN_CTF_SET, PointerParam,
                      pointerClass='SetOfCTFTomoSeries',
                      label="CTF tomo series",
                      important=True)

        form.addParam(IN_TS_SET, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      important=True,
                      label='Tilt-series',
                      help='Used to get the tilt angles.')

        form.addParam(REF_VOL, PointerParam,
                      pointerClass='Volume',
                      important=True,
                      label="Reference volume")

        form.addParam('invert_contrast', BooleanParam,
                      default=True,
                      label='Invert reference contrast?',
                      important=True,
                      help='The contrast of the template has to be the same as of the tomogram. If the '
                           'tomogram has features in black (which is typically for cryoET) then the template '
                           'has to have the same representation. For example, Relion outputs inverted '
                           'contrast (features are white) and such maps have to be inverted prior running '
                           'the Pytom_tm.')

        group = form.addGroup('Mask')
        group.addParam(IN_MASK, PointerParam,
                       pointerClass=VolumeMask,
                       important=True,
                       label='Reference mask')

        group.addParam('non_spherical_mask', BooleanParam,
                       label='Is it a spherical mask?',
                       default=False,
                       help="Flag to set when the mask is not spherical. It adds the required "
                            "computations for non-spherical masks and roughly doubles computation time."
                       )

        form.addSection(label='Angular Search')
        form.addParam('angular_search', FloatParam,
                      label='Angular Search (deg)',
                      important=True,
                      default=7.,
                      allowsNull=True,
                      help="Angular increment of template search. "
                           "If empty it will be computed from the maximum resolution (2*voxel_size or low_pass)"
                           "and the particle box size.")

        form.addParam('z_axis_rotational_symmetry', IntParam,
                      label='Z Axis Rotational Symmetry',
                      default=1,
                      allowsNull=False,
                      validators=[GT(0)],
                      help="Integer value indicating the rotational symmetry of the template around "
                           "the z-axis. The length of the rotation search will be shortened through "
                           "division by this value. Only works for template symmetry around the z-axis.")

        form.addSection(label='Volume control')
        form.addParam('volume_split', StringParam,
                      label='Volume Split',
                      default='1 1 1',
                      help="Split the volume into smaller parts for the search, "
                           "can be relevant if the volume does not fit into GPU memory. "
                           "Format is x y z, e.g. --volume-split 1 2 1")

        form.addParam('defocus_handedness', EnumParam,
                      choices=[DEFOCUS_HAND_NEG, DEFOCUS_HAND_OFF, DEFOCUS_HAND_POS],
                      display=EnumParam.DISPLAY_HLIST,
                      label='Defocus Handedness',
                      expertLevel=LEVEL_ADVANCED,
                      default=DEFOCUS_HAND_OFF,
                      condition='volume_split != "1 1 1"',
                      help="Specify the defocus handedness for defocus gradient correction of the "
                           "CTF in each subvolumes. The more subvolumes in x and z, "
                           "the finer the defocus gradient will be corrected, at the cost of "
                           "increased computing time. It will only have effect for very clean and "
                           "high-resolution data, such as isolated macromolecules. "
                           "A value of 0 means no defocus gradient correction (default), 1 means "
                           "correction assuming correct handedness (as specified in Pyle and "
                           "Zianetti (2021)), -1 means the handedness will be inverted. If uncertain "
                           "better to leave off as an inverted correction might hamper results."
                      )

        search_indices_group = form.addGroup('Search Indices')
        x_index = search_indices_group.addLine('X axis',
                                               help="Start and end indices of the search along the x-axis")

        x_index.addParam('xmin',
                         IntParam,
                         label='X min',
                         allowsNull=True)

        x_index.addParam('xmax',
                         IntParam,
                         label='X max',
                         allowsNull=True)

        y_index = search_indices_group.addLine('Y axis',
                                               help="Start and end indices of the search along the y-axis")

        y_index.addParam('ymin',
                         IntParam,
                         label='Y min',
                         allowsNull=True)

        y_index.addParam('ymax',
                         IntParam,
                         label='Y max',
                         allowsNull=True)

        z_index = search_indices_group.addLine('Z axis',
                                               help="Start and end indices of the search along the z-axis")

        z_index.addParam('zmin',
                         IntParam,
                         label='Z min',
                         allowsNull=True)

        z_index.addParam('zmax',
                         IntParam,
                         label='Z max',
                         allowsNull=True)

        form.addParam(TOMO_MASKS, PointerParam,
                      pointerClass=SetOfTomoMasks,
                      label='Tomogram masks (segmentations)',
                      allowsNull=True,
                      help="Here you can provide a set of masks for matching with dimensions (in pixels) "
                           "equal to the tomogram. If a subvolume only has values <= 0 for this mask it "
                           "will be skipped."
                      )

        form.addSection(label='Filter Control')
        form.addParam('low_pass', FloatParam,
                      label='Low Pass filter ',
                      allowsNull=True,
                      help="Apply a low-pass filter to the tomogram and template. Generally desired "
                           "if the template was already filtered to a certain resolution. "
                           "Value is the resolution in A."
                      )
        form.addParam('high_pass', FloatParam,
                      label='High Pass filter ',
                      allowsNull=True,
                      help="Apply a high-pass filter to the tomogram and template to reduce "
                           "correlation with large low frequency variations. Value is a resolution in A, "
                           "e.g. 500 could be appropriate as the CTF is often incorrectly modelled "
                           "up to 50nm."
                      )
        form.addParam('spectral_whitening', BooleanParam,
                      label='Spectral Whitening',
                      default=False,
                      help="Calculate a whitening filtering from the power spectrum of the tomogram; "
                           "apply it to the tomogram patch and template. Effectively puts more weight on "
                           "high resolution features and sharpens the correlation peaks.")
        form.addParam('per_tilt_weighting', BooleanParam,
                      label='per-tilt-weighting',
                      default=True,
                      expertLevel=LEVEL_ADVANCED,
                      help="Flag to activate per-tilt-weighting. The base functionality creates a fanned wedge where each tilt is "
                           "weighted by cos(tilt_angle)."
                      )

        form.addSection(label='Additional Parameters')
        form.addParam('random_phase_correction', BooleanParam,
                      label='Random Phase Correction',
                      default=False,
                      help="Run template matching simultaneously with a phase randomized version of "
                           "the template, and subtract this 'noise' map from the final score map. "
                           "For this method please see STOPGAP as a reference: "
                           "https://doi.org/10.1107/S205979832400295X ."
                      )
        form.addParam('rng_seed', IntParam,
                      label='Phase randomization range seed',
                      allowsNull=True,
                      condition='random_phase_correction',
                      help="Specify a seed for the random number generator used for phase "
                           "randomization for consistent results! If empty a random number will be provided.",
                      )
        form.addParam('extraParams', StringParam,
                      label='Additional parameters',
                      help="In this box command-line arguments may be provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options and/or expert use of the program, e.g: \n"
                           "--half-precision \n"
                           "--warp-xml-file \n")

        form.addHidden(GPU_LIST, StringParam,
                       default='0',
                       label="Choose GPU IDs",
                       help="")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # FRANCESCA
        # import os
        # fname = "/home/francesa/test_FC.txt"
        # if os.path.exists(fname):
        #     os.remove(fname)
        # with open(fname, "a+") as fjj:
        #     fjj.write(f'FRANCESCA--------->onDebugMode PID {os.getpid()}')
        #     print(f'FRANCESCA--------->onDebugMode PID {os.getpid()}')
        # import time
        # time.sleep(10)
        # # FRANCESCA_END

        self._initialize()
        closeSetStepDeps = []

        cRId = self._insertFunctionStep(self.convertReferenceStep,
                                        prerequisites=[],
                                        needsGPU=False)

        for tsId in self.tomoDict.keys():
            cInputId = self._insertFunctionStep(self.convertInputStep, tsId,
                                                prerequisites=cRId,
                                                needsGPU=False)
            tmId = self._insertFunctionStep(self.templateMatchingStep, tsId,
                                            prerequisites=cInputId,
                                            needsGPU=True)
            cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                              prerequisites=tmId,
                                              needsGPU=False)
            closeSetStepDeps.append(cOutId)
        self._insertFunctionStep(self.closeOutputSetStep,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

        # -------------------------- STEPS functions ------------------------------

    def _initialize(self):
        tsSet = self._getFormAttrib(IN_TS_SET)
        tomoSet = self._getFormAttrib(IN_TOMOS)
        ctfSet = self._getFormAttrib(IN_CTF_SET)
        tomoMasks = self._getFormAttrib(TOMO_MASKS)
        self.samplingRate = tomoSet.getSamplingRate()
        # self.refName = self._genConvertedOrLinkedRefName(REF_VOL)
        # self.maskName = self._genConvertedOrLinkedRefName(IN_MASK)
        # self.tomosSRate = tomoSet.getSamplingRate()
        # self.tomosBinning = self._getTomogramsBinning()

        # Compute matching TS id among coordinates, the tilt-series and the CTFs, they all could be a subset

        if tomoMasks:
            presentTsIds = getTsIdsIntersection(tsSet, tomoSet, ctfSet, tomoMasks)
            self.tsDict, self.tomoDict, self.ctfDict, self.tomoMaskDict = (
                getTsIdsDicts(tsSet, tomoSet, ctfSet, tomoMasks, present_ts_ids=presentTsIds))

            # self.tomoMaskDict = {tsId: tomoMasks.clone() for tomoMasks in tomoMasks.iterItems() if (tsId:=tomoMasks.getTsId()) in presentTsIds}

        else:
            presentTsIds = getTsIdsIntersection(tsSet, tomoSet, ctfSet)

            self.tsDict, self.tomoDict, self.ctfDict = (
                getTsIdsDicts(tsSet, tomoSet, ctfSet, present_ts_ids=presentTsIds))

    def convertReferenceStep(self):
        logger.info(cyanStr(f"Converting the reference in the required format...'"))
        try:
            # Convert or link the reference
            ref = self._getFormAttrib(REF_VOL)
            inRefFile = ref.getFileName()
            refFile = self.getReferenceFileName()
            samplingRate = ref.getSamplingRate()

            if self.invert_contrast.get():
                invertContrast(inRefFile, refFile, samplingRate)
            else:  # if the contrast inversion is not needed, only convert to .mrc
                convertOrLink(inRefFile, refFile, samplingRate)

            # Convert or link the mask
            mask = self._getFormAttrib(IN_MASK)
            inMaskFile = mask.getFileName()
            outMaskFile = self.getMaskFileName()
            samplingRate = mask.getSamplingRate()
            convertOrLink(inMaskFile, outMaskFile, samplingRate)

        except Exception as e:
            raise Exception(f'Reference conversion failed with the exception -> {e}')

    def convertInputStep(self, tsId: str):

        try:
            tsDir = self._getCurrentTomoDir(tsId)
            tsTmpDir = self._getCurrentTomoTmpDir(tsId)
            makePath(tsDir, tsTmpDir)

            tomo = self.tomoDict[tsId]
            inTomoFile = tomo.getFileName()
            outTomoFile = self._getConvertedOrLinkedName(tsId, suffix=TOMO_SUFFIX)
            convertOrLink(inTomoFile, outTomoFile, samplingRate=tomo.getSamplingRate())

            if self.tomoMaskDict:
                tomomask = self.tomoMaskDict[tsId]
                inTomoMaskFile = tomomask.getFileName()
                outTomoMaskFile = self._getConvertedOrLinkedName(tsId, suffix=MASK_SUFFIX)
                convertOrLink(inTomoMaskFile, outTomoMaskFile, samplingRate=tomo.getSamplingRate())

            ts = self.tsDict[tsId]
            ctf = self.ctfDict[tsId]

            presentAcqOrders = getCommonTsAndCtfElements(ts, ctf)
            if len(presentAcqOrders) == 0:
                raise Exception(f'tsId = {tsId} -> No common acquisition orders found between the '
                                f'tilt-series and the CTF.')

            logger.info(cyanStr(f"tsId = {tsId} -> present acquisition orders in both "
                                f"the tilt-series and the CTF are {presentAcqOrders}.'"))

            outDefocus = self._getInputFileName(tsId, DEFOCUS_EXT)
            genDefocusFileFromScipion(ctf, ts, outDefocus)

            outTilt_Angles = self._getInputFileName(tsId, TILT_ANGLES_EXT)
            ts.generateTltFile(outTilt_Angles, presentAcqOrders)

            outDosePath = self._getInputFileName(tsId, DOSE_EXT, suffix=DOSE_SUFFIX)
            self.generateDoseFile(ts, outDosePath, presentAcqOrders)

        except Exception as e:
            self.failedTsIds.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def templateMatchingStep(self, tsId: str):
        if tsId in self.failedTsIds:
            return
        try:
            logger.info(cyanStr(f'===> tsId = {tsId}: performing the template matching...'))
            Plugin.runPytom(self, self._program, self._generateArguments(tsId))

        except Exception as e:
            self.failedTsIds.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> pytom execution failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def createOutputStep(self, tsId: str):
        if tsId in self.failedTsIds:
            return
        try:

            tomo = self.tomoDict[tsId]
            scoresMap = self._getOutputFileName(tsId, suffix=SCORE_SUFFIX)
            jsonFile = self._getOutputFileName(tsId, ext=JSON_EXT, suffix=JSON_SUFFIX)
            setMRCSamplingRate(scoresMap, tomo.getSamplingRate())  # Update the apix value in file header
            scoreTomoSet = self.createOutputSet()
            # Create the corresponding scoreTomo
            scoreTomo = PytomScoreTomogram()
            scoreTomo.copyInfo(tomo)
            scoreTomo.setFileName(scoresMap)
            scoreTomo.setTomoFile(tomo.getFileName())
            scoreTomo.setJsonFile(jsonFile)

            # Append to the set and store
            scoreTomoSet.append(scoreTomo)
            scoreTomoSet.write()
            self._store(scoreTomoSet)

        except Exception as e:
            logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with exception {e}. Skipping... '))
            logger.error(traceback.format_exc())

    def closeOutputSetStep(self):
        scoreTomoSet = getattr(self, self._possibleOutputs.scoreTomograms.name, None)
        if scoreTomoSet:
            self._closeOutputSet()
        else:
            raise Exception('No Pytom scored tomograms were generated. Maybe the tomograms are too large '
                            'for the GPU/s used. Consider to bin them before and/or make tiles from the '
                            'tomogram using the parameter Volume Split.')
        if self.failedTsIds:
            self.failedTsIdsStr.set(str(self.failedTsIds))
            self._store(self.failedTsIdsStr)

    # --------------------------- INFO functions ------------------------------

    def _validate(self) -> List[str]:
        valMsg = []
        lpf = self.low_pass.get()
        hpf = self.high_pass.get()
        x_min = self.xmin.get()
        x_max = self.xmax.get()
        y_min = self.ymin.get()
        y_max = self.ymax.get()
        z_min = self.zmin.get()
        z_max = self.zmax.get()

        if not self.validate_volume_split(self.volume_split.get()):
            valMsg.append('Volume split must be a list of three integers (e.g. 1 1 1)')

        if lpf is not None and hpf is not None:
            if lpf <= 0 or hpf <= 0:
                valMsg.append('Low and high pass filter must be >= 0')
            elif lpf <= hpf:
                valMsg.append('Low pass filter value must be greater than high pass filter value')

        iter = [(x_min, x_max), (y_min, y_max), (z_min, z_max)]
        for pair in iter:
            self.check_search_values(pair[0], pair[1], valMsg)

        return valMsg

    def _summary(self) -> List[str]:
        msg = []
        if self.isFinished():
            msg.append('*Pytom_TM is composed of 2 steps*. To extract the coordinates from the scored '
                       'tomograms calculated, call the protocol *pytom - extract coordinates*.')
            failedStrs = self.failedTsIdsStr.get()
            if failedStrs:
                msg.append(f'The following tsIds were not possible to be processed: *{failedStrs}*')
        return msg

        # --------------------------- UTILS functions ------------------------------

    def check_search_values(self, val1: Optional[int], val2: Optional[int], errorList: List[str]) -> None:
        if val1 is None and val2 is None:
            return
        if val1 is None or val2 is None:
            errorList.append('If min value is filled, max value must also be filled and viceversa.')
            return
        if val1 > val2:
            errorList.append('min value must be less than max value.')

    @staticmethod
    def validate_volume_split(text: str) -> bool:
        parts = text.split()

        if len(parts) != 3:  # check if it's a three elements list
            return False
        try:
            [int(x) for x in parts]
            return True
        except ValueError:
            return False

    def _getFormAttrib(self, attribName: str, returnPointer: bool = False) -> Optional[Union[SetOfTiltSeries,
    SetOfTomograms, SetOfCTFTomoSeries, Volume, Pointer]]:
        inTsPointer = getattr(self, attribName, None)
        if not inTsPointer:
            return None
        else:
            return inTsPointer if returnPointer else inTsPointer.get()

    # if IN_TS_SET is not an input, get it yourself via:
    def _getTsFromRelations(self) -> Optional[SetOfTiltSeries]:
        inCTFs = self._getFormAttrib(IN_CTF_SET)
        return getObjFromRelation(inCTFs, self, SetOfTiltSeries)

    def _getCurrentTomoDir(self, tsId: str) -> str:
        return self._getExtraPath(tsId)

    def getReferenceFileName(self) -> str:
        return self._getTmpPath(f'Reference{MRC_EXT}')

    def getMaskFileName(self) -> str:
        return self._getTmpPath(f'Mask{MRC_EXT}')

    def _getCurrentTomoTmpDir(self, tsId: str) -> str:
        return self._getTmpPath(tsId)

    def _getConvertedOrLinkedName(self, tsId: str, suffix: str = '') -> str:
        return join(self._getCurrentTomoTmpDir(tsId), f'{tsId}{suffix}{MRC_EXT}')

    def _getInputFileName(self, tsId: str, ext: str, suffix: str = '') -> str:
        return join(self._getCurrentTomoDir(tsId), f'{tsId}{suffix}{ext}')

    def _getOutputFileName(self, tsId: str, ext: str = MRC_EXT, suffix: str = '') -> str:
        return self._getInputFileName(tsId, ext=ext, suffix=suffix)

    def generateDoseFile(self,
                         ts: TiltSeries,
                         dosePath: str,
                         presentAcqOrders: typing.Set[int]) -> None:

        doseList = []
        for ti in ts.iterItems(orderBy=TiltImage.TILT_ANGLE_FIELD):
            if ti.getAcquisitionOrder() in presentAcqOrders:
                doseList.append(ti.getAcquisition().getAccumDose())

        with open(dosePath, 'w') as f:
            f.writelines(f"{dose:0.2f}\n" for dose in doseList)
            # For parallel processing, ensure that the file is completely written and persists on disk
            f.flush()  # Empty python buffer
            os.fsync(f.fileno())  # Empty system buffer

    def getAngularStep(self) -> float:
        angular_search = self.angular_search.get()
        if angular_search:
            return angular_search
        else:
            low_pass = self.low_pass.get()
            x, y, z, _ = MRCImageReader.getDimensions(self.getReferenceFileName())
            particle_diameter = max(x, y, z)
            max_res = max(
                2 * self.samplingRate, low_pass if low_pass is not None else 0
            )
            return np.rad2deg(max_res / particle_diameter)

    def _generateArguments(self, tsId: str) -> str:
        x_min = self.xmin.get()
        x_max = self.xmax.get()
        y_min = self.ymin.get()
        y_max = self.ymax.get()
        z_min = self.zmin.get()
        z_max = self.zmax.get()

        lpf = self.low_pass.get()
        hpf = self.high_pass.get()

        ts = self.tsDict[tsId]
        acquisition = ts.getAcquisition()

        # ctf = self.ctfDict[tsId]
        # phaseShift = ctf.getPhaseShift()

        tomo = self.tomoDict[tsId]
        ctfCorrected = tomo.ctfCorrected()

        gpu = ' '.join([str(el) for el in self.getGpuList()])

        cmd = [
            f'--template {self.getReferenceFileName()}',
            f'--tomogram {self._getConvertedOrLinkedName(tsId, suffix=TOMO_SUFFIX)}',
            f'--destination {self._getCurrentTomoDir(tsId)}',
            f'--mask {self.getMaskFileName()}',
            f'--angular-search {self.getAngularStep()}',
            f'--z-axis-rotational-symmetry {self.z_axis_rotational_symmetry.get()}',
            f'--volume-split {self.volume_split.get()}',
            f'--tilt-angles {self._getInputFileName(tsId, TILT_ANGLES_EXT)}',
            f'--voxel-size-angstrom {self.samplingRate:.3f}',
            f'--dose-accumulation {self._getInputFileName(tsId, DOSE_EXT, suffix=DOSE_SUFFIX)}',
            f'--defocus {self._getInputFileName(tsId, DEFOCUS_EXT)}',
            f'--amplitude-contrast {acquisition.getAmplitudeContrast()}',
            f'--spherical-aberration {acquisition.getSphericalAberration()}',
            f'--voltage {acquisition.getVoltage()}',
            f'--gpu-ids {gpu}',
            '--log info'

        ]
        if x_min:
            cmd.append(f'--search-x {x_min} {x_max}')
        if y_min:
            cmd.append(f'--search-y {y_min} {y_max}')
        if z_min:
            cmd.append(f'--search-z {z_min} {z_max}')

        if self.tomoMaskDict:
            cmd.append(f'--tomogram-mask {self._getConvertedOrLinkedName(tsId, suffix=MASK_SUFFIX)}')

        if self.non_spherical_mask.get():
            cmd.append('--non-spherical-mask')

        if self.per_tilt_weighting.get():
            cmd.append('--per-tilt-weighting')
        # TODO
        # if lpf:
        #     cmd.append(f'--low-pass {lpf:.2f}')
        #
        # if hpf:
        #     cmd.append(f'--high-pass {hpf:.2f}')

        # TODO
        # if phaseShift:
        #     cmd.append(f'--phase-shift {phaseShift}')

        if ctfCorrected:
            cmd.append('--tomogram-ctf-model phase-flip')

        if self.volume_split.get() != '1 1 1':
            cmd.append(f'--defocus-handedness {self.defocus_handedness.get()}')

        if self.spectral_whitening.get():
            cmd.append(f'--spectral-withening')

        if self.random_phase_correction.get():
            cmd.append('--random-phase-correction')

            rng_seed = self.rng_seed.get()
            rng_seed = rng_seed if rng_seed else self.seed_from_name(tsId)
            cmd.append(f'--rng-seed {rng_seed}')

        return ' '.join(cmd)

    @staticmethod
    def seed_from_name(tsId: str, base_seed: int = BASE_SEED) -> int:
        """Derive a deterministic per-tomogram seed from a base seed and a name.

           Hashes ``name`` (e.g. a tomogram filename) with SHA-256 and combines it
           with ``base_seed`` to produce a reproducible integer seed for the random
           number generator used in phase randomization. This ensures each
           tomogram gets its own independent noise map, while remaining
           deterministic across runs regardless of processing order.

           Args:
               base_seed: Fixed base seed for reproducibility across runs.
               tsId: Unique identifier for the tomogram (e.g. filename), used to
                   derive a distinct seed.

           Returns:
               An integer seed in the range [0, 2**31 - 2], suitable for
               ``numpy.random.default_rng`` or similar RNG constructors.
           """
        h = int(hashlib.sha256(tsId.encode()).hexdigest(), 16)
        return (base_seed + h) % (2 ** 31 - 1)

    def createOutputSet(self) -> SetOfPytomScoreTomograms:
        scoreTomoSet = getattr(self, self._possibleOutputs.scoreTomograms.name, None)
        if scoreTomoSet:
            scoreTomoSet.enableAppend()
        else:
            inTomosPointer = self._getFormAttrib(IN_TOMOS, returnPointer=True)
            inTomos = inTomosPointer.get()
            scoreTomoSet = SetOfPytomScoreTomograms.create(self._getPath(), template="scoreTomograms%s")
            scoreTomoSet.copyInfo(inTomos)
            scoreTomoSet.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{self._possibleOutputs.scoreTomograms.name: scoreTomoSet})
            self._defineSourceRelation(inTomosPointer, scoreTomoSet)

        return scoreTomoSet
