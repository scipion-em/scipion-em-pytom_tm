import json
import logging
import shutil
import traceback
from enum import Enum
from os.path import join, basename
from typing import List, Optional

import numpy as np
from emtable import Table

from pwem.convert import transformations
from pytom_tm import Plugin
from pytom_tm.constants import IN_TM_PROTOCOL, TOMO_MASKS, MASK_PYTOM_TM, MASK_OTHER, MASK_SUFFIX, IN_TOMOS, \
    SCORE_SUFFIX, ANGLES_SUFFIX
from pytom_tm.objects import SetOfPytomScoreTomograms
from pytom_tm.protocols.protocol_base import ProtPytomBase
from pyworkflow import BETA
from pyworkflow.object import String, Set
from pyworkflow.protocol import PointerParam, IntParam, GT, FloatParam, GE, LE, StringParam, EnumParam
from pyworkflow.utils import Message, cyanStr, redStr, yellowStr
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SetOfCoordinates3D, SetOfTomoMasks, Tomogram, Coordinate3D
from tomo.utils import getTsIdsDicts, getTsIdsIntersection, check_sr_and_size, convertOrLink

logger = logging.getLogger(__name__)


class Pytom_extract_outputs(Enum):
    coordinates = SetOfCoordinates3D


class ProtPytomExtractCoordinates(ProtPytomBase):
    """Extract coordinates from Pytom score tomograms."""

    _label = 'extract coordinates from Pytom'
    _devStatus = BETA
    _possibleOutputs = Pytom_extract_outputs
    _program = 'pytom_extract_candidates.py'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.scoreTomoDict = {}
        self.tomoMasksDict = {}
        self.inTomoDict = {}
        self.failedTsIds = []
        self.failedTsIdsStr = String()
        self.sRate = -1

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TM_PROTOCOL, PointerParam,
                      pointerClass='ProtPytomTemplateMatching',
                      important=True,
                      label='Protocol Pytom TM')
        form.addParam('mask_choice', EnumParam,
                      choices=['No TomoMasks', 'Same as input protocol', 'Other TomoMasks'],
                      display=EnumParam.DISPLAY_HLIST,
                      label='TomoMasks (segmentations)',
                      default=MASK_PYTOM_TM,
                      help="Here you can provide a mask for the extraction with dimensions "
                           "(in pixels) equal to the tomogram. All values in the mask that are smaller or "
                           "equal to 0 will be removed, all values larger than 0 are considered regions "
                           "of interest. It can be used to extract annotations only within a specific "
                           "cellular region. If the job was run with a tomogram mask, this file will be "
                           "used instead of the job mask."
                      )
        form.addParam(TOMO_MASKS, PointerParam,
                      pointerClass=SetOfTomoMasks,
                      label='Tomogram masks (segmentations)',
                      allowsNull=True,
                      condition=f'mask_choice == {MASK_OTHER}',
                      help="All values in the mask that are smaller or "
                           "equal to 0 will be removed, all values larger than 0 are considered regions "
                           "of interest. It can be used to extract annotations only within a specific "
                           "cellular region."
                      )
        form.addParam('number_of_particles', IntParam,
                      label='Number of particles',
                      allowsNull=False,
                      default=1000,
                      validators=[GT(0)],
                      help="Maximum number of particles to extract from tomogram."
                      )

        form.addParam('number_false_positives', FloatParam,
                      label='Number of false positives',
                      important=True,
                      validators=[GE(0), LE(1)],
                      default=1.0,
                      help="Number of false positives to determine the false alarm rate. Here one "
                           "can increase the recall of the particle of interest at the expense "
                           "of more false positives. The default value of 1 is recommended for "
                           "particles that can be distinguished well from the background (high "
                           "specificity). The value can also be set between 0 and 1 to make "
                           "the cut-off more restrictive."  # TODO check description
                      )
        form.addParam('particle_diameter', FloatParam,
                      label='Particle diameter (angst)',
                      important=True,
                      validators=[GT(0)],
                      help="Particle diameter of the template in Angstrom. It is used during "
                           "extraction to remove areas around peaks to prevent double extraction. "
                           "Minimal peak-to-peak distance after extraction will be diameter/2."
                           "If not previously specified, this option is required. If "
                           "specified in pytom_match_template, this is optional and "
                           "can be used to overwrite it, which might be relevant for strongly "
                           "elongated particles--where the angular sampling should be "
                           "determined using its long axis but the extraction mask should use its "
                           "short axis."
                      )
        form.addParam('tophat_filter_con', IntParam,
                      label='Tophat filter connectivity',
                      validators=[GE(0)],
                      default=0,
                      help="Set kernel connectivity for ndimage binary structure used for the "
                           "tophat transform. Integer value in range 1-3. 1 is the most "
                           "restrictive, 3 the least restrictive. Generally recommended to "
                           "leave at 1. 0 value means that tophat filter (attempt to filter "
                           "only sharp correlation peaks with a tophat transform) is disabled.",
                      )
        form.addSection(label='Additional Parameters')
        form.addParam('extraParams', StringParam,
                      label='Additional parameters',
                      help="In this box command-line arguments may be provided that are not generated by the GUI. This "
                           "may be useful for testing developmental options and/or expert use of the program, e.g: \n"
                           "--cut-off \n"
                      )

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        pId = []

        for tsId in self.scoreTomoDict.keys():

            pId = self._insertFunctionStep(self.convertInputStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=False)
            pId = self._insertFunctionStep(self.extractCoordinatesStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=True)
            pId = self._insertFunctionStep(self.createOutputStep, tsId,
                                           prerequisites=pId,
                                           needsGPU=False)
            closeSetStepDeps.append(pId)
        self._insertFunctionStep(self.closeOutputSetStep,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

        # -------------------------- STEPS functions ------------------------------

    def _initialize(self):
        mask_choice = self.mask_choice.get()
        scoreTomos = self.getScoreTomos()
        tomoMasks = self.getTomoMasks()
        self.sRate = scoreTomos.getSamplingRate()

        input_protocol = self._getFormAttrib(IN_TM_PROTOCOL)
        inTomoSet = getattr(input_protocol, IN_TOMOS).get()
        presentTsIds = getTsIdsIntersection(scoreTomos, inTomoSet)

        self.scoreTomoDict, self.inTomoDict = getTsIdsDicts(scoreTomos, inTomoSet, present_ts_ids=presentTsIds)

        if tomoMasks:
            presentTsIds = None
            if mask_choice == MASK_OTHER:
                presentTsIds = getTsIdsIntersection(scoreTomos, tomoMasks, inTomoSet)
            self.scoreTomoDict, self.tomoMasksDict, self.inTomoDict = getTsIdsDicts(scoreTomos, tomoMasks,
                                                                                    inTomoSet,
                                                                                    present_ts_ids=presentTsIds)

    def convertInputStep(self, tsId: str):
        logger.info(cyanStr(f'tsId = {tsId}: converting the tomo mask...'))

        try:

            self.make_dirs(tsId)
            scoreTomo = self.scoreTomoDict[tsId]
            mask_choice = self.mask_choice.get()

            # modify output_dir field in json file to save .star file
            jsonIn = scoreTomo.getJsonFile()
            jsonOut = self.getJsonOut(jsonIn, tsId)
            outputDir = self._getCurrentTomoDir(tsId)
            self.copy_and_update_json(jsonIn, jsonOut, outputDir)

            # linking objects to the new output_dir: the extract_coordinates script looks for these objects
            # by reconstructing their own path starting from the output_dir path within the json file

            # linking score tomos
            inScoreTomoFile = scoreTomo.getFileName()
            outScoreTomo = self._getConvertedOrLinkedName(tsId, suffix=SCORE_SUFFIX)
            convertOrLink(inScoreTomoFile, outScoreTomo, samplingRate=scoreTomo.getSamplingRate())

            # linking angles
            inAngles = scoreTomo.getFileName().replace(SCORE_SUFFIX, ANGLES_SUFFIX)
            outAngles = self._getConvertedOrLinkedName(tsId, suffix=ANGLES_SUFFIX)
            convertOrLink(inAngles, outAngles, samplingRate=scoreTomo.getSamplingRate())

            if self.tomoMasksDict:
                if mask_choice == MASK_OTHER:
                    tomomask = self.tomoMasksDict[tsId]
                    msg = check_sr_and_size(tomomask, scoreTomo)
                    if msg:
                        self.failedTsIds.append(tsId)
                        logger.info(yellowStr(f'tsId = {tsId} -> {msg}'))
                        return

                inTomoMaskFile = tomomask.getFileName()
                outTomoMaskFile = self._getConvertedOrLinkedName(tsId, suffix=MASK_SUFFIX)
                convertOrLink(inTomoMaskFile, outTomoMaskFile, samplingRate=scoreTomo.getSamplingRate())

        except Exception as e:
            self.failedTsIds.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def extractCoordinatesStep(self, tsId: str):
        if tsId in self.failedTsIds:
            return
        try:
            logger.info(cyanStr(f'tsId = {tsId}: performing extract coordinates...'))
            Plugin.runPytom(self, self._program, self._generateArguments(tsId), useGpu=False)

        except Exception as e:
            self.failedTsIds.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> pytom extract coordinates failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def createOutputStep(self, tsId: str):
        if tsId in self.failedTsIds:
            return
        try:
            logger.info(cyanStr(f'tsId = {tsId}: generating the output...'))
            outputCoords = self.createOutputSet()
            self.starFile2Coords3D(tsId, outputCoords)

        except Exception as e:
            self.failedTsIds.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> generating the output failed with the exception -> {e}'))
            logger.error(traceback.format_exc())

    def closeOutputSetStep(self):
        self._closeOutputSet()
        outputSet = getattr(self, self._possibleOutputs.coordinates.name, None)

        if outputSet is None or (outputSet is not None and len(outputSet) == 0):
            raise Exception('No coordinates were extracted.')

        # --------------------------- INFO functions ------------------------------

    def _validate(self) -> List[str]:
        valMsg = []
        inTomoMasks = self._getFormAttrib(TOMO_MASKS)
        mask_choice = self.mask_choice.get()

        if mask_choice == MASK_PYTOM_TM and not self.getTMTomoMasks():
            valMsg.append('No TomoMasks were used in the introduced protocol. If you want to use TomoMasks choose '
                          'the "Other TomoMasks" option in parameter "TomoMasks (segmentations)".')
        if mask_choice == MASK_OTHER and not inTomoMasks:
            valMsg.append('No TomoMasks are available but "Other TomoMasks" was selected in parameter '
                          '"TomoMasks (segmentations)".')
        if not self.getScoreTomos():
            valMsg.append('No Score Tomograms were generated by the introduced protocol.')

        return valMsg

    # --------------------------- UTILS functions ------------------------------

    def getScoreTomos(self) -> Optional[SetOfPytomScoreTomograms]:
        protTM = self._getFormAttrib(IN_TM_PROTOCOL)
        scoreTomos = getattr(protTM, protTM._possibleOutputs.scoreTomograms.name, None)
        return scoreTomos

    def getTMTomoMasks(self) -> Optional[SetOfTomoMasks]:
        protTM = self._getFormAttrib(IN_TM_PROTOCOL)
        tomoMasksPointer = getattr(protTM, TOMO_MASKS, String())
        return tomoMasksPointer.get()

    def getTomoMasks(self) -> Optional[SetOfTomoMasks]:
        mask_choice = self.mask_choice.get()

        if mask_choice == MASK_PYTOM_TM:
            return self.getTMTomoMasks()
        if mask_choice == MASK_OTHER:
            return self._getFormAttrib(TOMO_MASKS)
        return None

    def _generateArguments(self, tsId: str) -> str:
        scoreTomo = self.scoreTomoDict[tsId]
        jsonIn = scoreTomo.getJsonFile()
        jsonOut = self.getJsonOut(jsonIn, tsId)
        outTomoMaskFile = self._getConvertedOrLinkedName(tsId, suffix=MASK_SUFFIX)

        cmd = [
            f'--job-file {jsonOut}',
            f'--number-of-particles {self.number_of_particles.get()}',
            f'--number-of-false-positives {self.number_false_positives.get()}',
            f'--particle-diameter {self.particle_diameter.get()}',
            '--relion5-compat',
            '--log info',
            # f'--tophat-bins',
            # f'--plot-bins'
        ]

        if self.tomoMasksDict:
            cmd.append(f'--tomogram-mask {outTomoMaskFile}')

        if self.tophat_filter_con.get() > 0:
            cmd.append('--tophat-filter')
            cmd.append(f'--tophat-connectivity {self.tophat_filter_con.get()}')

        return ' '.join(cmd)

    def starFile2Coords3D(self, tsId: str,
                          coordsSet: SetOfCoordinates3D) -> None:

        starFile = self._getExtraPath(tsId, f'{tsId}_tomo_particles.star')
        dataTable = Table()
        dataTable.read(starFile, tableName='particles')
        inTomo = self.inTomoDict[tsId]

        for row in dataTable:
            # Consider that there can be coordinates in the star file that does not belong to any of the tomograms
            # introduced
            coord = self.gen3dCoordFromStarRow(row, inTomo)
            coordsSet.append(coord)

    def gen3dCoordFromStarRow(self, row, inTomo: Tomogram) -> Coordinate3D:

        coordinate3d = Coordinate3D()
        x = row.get(RLN_CENTEREDCOORDINATEXANGST, 0) / self.sRate
        y = row.get(RLN_CENTEREDCOORDINATEYANGST, 0) / self.sRate
        z = row.get(RLN_CENTEREDCOORDINATEZANGST, 0) / self.sRate
        coordinate3d.setVolume(inTomo)

        coordinate3d.setX(float(x), BOTTOM_LEFT_CORNER)
        coordinate3d.setY(float(y), BOTTOM_LEFT_CORNER)
        coordinate3d.setZ(float(z), BOTTOM_LEFT_CORNER)

        rot = row.get(ROT, 0)
        tilt = row.get(TILT, 0)
        psi = row.get(PSI, 0)

        trMatrix = self.eulerAngles2matrix(rot, tilt, psi)
        coordinate3d.setMatrix(trMatrix)
        score = row.get(SCORE, 0)
        coordinate3d.setScore(score)
        coordinate3d.setBoxSize(self.particle_diameter.get())

        return coordinate3d

    @staticmethod
    def eulerAngles2matrix(tdrot, tilt, narot):
        # Relevant info:
        #   * Pytom's transformation system is ZXZ
        tdrot = np.deg2rad(float(tdrot))
        narot = np.deg2rad(float(narot))
        tilt = np.deg2rad(float(tilt))
        R = transformations.euler_matrix(tdrot, tilt, narot, axes='szxz')

        return R

    def createOutputSet(self) -> SetOfCoordinates3D:
        outCoords = getattr(self, self._possibleOutputs.coordinates.name, None)
        if outCoords:
            outCoords.enableAppend()
        else:
            inScoreTomos = self.getScoreTomos()
            outCoords = SetOfCoordinates3D.create(self._getPath(), template="coordinates%s")
            input_protocol = self._getFormAttrib(IN_TM_PROTOCOL)
            inTomoSet = getattr(input_protocol, IN_TOMOS).get()
            outCoords.setPrecedents(inTomoSet)
            outCoords.setSamplingRate(inScoreTomos.getSamplingRate())
            outCoords.setBoxSize(self.particle_diameter.get())
            outCoords.setStreamState(Set.STREAM_OPEN)

            self._defineOutputs(**{self._possibleOutputs.coordinates.name: outCoords})
            self._defineSourceRelation(inScoreTomos, outCoords)

        return outCoords

    def getJsonOut(self, jsonIn: str, tsId: str) -> str:
        fileName = basename(jsonIn)
        return join(self._getCurrentTomoDir(tsId), fileName)

    @staticmethod
    def copy_and_update_json(source_json: str, dest_json: str, new_output_dir: str) -> None:
        """
        Copies a JSON file and updates the 'output_dir' field.

        :param source_json: path to the original JSON file
        :param dest_json: path where the modified copy will be saved
        :param new_output_dir: new value for the 'output_dir' field
        """
        # 1. Copy the file as-is (optional, in case you want to keep the original untouched)
        shutil.copy(source_json, dest_json)

        # 2. Load the JSON content
        with open(dest_json, "r", encoding="utf-8") as f:
            data = json.load(f)

        # 3. Update the output_dir field
        data["output_dir"] = new_output_dir

        # 4. Write the changes back to the file
        with open(dest_json, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)

TOMO_NAME = 'rlnTomoName'
RLN_CENTEREDCOORDINATEXANGST = 'rlnCenteredCoordinateXAngst'
RLN_CENTEREDCOORDINATEYANGST = 'rlnCenteredCoordinateYAngst'
RLN_CENTEREDCOORDINATEZANGST = 'rlnCenteredCoordinateZAngst'
ROT = 'rlnAngleRot'
TILT = 'rlnAngleTilt'
PSI = 'rlnAnglePsi'
SCORE = 'rlnLCCmax'
